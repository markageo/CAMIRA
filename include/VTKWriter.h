#ifndef VTK_WRITER
#define VTK_WRITER

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>  
#include <functional>
#include <algorithm>
#include <map>
#include <type_traits>


namespace VTK
{

// Settings
#define MODE_ASCII          "ascii"
#define MODE_BINARY         "binary"
#define TYPE_FLOAT          "float"
#define TYPE_DOUBLE         "double"

// Defaults
#define DEFAULT_WRITE_MODE          MODE_ASCII       
#define DEFAULT_ASCII_PRECISION     6


// Data writing mdoes
enum WriteModes {
    ASCII, BINARY
};
std::array<std::string, 2> writeModeStrings{MODE_ASCII, MODE_BINARY};


// Type aliases
using sizeType = int long;
template<typename T> using gridVectorType = std::vector<const T *>;
template<typename T> using scalarMapType = std::map<std::string, const T *>;
template<typename T> using vectorMapType = std::map<std::string, std::vector<const T *>>;


// Utility functions
namespace
{
    // Swap between big endian and little endian representation
    // Credit to Michael Klimenko
    // https://mklimenko.github.io/english/2018/08/22/robust-endian-swap/
    template <typename T>
    T SwapEndian(const T &val) {
        union U {
            T val;
            std::array<std::uint8_t, sizeof(T)> raw;
        };
        U src, dst;
        src.val = val;
        std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
        return dst.val;
    }

}   // end anonymous namespace


/*-------------------------------------------------------------------------------------*\
                                   Configuration Class
\*-------------------------------------------------------------------------------------*/

// Configuration parameters for FieldWriter class
class VTKWriterConfig
{
    public:

        VTKWriterConfig(sizeType, sizeType, sizeType);

        void SetWriteMode(WriteModes);
        void SetASCIIPrecision(int);

        const sizeType &dim(sizeType) const;
        const std::string &WriteMode() const;
        const int &ASCIIPrecision() const;

    private:

        sizeType m_dims[3]; 
        std::string m_writeMode;
        int m_ASCIIPrecision;

};


// ------------------------------------- Class Member Definitions ------------------------------------- //

// Constructor sets default values
VTKWriterConfig::VTKWriterConfig(sizeType dimX,sizeType dimY, sizeType dimZ) :
    m_dims{dimX, dimY, dimZ},
    m_writeMode(DEFAULT_WRITE_MODE), 
    m_ASCIIPrecision(DEFAULT_ASCII_PRECISION)
    {};


// Dimensions
const sizeType &VTKWriterConfig::dim(sizeType i) const
{ return m_dims[i]; }

// Write mode set
void VTKWriterConfig::SetWriteMode(WriteModes writeMode) 
{ m_writeMode = writeModeStrings[ writeMode ]; }

// Write mode get
const std::string &VTKWriterConfig::WriteMode() const
{ return m_writeMode; }

// ASCII precision set
void VTKWriterConfig::SetASCIIPrecision(int ASCIIPrecision) 
{ m_ASCIIPrecision = ASCIIPrecision; }

// ASCII precision get
const int &VTKWriterConfig::ASCIIPrecision() const
{ return m_ASCIIPrecision; }


/*-------------------------------------------------------------------------------------*\
                                      Writer Class
\*-------------------------------------------------------------------------------------*/

// Write vector and scalar fields to file in legacy .vtk format
template<typename T>
class VTKWriter
{
    static_assert(std::is_same<T, float>::value ||
                  std::is_same<T, double>::value,
                  "Data must be either double or float type!");

    public:
        VTKWriter(const gridVectorType<T> &, const scalarMapType<T> &, const vectorMapType<T> &, const VTKWriterConfig &);
        VTKWriter(const gridVectorType<T> &, const scalarMapType<T> &, const VTKWriterConfig &);
        VTKWriter(const gridVectorType<T> &, const vectorMapType<T> &, const VTKWriterConfig &);
        VTKWriter(const gridVectorType<T> &, const VTKWriterConfig &);

        int WriteData(const std::string &, const std::string &);

    private:

        VTKWriterConfig m_config;
        std::ofstream m_outputFileStream;
        std::string m_outputDataType;
        std::string m_fileExtension = ".vtk"; 

        gridVectorType<T> m_gridVector;
        scalarMapType<T> m_scalarMap;
        vectorMapType<T> m_vectorMap;

        std::function<void(T)> DataProperWriteFunc;
        void WriteDataArray(const T *, const sizeType &) const ;
        void WriteDataArray(const std::vector<const T *> , const sizeType &) const;
        void WriteDatasetRectilinearGrid();
        void WriteDataAttributeScalar(const std::string &, const T *);
        void WriteDataAttributeVector(const std::string &, const std::vector<const T *> &);
};


// ------------------------------------- Class Member Definitions ------------------------------------- //

// Constructors
template<typename T>
VTKWriter<T>::VTKWriter(const gridVectorType<T> &gridVector, 
                        const scalarMapType<T> &scalarMap, 
                        const vectorMapType<T> &vectorMap, 
                        const VTKWriterConfig &config) : 
    m_config(config),
    m_gridVector(gridVector), 
    m_scalarMap(scalarMap),
    m_vectorMap(vectorMap) 
    {
        // Set the writing function
        if (config.WriteMode() == MODE_ASCII) {
            DataProperWriteFunc =  [&](T data) { 
                m_outputFileStream << data << " "; 
            };
            m_outputFileStream << std::fixed << std::setprecision(config.ASCIIPrecision());

        } else if (config.WriteMode() == MODE_BINARY) {
            DataProperWriteFunc = [&](T data) {
                T tmp = SwapEndian(data);  // Legacy VTK uses big endian
                m_outputFileStream.write(reinterpret_cast<const char *>(&tmp), sizeof(tmp)); 
            };
        } 

        // Set the datatype
        if constexpr        ( std::is_same<T, double>::value ) {
            m_outputDataType = TYPE_DOUBLE;

        } else if constexpr ( std::is_same<T, float>::value ) {
            m_outputDataType = TYPE_FLOAT;
            
        }
    }


template<typename T>
VTKWriter<T>::VTKWriter(const gridVectorType<T> &gridVector, 
                        const scalarMapType<T> &scalarMap, 
                        const VTKWriterConfig &config) : 
    VTKWriter(gridVector, scalarMap, {}, config) {}


template<typename T>
VTKWriter<T>::VTKWriter(const gridVectorType<T> &gridVector, 
                        const vectorMapType<T> &vectorMap, 
                        const VTKWriterConfig &config) : 
    VTKWriter(gridVector, {}, vectorMap, config) {}


template<typename T>
VTKWriter<T>::VTKWriter(const gridVectorType<T> &gridVector, 
                        const VTKWriterConfig &config) : 
    VTKWriter(gridVector, {}, {}, config) {}



// Writes fields to file
// return 0 for success
// return -1 for failure
template<typename T>
int VTKWriter<T>::WriteData(const std::string &filename, 
                            const std::string &title)
{
    // Add .vtk file extension to file if it is not already there
    std::string ext = "";
    if (filename.substr(filename.find_last_of(".")) != m_fileExtension) {
        ext = m_fileExtension;
    }

    // Open new output file stream
    if (m_outputFileStream.is_open()) {
        m_outputFileStream.close();
    }
    m_outputFileStream.open(filename+ext, std::ofstream::binary);

    if (!m_outputFileStream) {
        std::cout << "ERROR: Could not open file." << std::endl;
        return -1;
    }

    // Header
    m_outputFileStream << "# vtk DataFile Version 3.0" << "\n";

    // Title
    m_outputFileStream << title << "\n";

    // Data type
    m_outputFileStream << m_config.WriteMode() << "\n";

    // Dataset
    WriteDatasetRectilinearGrid();
    m_outputFileStream << "\n";     

    // Dataset attributes
    m_outputFileStream << "POINT_DATA"
                       << " " << m_config.dim(0)*m_config.dim(1)*m_config.dim(2)
                       << "\n";
                    
    // Scalars
    for (const auto& [scalarFieldName, scalarField] : m_scalarMap) {
        WriteDataAttributeScalar(scalarFieldName, scalarField);
        m_outputFileStream << "\n";
    }

    // Vectors
    for (const auto& [vectorFieldName, vectorField] : m_vectorMap) {
        WriteDataAttributeVector(vectorFieldName, vectorField);
        m_outputFileStream << "\n";
    }

    m_outputFileStream.flush();
    m_outputFileStream.close();
    return 0;
}



// Iterate raw data pointer and write 1D data
template<typename T>
void VTKWriter<T>::WriteDataArray(const T *dataPtr, 
                                  const sizeType &iterLength) const {
    for (sizeType i = 0; i != iterLength; i++) {
        DataProperWriteFunc( dataPtr[i] );
    }
}

// Iteratre raw data pointer and write 3D data
template<typename T>
void VTKWriter<T>::WriteDataArray(const std::vector<const T *> dataPtrVec, 
                                  const sizeType &iterLength) const {
    auto dataPtr1 = dataPtrVec[0];
    auto dataPtr2 = dataPtrVec[1];
    auto dataPtr3 = dataPtrVec[2];
    for (sizeType i = 0; i != iterLength; i++) {
        DataProperWriteFunc( dataPtr1[i] );
        DataProperWriteFunc( dataPtr2[i] );
        DataProperWriteFunc( dataPtr3[i] );
    }
}


// Write rectilinear grid data to file
template<typename T>
void VTKWriter<T>::WriteDatasetRectilinearGrid()
{

    m_outputFileStream << "DATASET RECTILINEAR_GRID" << "\n";
    m_outputFileStream << "DIMENSIONS"
                       << " " << m_config.dim(0)
                       << " " << m_config.dim(1)
                       << " " << m_config.dim(2)
                       << "\n";

    m_outputFileStream << "X_COORDINATES"
                       << " " << m_config.dim(0)
                       << " " << m_outputDataType
                       << "\n";
    WriteDataArray(m_gridVector[0], m_config.dim(0));
    m_outputFileStream << "\n";
                       
    m_outputFileStream << "Y_COORDINATES"
                       << " " << m_config.dim(1)
                       << " " << m_outputDataType
                       << "\n";
    WriteDataArray(m_gridVector[1], m_config.dim(1));
    m_outputFileStream << "\n";

    m_outputFileStream << "Z_COORDINATES"
                       << " " << m_config.dim(2)
                       << " " << m_outputDataType
                       << "\n";
    WriteDataArray(m_gridVector[2], m_config.dim(2));
    m_outputFileStream << "\n";
}


// Write scalar data attribute to file
template<typename T>
void VTKWriter<T>::WriteDataAttributeScalar(const std::string &scalarFieldName, 
                                            const T *pScalarField)
{
    m_outputFileStream << "SCALARS" 
                           << " " << scalarFieldName
                           << " " << m_outputDataType
                           << " " << "1"
                           << "\n";

    m_outputFileStream << "LOOKUP_TABLE default" << "\n";
    WriteDataArray(pScalarField, m_config.dim(0)*m_config.dim(1)*m_config.dim(2));
}


// Write vector data attribute to file
template<typename T>
void VTKWriter<T>::WriteDataAttributeVector(const std::string &vectorFieldName, 
                                            const std::vector<const T *> &vectorField)
{
    m_outputFileStream << "VECTORS" 
                           << " " << vectorFieldName
                           << " " << m_outputDataType
                           << "\n";
    WriteDataArray(vectorField, m_config.dim(0)*m_config.dim(1)*m_config.dim(2));
}


} // end namespace VTK


#endif // VTK_WRITER
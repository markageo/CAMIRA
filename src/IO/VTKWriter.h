#ifndef VTK_WRITER
#define VTK_WRITER

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>  
#include <functional>
#include <algorithm>
#include <type_traits>


namespace VTK
{
    
#define DEFAULT_ASCII_PRECISION     6   // Maybe name this something else since it is now in a header

// Data writing mdoes
enum class WriteModes {
    ASCII, BINARY
};

// Point or cell data
enum class GridTypes {
    CELL_DATA, POINT_DATA
};


template< typename T >
struct ScalarData {
    std::string fieldName;
    GridTypes gridType;
    const T * dataPointer; 
};

template< typename T >
struct VectorData {
    std::string fieldName;
    GridTypes gridType;
    std::vector<const T*> dataPointers; 
};

// Type aliases
using sizeType = long;
template<typename T> using gridVectorType = std::vector<const T *>;
template<typename T> using scalarMapType = std::vector< ScalarData<T> >;
template<typename T> using vectorMapType = std::vector< VectorData<T> >;


// Utility functions
namespace INTERNAL
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

}   // end namespace INTERNAL





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

        const std::array<sizeType, 3> &dims() const;
        sizeType dim(sizeType) const;
        WriteModes WriteMode() const;
        const int &ASCIIPrecision() const;

    private:

        std::array<sizeType, 3> m_dims; 
        VTK::WriteModes m_writeMode;
        int m_ASCIIPrecision;

};


// ------------------------------------- Class Member Definitions ------------------------------------- //

// Constructor sets default values
VTKWriterConfig::VTKWriterConfig(sizeType dimX, sizeType dimY, sizeType dimZ) :
    m_dims{dimX, dimY, dimZ},
    m_writeMode( WriteModes::ASCII ), 
    m_ASCIIPrecision( DEFAULT_ASCII_PRECISION )
    {};


// Dimensions
sizeType VTKWriterConfig::dim(sizeType i) const
{ return m_dims[ static_cast<size_t>(i) ]; }

const std::array<sizeType, 3> &VTKWriterConfig::dims() const
{ return m_dims; } 


// Write mode set
void VTKWriterConfig::SetWriteMode(WriteModes writeMode) 
{ m_writeMode = writeMode; }

// Write mode get
WriteModes VTKWriterConfig::WriteMode() const
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
        std::string m_writeMode;
        std::string m_outputDataType;
        std::string m_gridType;
        std::string m_fileExtension = ".vtk"; 

        std::array<sizeType, 3> m_gridDims;
        sizeType m_nPointData, m_nCellData;

        gridVectorType<T> m_gridVector;
        scalarMapType<T> m_scalarMapPointData, m_scalarMapCellData;
        vectorMapType<T> m_vectorMapPointData, m_vectorMapCellData;

        std::function<void(T)> DataProperWriteFunc;
        void WriteDataArray(const T *, const sizeType &) const ;
        void WriteDataArray(const std::vector<const T *> , const sizeType &) const;
        void WriteDatasetRectilinearGrid();
        void WriteDataAttributeScalar(const std::string &, const T *, const sizeType);
        void WriteDataAttributeVector(const std::string &, const std::vector<const T *> &, const sizeType);
        void WriteFieldData( const scalarMapType<T>, const vectorMapType<T>, const sizeType);
};


// ------------------------------------- Class Member Definitions ------------------------------------- //

// Constructors
template<typename T>
VTKWriter<T>::VTKWriter(const gridVectorType<T> &gridVector, 
                        const scalarMapType<T> &scalarMap, 
                        const vectorMapType<T> &vectorMap, 
                        const VTKWriterConfig &config) : 
    m_config(config),
    m_gridDims( config.dims() ),
    m_nPointData( m_gridDims[0] * m_gridDims[1] * m_gridDims[2] ),
    m_nCellData( (m_gridDims[0]-1) * (m_gridDims[1]-1) * (m_gridDims[2]-1) ),
    m_gridVector( gridVector )
    {
        using enum WriteModes;
        using enum GridTypes;

        switch ( config.WriteMode() ) {
            case ASCII:
                DataProperWriteFunc =  [&](T data) { 
                    m_outputFileStream << data << " "; 
                };
                m_outputFileStream << std::fixed << std::setprecision(config.ASCIIPrecision());
                m_writeMode = "ASCII";
                break;

            case BINARY:
                DataProperWriteFunc = [&](T data) {
                    T tmp = INTERNAL::SwapEndian(data);  // Legacy VTK uses big endian
                    m_outputFileStream.write(reinterpret_cast<const char *>(&tmp), sizeof(tmp)); 
                };
                m_writeMode = "BINARY";
                break;
        }


        if constexpr        ( std::is_same<T, double>::value ) {
            m_outputDataType = "DOUBLE";
        } else if constexpr ( std::is_same<T, float>::value ) {
            m_outputDataType = "FLOAT";
        }


        for ( const auto &scalarData : scalarMap ) {
            switch ( scalarData.gridType ) {
                case GridTypes::POINT_DATA:
                    m_scalarMapPointData.push_back( scalarData );
                    break;

                case GridTypes::CELL_DATA:  
                    m_scalarMapCellData.push_back( scalarData );
                    break;
            }
        }


        for ( const auto &vectorData : vectorMap ) {
            switch ( vectorData.gridType ) {
                case GridTypes::POINT_DATA:
                    m_vectorMapPointData.push_back( vectorData );
                    break;

                case GridTypes::CELL_DATA:  
                    m_vectorMapCellData.push_back( vectorData );
                    break;
            }
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

    if ( !m_outputFileStream ) {
        std::cout << "VTK Writer ERROR: Could not open/create file. Output will not be written.\n" << std::endl;
        return -1;
    }

    // Header
    m_outputFileStream << "# vtk DataFile Version 3.0" << "\n";

    // Title
    m_outputFileStream << title << "\n";

    // Data type
    m_outputFileStream << m_writeMode << "\n";

    // Dataset
    WriteDatasetRectilinearGrid();
    m_outputFileStream << "\n";     

    // Point data
    m_outputFileStream << "POINT_DATA"
                       << " " << m_nPointData
                       << "\n";
    WriteFieldData(m_scalarMapPointData, m_vectorMapPointData, m_nPointData);

    // Cell data
    m_outputFileStream << "CELL_DATA"
                       << " " << m_nCellData
                       << "\n";
    WriteFieldData(m_scalarMapCellData, m_vectorMapCellData, m_nCellData);

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
    auto WriteAxis = [&] ( const std::string &axisName, const size_t dim ) {
        m_outputFileStream << axisName
                           << " " << m_gridDims[ dim ]
                           << " " << m_outputDataType
                           << "\n";
        WriteDataArray(m_gridVector[ dim ], m_gridDims[ dim ]);
        m_outputFileStream << "\n";
    };

    m_outputFileStream << "DATASET RECTILINEAR_GRID" << "\n";
    m_outputFileStream << "DIMENSIONS"
                       << " " << m_gridDims[0]
                       << " " << m_gridDims[1]
                       << " " << m_gridDims[2]
                       << "\n";
    WriteAxis( "X_COORDINATES", 0 );
    WriteAxis( "Y_COORDINATES", 1 );
    WriteAxis( "Z_COORDINATES", 2 );
}


// Write scalar data attribute to file
template<typename T>
void VTKWriter<T>::WriteDataAttributeScalar(const std::string &scalarFieldName, 
                                            const T *pScalarField,
                                            const sizeType nPoints)
{
    m_outputFileStream << "SCALARS" 
                           << " " << scalarFieldName
                           << " " << m_outputDataType
                           << " " << "1"
                           << "\n";
    m_outputFileStream << "LOOKUP_TABLE default" << "\n";

    WriteDataArray( pScalarField, nPoints );
}


// Write vector data attribute to file
template<typename T>
void VTKWriter<T>::WriteDataAttributeVector(const std::string &vectorFieldName, 
                                            const std::vector<const T *> &vectorField,
                                            const sizeType nPoints)
{
    m_outputFileStream << "VECTORS" 
                           << " " << vectorFieldName
                           << " " << m_outputDataType
                           << "\n";

    WriteDataArray(vectorField, nPoints );
}


template<typename T>
void VTKWriter<T>::WriteFieldData( const scalarMapType<T> scalarMap,
                     const vectorMapType<T> vectorMap,
                     const sizeType nPoints)
{
    // Scalars
    for (const auto& scalarData : scalarMap) {
        WriteDataAttributeScalar(scalarData.fieldName, scalarData.dataPointer, nPoints);
        m_outputFileStream << "\n\n";
    }

    // Vectors
    for (const auto& vectorData : vectorMap) {
        WriteDataAttributeVector(vectorData.fieldName, vectorData.dataPointers, nPoints);
        m_outputFileStream << "\n\n";
    }
}
    


} // end namespace VTK


#endif // VTK_WRITER
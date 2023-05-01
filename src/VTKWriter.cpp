#include "VTKWriter.h"

// Settings
#define MODE_ASCII          "ascii"
#define MODE_BINARY         "binary"
#define TYPE_FLOAT          "float"
#define TYPE_DOUBLE         "double"

// Defaults
#define DEFAULT_WRITE_MODE          MODE_ASCII       
#define DEFAULT_ASCII_PRECISION     6

using namespace VTK;

// Helper functions
namespace
{
    // Return input string with all characters in lower case
    std::string StringLower(const std::string &str) 
    {
        std::string strLower;
        for (std::string::const_iterator istr = str.begin(); istr != str.end(); ++istr) {
            strLower += std::tolower(*istr);
        }
        return strLower;
    }
}


/*-------------------------------------------------------------------------------------*\
                                    Constructors
\*-------------------------------------------------------------------------------------*/

VTKWriter::VTKWriter(const gridVectorType &gridVector, const scalarMapType &scalarMap, const vectorMapType &vectorMap, const VTKWriterConfig &config) : 
    m_config(config),
    m_gridVector(gridVector), 
    m_scalarMap(scalarMap),
    m_vectorMap(vectorMap) 
    {
        if (config.WriteMode() == MODE_ASCII) {
            DataProperWriteFunc =  [&](double data) {   // There should be a float version !?!?!?!?!?
                m_outputFileStream << data << " "; 
            };
            m_outputFileStream << std::fixed << std::setprecision(config.ASCIIPrecision());
            m_outputDataType = TYPE_DOUBLE; 

        } else if (config.WriteMode() == MODE_BINARY) {

            if (config.DataType() == DOUBLE) {
                DataProperWriteFunc = [&](double data) {
                    double tmp = SwapEndian(data);  // Legacy VTK uses big endian
                    m_outputFileStream.write(reinterpret_cast<const char *>(&tmp), sizeof(tmp)); 
                };
                m_outputDataType = TYPE_DOUBLE;

            } else if (config.DataType() == FLOAT) {
                DataProperWriteFunc = [&](float data) {
                    float tmp = SwapEndian(data);  // Legacy VTK uses big endian
                    m_outputFileStream.write(reinterpret_cast<const char *>(&tmp), sizeof(tmp)); 
                };  
                m_outputDataType = TYPE_FLOAT;

            }
        } 
    }

VTKWriter::VTKWriter(const gridVectorType &gridVector, const scalarMapType &scalarMap, const VTKWriterConfig &config) : 
    VTKWriter(gridVector, scalarMap, {}, config) {}

VTKWriter::VTKWriter(const gridVectorType &gridVector, const vectorMapType &vectorMap, const VTKWriterConfig &config) : 
    VTKWriter(gridVector, {}, vectorMap, config) {}


/*-------------------------------------------------------------------------------------*\
                               Public Member Functions
\*-------------------------------------------------------------------------------------*/


// Writes fields to file
// return 0 for success
// return -1 for failure
int VTKWriter::WriteData(const std::string &filename, const std::string &title)
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
        std::cout << "ERROR: Cout not open file." << std::endl;
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

/*-------------------------------------------------------------------------------------*\
                                Private Member Functions
\*-------------------------------------------------------------------------------------*/


// Iterate raw data pointer and write 1D data
void VTKWriter::WriteDataArray(const void *voidDataPtr, const sizeType &iterLength) const {
    if (m_config.DataType() == DOUBLE) {
        auto dataPtr = reinterpret_cast<const double *>(voidDataPtr);
        for (sizeType i = 0; i != iterLength; i++) {
            DataProperWriteFunc( dataPtr[i] );
        }
    } else if (m_config.DataType() == FLOAT) {
        auto dataPtr = reinterpret_cast<const float *>(voidDataPtr);
        for (sizeType i = 0; i != iterLength; i++) {
            DataProperWriteFunc( dataPtr[i] );
        }
    }
}

// Iteratre raw data pointer and write 3D data
void VTKWriter::WriteDataArray(const std::vector<const void *> voidDataPtrVec, const sizeType &iterLength) const {
    if (m_config.DataType() == DOUBLE) {
        auto dataPtr1 = reinterpret_cast<const double *>(voidDataPtrVec[0]);
        auto dataPtr2 = reinterpret_cast<const double *>(voidDataPtrVec[1]);
        auto dataPtr3 = reinterpret_cast<const double *>(voidDataPtrVec[2]);
        for (sizeType i = 0; i != iterLength; i++) {
            DataProperWriteFunc( dataPtr1[i] );
            DataProperWriteFunc( dataPtr2[i] );
            DataProperWriteFunc( dataPtr3[i] );
        }
    } else if (m_config.DataType() == FLOAT) {
        auto dataPtr1 = reinterpret_cast<const float *>(voidDataPtrVec[0]);
        auto dataPtr2 = reinterpret_cast<const float *>(voidDataPtrVec[1]);
        auto dataPtr3 = reinterpret_cast<const float *>(voidDataPtrVec[2]);
        for (sizeType i = 0; i != iterLength; i++) {
            DataProperWriteFunc( dataPtr1[i] );
            DataProperWriteFunc( dataPtr2[i] );
            DataProperWriteFunc( dataPtr3[i] );
        }
    }
}


// Write rectilinear grid data to file
void VTKWriter::WriteDatasetRectilinearGrid()
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
void VTKWriter::WriteDataAttributeScalar(const std::string &scalarFieldName, const void *pScalarField)
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
void VTKWriter::WriteDataAttributeVector(const std::string &vectorFieldName, const std::vector<const void *> &vectorField)
{
    m_outputFileStream << "VECTORS" 
                           << " " << vectorFieldName
                           << " " << m_outputDataType
                           << "\n";
    WriteDataArray(vectorField, m_config.dim(0)*m_config.dim(1)*m_config.dim(2));
}


// Swap between big endian and little endian representation
// Credit to: Michael Klimenko
// https://mklimenko.github.io/english/2018/08/22/robust-endian-swap/
template <typename T>
T VTKWriter::SwapEndian(const T &val) {
    union U {
        T val;
        std::array<std::uint8_t, sizeof(T)> raw;
    };
    U src, dst;
    src.val = val;
    std::reverse_copy(src.raw.begin(), src.raw.end(), dst.raw.begin());
    return dst.val;
}



/*-------------------------------------------------------------------------------------*\
                                Configuration Class
\*-------------------------------------------------------------------------------------*/

// Constructor sets default values
VTKWriterConfig::VTKWriterConfig(const sizeType &dimX, const sizeType &dimY, const sizeType &dimZ, dataType type) :
    m_dims{dimX, dimY, dimZ},
    m_writeMode(DEFAULT_WRITE_MODE), 
    m_ASCIIPrecision(DEFAULT_ASCII_PRECISION),
    m_dataType(type) 
    {};

// Data type get
const dataType &VTKWriterConfig::DataType() const
{ return m_dataType; }

// Dimensions
const sizeType &VTKWriterConfig::dim(const sizeType &i) const
{ return m_dims[i]; }

// Write mode set
void VTKWriterConfig::SetWriteMode(const std::string &writeMode) 
{
    std::string writeModeLower = StringLower(writeMode);
    if (writeModeLower != MODE_ASCII && writeModeLower != MODE_BINARY) {
        std::cout << "WARNING: Invalid write mode given, writing to ASCII..." << std::endl;
        m_writeMode = MODE_ASCII;
    } else {
        m_writeMode = writeModeLower;
    }
}

// Write mode get
const std::string &VTKWriterConfig::WriteMode() const
{ return m_writeMode; }

// ASCII precision set
void VTKWriterConfig::SetASCIIPrecision(const unsigned &ASCIIPrecision) 
{ m_ASCIIPrecision = ASCIIPrecision; }

// ASCII precision get
const int &VTKWriterConfig::ASCIIPrecision() const
{ return m_ASCIIPrecision; }

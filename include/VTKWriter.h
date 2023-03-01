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

using sizeType = int long;
using gridVectorType = std::vector<void *> ;
using scalarMapType = std::map<std::string, void *>;
using vectorMapType = std::map<std::string, std::vector<void *>>;

typedef enum {
    DOUBLE, FLOAT
} dataType;


// Configuration parameters for FieldWriter class
class VTKWriterConfig
{
    public:

        VTKWriterConfig(const sizeType &, const sizeType &, const sizeType &, dataType);

        void SetWriteMode(const std::string &);
        void SetASCIIPrecision(const unsigned &);

        const sizeType &dim(const sizeType &) const;
        const dataType &DataType() const;
        const std::string &WriteMode() const;
        const int &ASCIIPrecision() const;

    private:

        sizeType m_dims[3]; 
        std::string m_writeMode;
        int m_ASCIIPrecision;
        dataType m_dataType;
        

};


// Write vector and scalar fields to file in legacy .vtk format
class VTKWriter
{
    public:

        VTKWriter(const gridVectorType &, const scalarMapType &, const vectorMapType &, const VTKWriterConfig &);
        VTKWriter(const gridVectorType &, const scalarMapType &, const VTKWriterConfig &);
        VTKWriter(const gridVectorType &, const vectorMapType &, const VTKWriterConfig &);


        int WriteData(const std::string &, const std::string &);

    private:

        VTKWriterConfig m_config;
        std::ofstream m_outputFileStream;
        std::string m_outputDataType;
        std::string m_fileExtension = ".vtk"; 

        gridVectorType m_gridVector;
        scalarMapType m_scalarMap;
        vectorMapType m_vectorMap;

        std::function<void(double)> DataProperWriteFunc;
        void WriteDataArray(const void *, const sizeType &) const ;
        void WriteDataArray(const std::vector<void *> , const sizeType &) const;
        void WriteDatasetRectilinearGrid();
        void WriteDataAttributeScalar(const std::string &, const void *);
        void WriteDataAttributeVector(const std::string &, const std::vector<void *> &);
        template <typename T> T SwapEndian(const T &);

};


} // end namespace VTKwriter

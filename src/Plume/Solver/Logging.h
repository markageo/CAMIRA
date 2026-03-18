#ifndef CAMIRA_PLUME_LOGGING
#define CAMIRA_PLUME_LOGGING

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Macros.h"
#include "Core/IO/VTKWriter.h"
#include "Core/IO/IOTools.h"
#include "Plume/InputProcessing/InputProcessing.h"


#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

// Set the output precision based on code precision
#ifdef CAMIRA_DOUBLE_PRECISION
    #define CAMIRA_FILE_WRITE_PRECISION 15
#else
    #define CAMIRA_FILE_WRITE_PRECISION 7;
#endif



// Postprocesses the solution fields and writes to file
class ConcentrationFieldWriter
{
    public:
        ConcentrationFieldWriter( const Tensor3D &concentrationField, 
                                  const Mesh &mesh,
                                  const InputData &inputData ) :
            m_concentrationField( concentrationField ),
            m_mesh( mesh ),
            m_baseFilename( IOTOOLS::RemoveFileExtension( inputData.fieldOutputFilename, ".vtk" ) )
            {

                switch ( inputData.outputFormatType ) {
                    case InputData::OutputFormatType::BINARY:
                        m_writeMode = VTK::WriteModes::BINARY;
                        break;

                    case InputData::OutputFormatType::ASCII:
                        m_writeMode = VTK::WriteModes::ASCII;
                        break;
                }

            };

            void WriteData( floatType currentTime, 
                            intType timeStep )
            { 
                SetWriter();

                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << currentTime;
                const std::string currentTimeString = oss.str();

                const std::string message  = "CFD solution at time " + currentTimeString + ", timestep = " + std::to_string(timeStep);
                const std::string filename = m_baseFilename  + "_timeStep" + std::to_string(timeStep) + ".vtk";
                m_vtkWriter->WriteData( filename, message ); 
            }


    private:
        const Tensor3D &m_concentrationField;
        const Mesh &m_mesh;
        std::unique_ptr< VTK::VTKWriter<floatType> > m_vtkWriter;
        const std::string m_baseFilename;
        VTK::WriteModes m_writeMode;

        void SetWriter()
        {
            using enum Axis::ENUMDATA;
            VTK::VTKWriterConfig config( m_mesh.cellFaces[X].size(), 
                                         m_mesh.cellFaces[Y].size(), 
                                         m_mesh.cellFaces[Z].size() );
                config.SetWriteMode( m_writeMode );
                
            VTK::gridVectorType<CAMIRA::floatType> gridVector = { m_mesh.cellFaces[X].data(), 
                                                                  m_mesh.cellFaces[Y].data(), 
                                                                  m_mesh.cellFaces[Z].data() };

            VTK::scalarCollectionType<floatType> scalarMap = { {"Concentration" , VTK::GridTypes::CELL_DATA , m_concentrationField.data()} };

            VTK::vectorCollectionType<floatType> vectorMap = { };
        
            m_vtkWriter = std::make_unique<VTK::VTKWriter<floatType>>(gridVector, scalarMap, vectorMap, config);
        }


};



}   // end namespace PLUME

}   // end namespace CAMIRA


#endif // CAMIRA_PLUME_LOGGING
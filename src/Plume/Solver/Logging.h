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
        ConcentrationFieldWriter( const Mesh &mesh,
                                  const InputData &inputData ) :
            m_mesh( mesh ),
            m_instantaneousFieldBaseFilename( IOTOOLS::RemoveFileExtension( inputData.fieldOutputFilename, ".vtk" ) )
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

            void WriteInstantaneousField( const Tensor3D &field,
                                          floatType currentTime, 
                                          intType timeStep )
            { 
                SetWriter( field );

                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << currentTime;
                const std::string currentTimeString = oss.str();

                const std::string message  = "Concentration field at time " + currentTimeString + ", timestep = " + std::to_string(timeStep);
                const std::string filename = m_instantaneousFieldBaseFilename  + "_timeStep" + std::to_string(timeStep) + ".vtk";
                m_vtkWriter->WriteData( filename, message ); 
            }


            void WriteTimeAveragedField( const Tensor3D &field,
                                        const std::string &filename,
                                        const intType startTimeStep,
                                        const intType endTimeStep,
                                        const floatType timeStepSize )
            {
                SetWriter( field );

                const floatType startTime = static_cast<floatType>(startTimeStep) * timeStepSize,
                                endTime   = static_cast<floatType>(endTimeStep)   * timeStepSize;

                std::ostringstream oss;
                oss << std::fixed << std::setprecision(3) << startTime;
                const std::string startTimeString = oss.str();

                oss.str();
                oss << std::fixed << std::setprecision(3) << endTime;
                const std::string endTimeString = oss.str();

                const std::string message  = "Time averaged concentration field between times " 
                                           + startTimeString + "(" + std::to_string(startTimeStep) + ") and " 
                                           + endTimeString   + "(" + std::to_string(endTimeStep)   + ")";
                m_vtkWriter->WriteData( filename, message ); 
            }


    private:
        const Mesh &m_mesh;
        std::unique_ptr< VTK::VTKWriter<floatType> > m_vtkWriter;
        const std::string m_instantaneousFieldBaseFilename;
        VTK::WriteModes m_writeMode;

        void SetWriter( const Tensor3D &field )
        {
            using enum Axis::ENUMDATA;
            VTK::VTKWriterConfig config( m_mesh.cellFaces[X].size(), 
                                         m_mesh.cellFaces[Y].size(), 
                                         m_mesh.cellFaces[Z].size() );
                config.SetWriteMode( m_writeMode );
                
            VTK::gridVectorType<CAMIRA::floatType> gridVector = { m_mesh.cellFaces[X].data(), 
                                                                  m_mesh.cellFaces[Y].data(), 
                                                                  m_mesh.cellFaces[Z].data() };

            VTK::scalarCollectionType<floatType> scalarMap = { {"Concentration" , VTK::GridTypes::CELL_DATA , field.data()} };

            VTK::vectorCollectionType<floatType> vectorMap = { };
        
            m_vtkWriter = std::make_unique<VTK::VTKWriter<floatType>>(gridVector, scalarMap, vectorMap, config);
        }


};



}   // end namespace PLUME

}   // end namespace CAMIRA


#endif // CAMIRA_PLUME_LOGGING
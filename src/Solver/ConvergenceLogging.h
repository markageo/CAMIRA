#ifndef SOLVER_LOGGING
#define SOLVER_LOGGING

#include "../Types.h"
#include "../Tools/SweepTransformations.h"
#include "../Tools/FVTools.h"
#include "../IO/IOTools.h"

#include "Solver.h"
#include "../Tools/FieldProbe.h"
#include "../IO/VTKWriter.h"
#include "../Macros.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <memory>

namespace CFD
{

class ConvergenceFile
{

    public:
        ConvergenceFile( const std::string &filename, 
                         const int precision = 6, 
                         const int columnWidth = 14 ) :
            m_fileStream( filename ),
            m_precision( precision ),
            m_columnWidth( columnWidth )
        {           
            if ( columnWidth < precision + 8 )
                m_columnWidth = precision + 8; 

            if ( !m_fileStream ) {
                std::cout << "Convergence Log File ERROR: Could not open/create file! Convergence history will not be written. \n" << std::endl;
    }

            m_fileStream << std::setprecision( m_precision ) << std::scientific;
        }

        template <class ...Args> 
        void WriteLine( Args... args )
        {
            const int nArgs = sizeof...(args);
            int i = 0;
        
            // Fold expression to loop through args. This is needed since each arg can have a different type
            ( [&] {

                i++;
                m_fileStream << std::left << std::setw( m_columnWidth ) << args;
                if ( i != nArgs )
                    m_fileStream << ", ";

            } (), ... );

            m_fileStream << std::endl;  // flush to write immediately to file
        }

        void WriteCommentLine( const std::string &comment )
        {
            m_fileStream << "# " << comment << "\n";
        }

    private:
        std::ofstream m_fileStream;
        int m_precision;
        int m_columnWidth;
};





class ResidualLogFile
{

    public:
        ResidualLogFile( const std::string &filename, 
                         const AxisTransformationMap &axisTransformation,
                         const int precision = 6 ) :
            m_AT( axisTransformation ),
            m_convergenceFile( filename, precision )
        {            
            // Comment description of file
            m_convergenceFile.WriteCommentLine( "Residuals convergence history" );

            // Write header
            m_convergenceFile.WriteLine( "Iteration", 
                                          "U_residual", 
                                          "V_residual",
                                          "W_residual",
                                          "P_residual", 
                                          "Global_mass_residual");
            
        }

        void WriteData( const FieldData<floatType> &residuals, 
                        const floatType massFluxResidual, 
                        const intType nIterations )
        {
            using enum Axis::ENUMDATA;

            m_convergenceFile.WriteLine( nIterations,
                                         residuals.U[ m_AT.CodeAxis(X) ], 
                                         residuals.U[ m_AT.CodeAxis(Y) ],
                                         residuals.U[ m_AT.CodeAxis(Z) ], 
                                         residuals.P                    , 
                                         massFluxResidual );
        }


    private:
        AxisTransformationMap m_AT;
        ConvergenceFile m_convergenceFile;
};





class ProbeLogFile
{

    public:
        ProbeLogFile( const std::string &filename, 
                      const AxisTransformationMap &axisTransformation,
                      const FieldProbe &fieldProbe,
                      const int precision = 6 ) :
            m_AT( axisTransformation ),
            m_convergenceFile( filename, precision )
        {            
            using enum Axis::ENUMDATA;

            // Comment description of file
            std::string fileDescription = "Field probe at coordinates: (" 
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(X) ) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(Y) ) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(Z) ) ) + ").";
            m_convergenceFile.WriteCommentLine( fileDescription );

            // Write header
            m_convergenceFile.WriteLine( "Iteration", 
                                         "X_velocity", 
                                         "Y_velocity",
                                         "Z_velocity",
                                         "Pressure");
            
        }

        void WriteData( const FieldData<floatType> &probeValues,  
                        const intType nIterations )
        {
            using enum Axis::ENUMDATA;

            m_convergenceFile.WriteLine( nIterations,
                                         probeValues.U[ m_AT.CodeAxis(X) ], 
                                         probeValues.U[ m_AT.CodeAxis(Y) ],
                                         probeValues.U[ m_AT.CodeAxis(Z) ], 
                                         probeValues.P );
        }


    private:
        AxisTransformationMap m_AT;
        ConvergenceFile m_convergenceFile;
};





class ConsoleLog
{

    public:
        ConsoleLog( const AxisTransformationMap &axisTransformation,
                    const int precision = 6 ) :
            m_AT( axisTransformation ),
            m_precision( precision )
        {            
            std::cout << std::setprecision( m_precision ) << std::scientific;
        }

        void WriteHeader()
        {
            std::cout << "Convergence Residuals:" << "\n";
        }

        void WriteResiduals( const FieldData<floatType> &residuals, 
                             const floatType massFluxResidual,
                             const intType nIterations )
        {
            using enum Axis::ENUMDATA;
            std::cout << "iteration: " << std::left << std::setw(5) << nIterations << ", "
                      << "U: " << residuals.U[ m_AT.CodeAxis(X) ] << ",   "
                      << "V: " << residuals.U[ m_AT.CodeAxis(Y) ] << ",   "
                      << "W: " << residuals.U[ m_AT.CodeAxis(Z) ] << ",   "
                      << "Local Coninuity: " << residuals.P             << ",   "
                      << "Global Mass Flux: " << massFluxResidual << "\n";
        }


    private:
        AxisTransformationMap m_AT;
        int m_precision;
};



// Postprocesses the fields and writes to file
class FieldWriter
{
    public:
        FieldWriter( const FieldData<Tensor3D> &fields, 
                     const Mesh &mesh,
                     const BoundaryConditionData &bcData,
                     const AxisTransformationMap &axisTransformation,
                     const std::string &baseFilename ) :
            m_fields( fields ),
            m_axisTransformation( axisTransformation ),
            m_transformedMesh( mesh ),
            m_transformedBcData( bcData ),
            m_baseFilename( IOTOOLS::RemoveFileExtension( baseFilename, ".vtk" ) )
            {
                // Only needs to be transformed once
                TransformMeshToUserCoordinates( m_transformedMesh, m_axisTransformation );
                TransformBCDataToUserCoordinates( m_transformedBcData, m_axisTransformation );
            };

            void WriteData( intType iterationNumber )
            { 
                TransformData();
                SetWriter();
                std::string message = "CFD solution at iteration " + std::to_string( iterationNumber );
                m_vtkWriter->WriteData( AppendFilename( iterationNumber ), message ); 
            }


    private:
        const FieldData<Tensor3D> &m_fields;
        const AxisTransformationMap &m_axisTransformation;
        FieldData<Tensor3D> m_transformedFields;
        FieldData<Tensor3D> m_transformedVertexFields;
        Mesh m_transformedMesh;
        BoundaryConditionData m_transformedBcData;
        std::unique_ptr< VTK::VTKWriter<floatType> > m_vtkWriter;
        const std::string m_baseFilename;

        std::string AppendFilename( intType iterationNumber )
        { return m_baseFilename  + "_iter" + std::to_string(iterationNumber) + ".vtk"; }

        void SetWriter()
        {
            using enum Axis::ENUMDATA;
            VTK::VTKWriterConfig config( m_transformedMesh.cellFaces[X].size(), 
                                         m_transformedMesh.cellFaces[Y].size(), 
                                         m_transformedMesh.cellFaces[Z].size() );
                config.SetWriteMode(VTK::WriteModes::BINARY);
                
            VTK::gridVectorType<CFD::floatType> gridVector = { m_transformedMesh.cellFaces[X].data(), 
                                                               m_transformedMesh.cellFaces[Y].data(), 
                                                               m_transformedMesh.cellFaces[Z].data() };

            VTK::scalarCollectionType<floatType> scalarMap = { {"Pressure", VTK::GridTypes::CELL_DATA, m_transformedFields.P.data()},
                                                               {"Pressure", VTK::GridTypes::POINT_DATA, m_transformedVertexFields.P.data()}};

            VTK::vectorCollectionType<floatType> vectorMap = { {"Velocity", VTK::GridTypes::CELL_DATA, { m_transformedFields.U[X].data(), 
                                                                                                         m_transformedFields.U[Y].data(), 
                                                                                                         m_transformedFields.U[Z].data()}},
                                                               {"Velocity", VTK::GridTypes::POINT_DATA, { m_transformedVertexFields.U[X].data(), 
                                                                                                          m_transformedVertexFields.U[Y].data(), 
                                                                                                          m_transformedVertexFields.U[Z].data()}} };
        
            m_vtkWriter = std::make_unique<VTK::VTKWriter<floatType>>(gridVector, scalarMap, vectorMap, config);
        }


        void TransformData()
        {
            ForAllFieldData([&](intType f) { 
                m_transformedFields[f] = FVT::RemoveGhostCells(m_fields[f], nGhost); 
            });
            TransformFieldToUserCoordinates( m_transformedFields, m_axisTransformation );
            m_transformedVertexFields = GetVertexFields(m_transformedFields, m_transformedMesh, m_transformedBcData);
        }

};


// Postprocesses the fields and writes to file
class ResidualFieldWriter
{
    public:
        ResidualFieldWriter( const FieldData<Tensor3D> &fields, 
                             const Mesh &mesh,
                             const AxisTransformationMap &axisTransformation,
                             const std::string &baseFilename ) :
            m_fields( fields ),
            m_axisTransformation( axisTransformation ),
            m_transformedMesh( mesh ),
            m_baseFilename(  IOTOOLS::RemoveFileExtension( baseFilename, ".vtk" ) )
            {
                // Only needs to be transformed once
                TransformMeshToUserCoordinates( m_transformedMesh, m_axisTransformation );
            };

            void WriteData( intType iterationNumber )
            { 
                TransformData();
                SetWriter();
                std::string message = "CFD residuals at iteration " + std::to_string( iterationNumber );
                m_vtkWriter->WriteData( AppendFilename( iterationNumber ), message ); 
            }


    private:
        const FieldData<Tensor3D> &m_fields;
        const AxisTransformationMap &m_axisTransformation;
        FieldData<Tensor3D> m_transformedFields;
        Mesh m_transformedMesh;
        std::unique_ptr< VTK::VTKWriter<floatType> > m_vtkWriter;
        const std::string m_baseFilename;

        std::string AppendFilename( intType iterationNumber )
        { return m_baseFilename  + "_iter" + std::to_string(iterationNumber) + ".vtk"; }

        void SetWriter()
        {
            using enum Axis::ENUMDATA;
            VTK::VTKWriterConfig config( m_transformedMesh.cellFaces[X].size(), 
                                         m_transformedMesh.cellFaces[Y].size(), 
                                         m_transformedMesh.cellFaces[Z].size() );
                config.SetWriteMode(VTK::WriteModes::BINARY);
                
            VTK::gridVectorType<CFD::floatType> gridVector = { m_transformedMesh.cellFaces[X].data(), 
                                                               m_transformedMesh.cellFaces[Y].data(), 
                                                               m_transformedMesh.cellFaces[Z].data() };

            VTK::scalarCollectionType<floatType> scalarMap = { {"Pressure", VTK::GridTypes::CELL_DATA, m_transformedFields.P.data()} };

            VTK::vectorCollectionType<floatType> vectorMap = { {"Velocity", VTK::GridTypes::CELL_DATA, { m_transformedFields.U[X].data(), 
                                                                                                         m_transformedFields.U[Y].data(), 
                                                                                                         m_transformedFields.U[Z].data()}} };
        
            m_vtkWriter = std::make_unique<VTK::VTKWriter<floatType>>(gridVector, scalarMap, vectorMap, config);
        }


        void TransformData()
        {
            m_transformedFields = m_fields;
            ForAllFieldData([&](intType f) { 
                m_transformedFields[f] = FVT::RemoveGhostCells(m_fields[f], nGhost); 
            });
            TransformFieldToUserCoordinates( m_transformedFields, m_axisTransformation );
        }

};


}   // end namespace CFD


#endif // CONVERGENCE_LOGGING
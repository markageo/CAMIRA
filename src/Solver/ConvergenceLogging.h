#ifndef SOLVER_LOGGING
#define SOLVER_LOGGING

#include "../Types.h"
#include "../Tools/SweepTransformations.h"
#include "../Tools/FVTools.h"

#include "Solver.h"
#include "FieldProbe.h"
#include "../IO/VTKWriter.h"

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

        void WriteResiduals( const FieldData<floatType> &residuals, 
                             const floatType massFluxResidual,
                             const intType nIterations )
        {
            using enum Axis::ENUMDATA;
            std::cout << "iteration: " << std::left << std::setw(5) << nIterations << ", "
                      << "U residual: " << residuals.U[ m_AT.CodeAxis(X) ] << ",   "
                      << "V residual: " << residuals.U[ m_AT.CodeAxis(Y) ] << ",   "
                      << "W residual: " << residuals.U[ m_AT.CodeAxis(Z) ] << ",   "
                      << "Coninuity residual: " << residuals.P             << ",   "
                      << "Global Mass residual: " << massFluxResidual << "\n";
        }


    private:
        AxisTransformationMap m_AT;
        int m_precision;
};



// Writes the raw field to file (including ghost nodes). This is intended for writing during the solution process, and is to
// be transformed and processed later.
class RawFieldWriter
{
    public:
        RawFieldWriter( const FieldData<array3D> &fields, 
                        const Mesh &mesh,
                        const std::string &baseFilename ) :
            m_baseFilename( baseFilename )
            {
                using enum Axis::ENUMDATA;

                // Add ghost cells (faces)
                Eigen::array<std::pair<int, int>, 1> paddings;
                paddings[0] = std::make_pair(nGhost, nGhost);
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    m_cellFaces[axis] = mesh.cellFaces[axis].pad( paddings );
                } );

                // Config for the writer
                VTK::VTKWriterConfig config(m_cellFaces[X].size(), m_cellFaces[Y].size(), m_cellFaces[Z].size());
                    config.SetWriteMode(VTK::WriteModes::BINARY);
                VTK::gridVectorType<CFD::floatType> gridVector = { m_cellFaces[X].data(), m_cellFaces[Y].data(), m_cellFaces[Z].data() };
                VTK::scalarMapType<CFD::floatType> scalarMap = { {"Pressure", VTK::GridTypes::CELL_DATA, fields.P.data()} };
                VTK::vectorMapType<CFD::floatType> vectorMap = { {"Velocity", VTK::GridTypes::CELL_DATA, {fields.U[X].data(), fields.U[Y].data(), fields.U[Z].data()}} };

                // Instantiate the writer
                m_vtkWriter = std::make_unique<VTK::VTKWriter<floatType>>(gridVector, scalarMap, vectorMap, config);

            };

            void WriteData( intType iterationNumber )
            { m_vtkWriter->WriteData( AppendFilename( iterationNumber ), "Raw solver field output with ghost cells." ); }


    private:
        EnumVector<Axis, array1D> m_cellFaces;  // The writer will have a reference to this
        std::unique_ptr< VTK::VTKWriter<floatType> > m_vtkWriter;
        const std::string m_baseFilename;

        std::string AppendFilename( intType iterationNumber )
        { return m_baseFilename + "_" + std::to_string(iterationNumber) + ".vtk"; }
};


}   // end namespace CFD


#endif // CONVERGENCE_LOGGING
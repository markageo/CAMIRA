#ifndef SOLVER_LOGGING
#define SOLVER_LOGGING

#include "Types.h"
#include "Solver.h"
#include "SweepTransformations.h"
#include "FieldProbe.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

namespace CFD
{

class ConvergenceFile
{
    using F = Fields::ENUMDATA;

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
    using F = Fields::ENUMDATA;

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
                                          "U residual", 
                                          "V residual",
                                          "W residual",
                                          "P residual", 
                                          "Global Mass Residual");
            
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
    using F = Fields::ENUMDATA;

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
            std::string fileDescription = "Probe: '" + fieldProbe.Name() + "', "
                                        + ", Location: (" 
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(X) ) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(Y) ) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate( m_AT.CodeAxis(Z) ) ) + ").";
            m_convergenceFile.WriteCommentLine( fileDescription );

            // Write header
            m_convergenceFile.WriteLine( "Iteration", 
                                         "X velocity", 
                                         "Y velocity",
                                         "Z velocity",
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
        EnumVector< Fields, F > m_userFields;
};





class ConsoleLog
{
    using F = Fields::ENUMDATA;

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
                      << "P residual: " << residuals.P                   << ",   "
                      << "Mass residual: " << massFluxResidual << "\n";
        }


    private:
        AxisTransformationMap m_AT;
        int m_precision;
};

}   // end namespace CFD


#endif // CONVERGENCE_LOGGING
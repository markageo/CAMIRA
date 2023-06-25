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

            } , ... )

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
            m_convergenceFile( filename, precision )
        {
            // Store the mapped fields
            EnumFor<Axis>( [&] ( Axis::ENUMDATA userAxis) {
                m_userFields[ AxisVelocity[ userAxis ] ] = AxisVelocity[ axisTransformation.CodeAxis( userAxis ) ];
            } );
            m_userFields[F::P] = F::P;
            
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

        void WriteData( const EnumVector<Fields, floatType> &residuals, 
                        const floatType massFluxResidual, 
                        const intType nIterations )
        {
            using enum Fields::ENUMDATA;

            m_convergenceFile.WriteLine( nIterations,
                                         residuals[ m_userFields[U] ], 
                                         residuals[ m_userFields[V] ],
                                         residuals[ m_userFields[W] ], 
                                         residuals[ m_userFields[P] ], 
                                         massFluxResidual );
        }


    private:
        ConvergenceFile m_convergenceFile;
        EnumVector< Fields, F > m_userFields;
};





class ProbeLogFile
{
    using F = Fields::ENUMDATA;

    public:
        ProbeLogFile( const std::string &filename, 
                      const AxisTransformationMap &axisTransformation,
                      const FieldProbe &fieldProbe,
                      const int precision = 6 ) :
            m_convergenceFile( filename, precision )
        {
            // Store the mapped fields
            EnumFor<Axis>( [&] ( Axis::ENUMDATA userAxis) {
                m_userFields[ AxisVelocity[ userAxis ] ] = AxisVelocity[ axisTransformation.CodeAxis( userAxis ) ];
            } );
            m_userFields[F::P] = F::P;
            
            // Comment description of file
            std::string fileDescription = "Probe: '" + fieldProbe.Name() + "', "
                                        + ", Location: (" 
                                        + std::to_string( fieldProbe.Coordinate(0) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate(1) ) + ", "
                                        + std::to_string( fieldProbe.Coordinate(2) ) + ").";
            m_convergenceFile.WriteCommentLine( fileDescription );

            // Write header
            m_convergenceFile.WriteLine( "Iteration", 
                                         "X velocity", 
                                         "Y velocity",
                                         "Z velocity",
                                         "Pressure");
            
        }

        void WriteData( const EnumVector<Fields, floatType> &probeValues,  
                                   const intType nIterations )
        {
            using enum Fields::ENUMDATA;

            m_convergenceFile.WriteLine( nIterations,
                                         probeValues[ m_userFields[U] ], 
                                         probeValues[ m_userFields[V] ],
                                         probeValues[ m_userFields[W] ], 
                                         probeValues[ m_userFields[P] ] );
        }


    private:
        ConvergenceFile m_convergenceFile;
        EnumVector< Fields, F > m_userFields;
};


}   // end namespace CFD


#endif // CONVERGENCE_LOGGING
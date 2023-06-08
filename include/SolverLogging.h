#ifndef SOVLER_LOGGING
#define SOLVER_LOGGING

#include "Types.h"
#include "Solver.h"
#include "SweepTransformations.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>

namespace CFD
{

class ConvergenceLogger
{
    using F = Fields::ENUMDATA;

    public:
        ConvergenceLogger( const std::string &filename, 
                           const AxisTransformationMap &axisTransformation,
                           const int precision = 6 ) :
            m_fileStream( filename ),
            m_precision( precision ),
            m_columnWidth( precision + 8 )
        {
            // Store the mapped fields
            EnumVector<Axis, F> axisField({ F::U, F::V, F::W });
            EnumFor<Axis>( [&] ( Axis::ENUMDATA userAxis) {
                m_userFields[ axisField[ userAxis ] ] = axisField[ axisTransformation.CodeAxis( userAxis ) ];
            } );
            m_userFields[F::P] = F::P;
            
            WriteFileHeader();

            m_fileStream << std::setprecision( m_precision ) << std::scientific;
            std::cout << std::setprecision( m_precision ) << std::scientific;
        }

        void WriteResidualsToFile( const EnumVector<Fields, floatType> &residuals, 
                                   const floatType massFluxResidual, 
                                   const intType nIterations )
        {
            m_fileStream << std::left << std::setw(m_columnWidth) << nIterations << ", "
                         << std::left << std::setw(m_columnWidth) << residuals[ m_userFields[F::U] ] << ", "
                         << std::left << std::setw(m_columnWidth) << residuals[ m_userFields[F::V] ] << ", "
                         << std::left << std::setw(m_columnWidth) << residuals[ m_userFields[F::W] ] << ", "
                         << std::left << std::setw(m_columnWidth) << residuals[ m_userFields[F::P] ] << ", "
                         << std::left << std::setw(m_columnWidth) << massFluxResidual                << "\n";
        }

        void WriteResidualsToScreen( const EnumVector<Fields, floatType> &residuals, 
                                     const floatType massFluxResidual,
                                     const intType nIterations )
        {
            std::cout << "iteration: " << std::left << std::setw(5) << nIterations << ", "
                      << "U residual: " << residuals[ m_userFields[F::U] ] << ",   "
                      << "V residual: " << residuals[ m_userFields[F::V] ] << ",   "
                      << "W residual: " << residuals[ m_userFields[F::W] ] << ",   "
                      << "P residual: " << residuals[ m_userFields[F::P] ] << ",   "
                      << "Mass residual: " << massFluxResidual << "\n";
        }


    private:
        std::ofstream m_fileStream;
        EnumVector< Fields, F > m_userFields;
        int m_precision;
        int m_columnWidth;

        void WriteFileHeader()
        {
            m_fileStream << std::left << std::setw(m_columnWidth) << "Iteration" << ", "
                         << std::left << std::setw(m_columnWidth) << "U_residual" << ", "
                         << std::left << std::setw(m_columnWidth) << "V_residual" << ", "
                         << std::left << std::setw(m_columnWidth) << "W_residual" << ", "
                         << std::left << std::setw(m_columnWidth) << "P_residual" << ", "
                         << std::left << std::setw(m_columnWidth) << "Mass_residual" << "\n";
        }

};


}   // end namespace CFD


#endif // SOLVER_LOGGING
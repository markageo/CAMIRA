// #ifndef LOGGING
// #define LOGGING

// #include "Types.h"
// #include "Solver.h"

// #include <fstream>
// #include <iomanip>
// #include <iostream>
// #include <sstream>
// #include <string>


// namespace CFD
// {

// // Write convergence data to file
// void WriteConvergenceData( const std::string &filename, 
//                            const ConvergenceData &convergenceData)
// {
//     using F = Fields::ENUMDATA;

//     std::ofstream fileStream(filename);

//     // Header
//     fileStream << "U_res" << ","
//                << "V_res" << ","
//                << "W_res" << ","
//                << "P_res" << "\n";            

//     // Data
//     for ( const auto &residuals : convergenceData.residuals ) {

//         fileStream << residuals[F::U] << ", "
//                    << residuals[F::V] << ", "
//                    << residuals[F::W] << ", "
//                    << residuals[F::P] << "\n";
//     }

// }

// }   // end namespace CFD


// #endif // LOGGING
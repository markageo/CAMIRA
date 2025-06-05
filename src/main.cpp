/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Core/Types.h"
#include "Core/Macros.h"
#include "Core/FVTools.h"
#include "Geometry/Geometry.h"
#include "ImmersedBoundary/ImmersedBoundary.h"
#include "IO/ArrayIO.h"
#include "IO/InputProcessing.h"
#include "IO/VTKWriter.h"
#include "CoordinateTransformations/AxisTransformationFunctions.h"
#include "FiniteVolume/FiniteVolume.h"
#include "Parallel/Parallel.h"
#include "Solver/Solver.h"

#include <fstream>
#include <iostream>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    CFD::AxisTransformationMap axisTransformation = CreateAxisTransformation( inputData.smootherSettings.planeSweepDirection,
                                                                              inputData.smootherSettings.lineSweepDirection );

    CFD::TransformUserInputData(inputData, axisTransformation);

    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    // std::cout << "Press enter to begin solve.";
    // std::cin.ignore();
    // std::cout << std::endl;

    omp_set_num_threads( inputData.parallelSettings.numberOfThreads );

    switch ( inputData.schemes.momentumInterpolation ) {

        using MI = CFD::MomentumInterpolation;

        case ( MI::Implicit ):

            if ( inputData.transient ) {
                CFD::SolveTransient< MI::Implicit >(inputData, axisTransformation);
            } else {
                CFD::SolveSteady< MI::Implicit >(inputData, axisTransformation);
            }
            break;

        case ( MI::SemiExplicit ):

            if ( inputData.transient ) {
                CFD::SolveTransient< MI::SemiExplicit >(inputData, axisTransformation);
            } else {
                CFD::SolveSteady< MI::SemiExplicit >(inputData, axisTransformation);
            }

            break;
    }


    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/


    // Display profiling info to terminal and write to file
    #ifdef CFD_PROFILING
        std::cout << PROF::prof;

        std::ofstream profilingFilestream( inputData.profilingFilename );
        profilingFilestream << PROF::prof;
    #endif
    
    return 0;
}
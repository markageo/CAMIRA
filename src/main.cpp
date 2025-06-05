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

    CAMIRA::InputData inputData = CAMIRA::InputDataFromCommandLine(argc, argv);

    CAMIRA::AxisTransformationMap axisTransformation = CreateAxisTransformation( inputData.smootherSettings.planeSweepDirection,
                                                                              inputData.smootherSettings.lineSweepDirection );

    CAMIRA::TransformUserInputData(inputData, axisTransformation);

    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    // std::cout << "Press enter to begin solve.";
    // std::cin.ignore();
    // std::cout << std::endl;

    omp_set_num_threads( inputData.parallelSettings.numberOfThreads );

    switch ( inputData.schemes.momentumInterpolation ) {

        using MI = CAMIRA::MomentumInterpolation;

        case ( MI::Implicit ):

            if ( inputData.transient ) {
                CAMIRA::SolveTransient< MI::Implicit >(inputData, axisTransformation);
            } else {
                CAMIRA::SolveSteady< MI::Implicit >(inputData, axisTransformation);
            }
            break;

        case ( MI::SemiExplicit ):

            if ( inputData.transient ) {
                CAMIRA::SolveTransient< MI::SemiExplicit >(inputData, axisTransformation);
            } else {
                CAMIRA::SolveSteady< MI::SemiExplicit >(inputData, axisTransformation);
            }

            break;
    }


    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/


    // Display profiling info to terminal and write to file
    #ifdef CAMIRA_PROFILING
        std::cout << PROF::prof;

        std::ofstream profilingFilestream( inputData.profilingFilename );
        profilingFilestream << PROF::prof;
    #endif
    
    return 0;
}
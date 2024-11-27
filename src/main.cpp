/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "Macros.h"

#include "Geometry/Geometry.h"
#include "ImmersedBoundary/ImmersedBoundary.h"
#include "IO/ArrayIO.h"
#include "IO/InputProcessing.h"
#include "IO/VTKWriter.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"
#include "FiniteVolume/FiniteVolume.h"
#include "Solver/Solver.h"

#include <fstream>
#include <iostream>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    TIC("Input Processing")
    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    CFD::AxisTransformationMap axisTransformation = CreateAxisTransformation( inputData.linearSolverSettings.planeSweepDirection,
                                                                              inputData.linearSolverSettings.lineSweepDirection );

    CFD::TransformUserInputData(inputData, axisTransformation);
    TOC()

    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    std::cout << "Press enter to begin solve.";
    std::cin.ignore();
    std::cout << std::endl;


    TIC("Solver");
    switch ( inputData.schemes.momentumInterpolation ) {

        using MI = CFD::MomentumInterpolation;
        using LI = CFD::Linearisation;

        case ( MI::Implicit ):
           switch ( inputData.schemes.linearisation ) {
                case ( LI::Picard ):
                    if ( inputData.transient ) {
                        CFD::SolveTransient< MI::Implicit, LI::Picard >(inputData, axisTransformation);
                    } else {
                        CFD::SolveSteady< MI::Implicit, LI::Picard >(inputData, axisTransformation);
                    }
                    break;

                case ( LI::Newton ):
                    if ( inputData.transient ) {
                        CFD::SolveTransient< MI::Implicit, LI::Newton >(inputData, axisTransformation);
                    } else {
                        CFD::SolveSteady< MI::Implicit, LI::Newton >(inputData, axisTransformation);
                    }
                    
                    break;
            }
            break;

        case ( MI::SemiExplicit ):
            switch ( inputData.schemes.linearisation ) {
                case ( LI::Picard ):
                    if ( inputData.transient ) {
                        CFD::SolveTransient< MI::SemiExplicit, LI::Picard >(inputData, axisTransformation);
                    } else {
                        CFD::SolveSteady< MI::SemiExplicit, LI::Picard >(inputData, axisTransformation);
                    }
                    break;

                case ( LI::Newton ):
                    if ( inputData.transient ) {
                        CFD::SolveTransient< MI::SemiExplicit, LI::Newton >(inputData, axisTransformation);
                    } else {
                        CFD::SolveSteady< MI::SemiExplicit, LI::Newton >(inputData, axisTransformation);
                    }
                    
                    break;
            }
            break;
    }
    TOC();

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
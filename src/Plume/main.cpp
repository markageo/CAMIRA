/*---------------------------------------------------------------------------*\
   3D Lagrangian Particle Dispersion Solver

   Mark George
\*---------------------------------------------------------------------------*/

#include "Core/Types.h"
#include "Core/Macros.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/IO/ArrayIO.h"
#include "Core/IO/VTKWriter.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/Solver/Solver.h"

#include <fstream>
#include <iostream>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    CAMIRA::PLUME::InputData inputData = CAMIRA::PLUME::InputDataFromCommandLine(argc, argv);



    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    omp_set_num_threads( inputData.parallelSettings.numberOfThreads );

    CAMIRA::PLUME::SolvePlume( inputData );



    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/

    


    return 0;
}
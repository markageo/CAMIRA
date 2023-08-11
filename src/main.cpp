/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "Macros.h"

#include "IO/InputProcessing.h"
#include "IO/VTKWriter.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"
#include "FiniteVolume/FiniteVolume.h"
#include "Solver/Solver.h"

#include <iostream>
#include <fmt/core.h>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    TIC("Input Processing")
    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    CFD::AxisTransformationMap axisTransformation = CFD::TransformUserInputData(inputData);
    TOC()

    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    TIC("Meshing");
    CFD::Mesh mesh(inputData);
    TOC();

    TIC("Field Allocation");
    CFD::FieldData<CFD::array3D> fields = CFD::InitialiseFields(mesh, inputData);
    TOC();

    TIC("Boundary Condition Processing")
    CFD::FieldData<CFD::BoundaryConditionData> bcData = SetBoundaryConditionData(inputData, mesh);
    TOC()

    TIC("Solver");
    switch (inputData.schemes.momentumInterpolation)
    {
        using MI = CFD::MomentumInterpolation;
    case (MI::Implicit):
        CFD::SweepSolve<MI::Implicit>(fields, mesh, bcData, inputData, axisTransformation);
        break;

    case (MI::SemiExplicit):
        CFD::SweepSolve<MI::SemiExplicit>(fields, mesh, bcData, inputData, axisTransformation);
        break;
    }
    TOC();


    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/


    // Display profiling information
    #ifdef CFD_PROFILING
        std::cout << PROF::prof;
    #endif
    
    return 0;
}
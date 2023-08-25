/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "Macros.h"

#include "Geometry/Geometry.h"
#include "ImmersedBoundary/ImmersedBoundary.h"
#include "IO/ArrayIO.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/draw_polyhedron.h>
#include <fstream>

#include "IO/InputProcessing.h"
#include "IO/VTKWriter.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"
#include "FiniteVolume/FiniteVolume.h"
#include "Solver/Solver.h"

#include <iostream>

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
    CFD::FieldData<CFD::Tensor3D> fields = CFD::InitialiseFields(mesh, inputData);
    TOC();

    TIC("Boundary Condition Processing")
    CFD::FieldData<CFD::BoundaryConditionData> bcData = SetBoundaryConditionData(inputData, mesh);
    TOC()

    TIC("Geometry Creation")

    typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
    typedef CGAL::Polyhedron_3<Kernel>                       Polyhedron;
    
    CFD::Polyhedron P = CFD::MakeGeometry( inputData );

    // CGAL::draw( P );

    CFD::CellIDTensor3D cellID =  TagCells( mesh, P);

    CFD::Tensor3D cellIDCasted = cellID.cast<CFD::floatType>();
    CFD::WriteArray( "cellID.dbg", cellIDCasted );

    CFD::IBData ibData = CFD::CreateImmersedBoundaryData( P, cellID, mesh );
    
    TOC()

    // TIC("Solver");
    // switch ( inputData.schemes.momentumInterpolation ) {

    //     using MI = CFD::MomentumInterpolation;
    //     using LI = CFD::Linearisation;

    //     case ( MI::Implicit ):
    //        switch ( inputData.schemes.linearisation ) {
    //             case ( LI::Picard ):
    //                 CFD::SweepSolve< MI::Implicit, LI::Picard >(fields, mesh, bcData, inputData, axisTransformation);
    //                 break;

    //             case ( LI::Newton ):
    //                 CFD::SweepSolve< MI::Implicit, LI::Newton >(fields, mesh, bcData, inputData, axisTransformation);
    //                 break;
    //         }
    //         break;

    //     case ( MI::SemiExplicit ):
    //         switch ( inputData.schemes.linearisation ) {
    //             case ( LI::Picard ):
    //                 CFD::SweepSolve< MI::SemiExplicit, LI::Picard >(fields, mesh, bcData, inputData, axisTransformation);
    //                 break;

    //             case ( LI::Newton ):
    //                 CFD::SweepSolve< MI::SemiExplicit, LI::Newton >(fields, mesh, bcData, inputData, axisTransformation);
    //                 break;
    //         }
    //         break;
    // }
    // TOC();

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
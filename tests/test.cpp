/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids
    TESTS

    Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "utils.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "FiniteVolumeFunctions.h"
#include "VTKWriter.h"
#include "Solver.h"

#include <iostream>
#include <fmt/core.h>


int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    /*-------------------------------------------------------------------------------------*\
                                          Test Functions
    \*-------------------------------------------------------------------------------------*/

    using AX = CFD::Axis::ENUMDATA;
    using F = CFD::Fields::ENUMDATA;

    TIC("Meshing");
    const CFD::Mesh mesh(inputData);
    TOC();

    TIC("Field Allocation");
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells);
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields2({F::U, F::V, F::W, F::P}, mesh.nCells);

    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> faceVelocities( {{F::U, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                                    {F::V, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                                                    {F::W, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
    TOC();

    TIC("Set Values");
    fields[F::U].setRandom();
    fields[F::V].setRandom();
    fields[F::W].setRandom();
    TOC();

    TIC("Face Velocity Update")
    CFD::UpdateFaceVelocities(faceVelocities, mesh, fields, inputData.boundaryConditions);
    TOC();

    TIC("Allocate and Initialise FV Coefficients")
    CFD::FVCoefficients fvCoeffs = CFD::InitialiseFVCoefficients(mesh, faceVelocities, inputData);
    TOC()

    // // Cell centers to file
    // UTIL::writeArray("debug/cell_centers_x.dat", mesh.cellCenters[AX::X]);
    // UTIL::writeArray("debug/cell_centers_y.dat", mesh.cellCenters[AX::Y]);
    // UTIL::writeArray("debug/cell_centers_z.dat", mesh.cellCenters[AX::Z]);

    // // Cell centers to file
    // UTIL::writeArray("debug/cell_faces_x.dat", mesh.cellFaces[AX::X]);
    // UTIL::writeArray("debug/cell_faces_y.dat", mesh.cellFaces[AX::Y]);
    // UTIL::writeArray("debug/cell_faces_z.dat", mesh.cellFaces[AX::Z]);

    // // Write extrapolation factors to a file
    // UTIL::writeArray("debug/interp_factors_x.dat", mesh.interpFactors[AX::X]);
    // UTIL::writeArray("debug/interp_factors_y.dat", mesh.interpFactors[AX::Y]);
    // UTIL::writeArray("debug/interp_factors_z.dat", mesh.interpFactors[AX::Z]);

    // // Write fields to a file
    // UTIL::writeArray("debug/U_cell_centers.dat", fields[F::U]);
    // UTIL::writeArray("debug/U_cell_faces.dat", faceVelocities[F::U]);

    // UTIL::writeArray("debug/V_cell_centers.dat", fields[F::V]);
    // UTIL::writeArray("debug/V_cell_faces.dat", faceVelocities[F::V]);

    // UTIL::writeArray("debug/W_cell_centers.dat", fields[F::W]);
    // UTIL::writeArray("debug/W_cell_faces.dat", faceVelocities[F::W]);

    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    // Data to pass to writer
    VTK::dataType VTKDataType = VTK::DOUBLE;
    if (std::is_same<CFD::floatType, float>::value)
    {
        VTKDataType = VTK::FLOAT;
    }
    VTK::VTKWriterConfig config(mesh.nCells[AX::X], mesh.nCells[AX::Y], mesh.nCells[AX::Z], VTKDataType);
        config.SetWriteMode("binary");
    VTK::gridVectorType gridVector = {mesh.cellCenters[AX::X].data(), mesh.cellCenters[AX::Y].data(), mesh.cellCenters[AX::Z].data()};
    VTK::scalarMapType scalarMap = {{"U", fields[F::U].data()}};
    VTK::vectorMapType vectorMap = {};

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    TIC("Writer");
    writer.WriteData("mesh.vtk", "3D Rectilinear Grid");
    TOC();

// Display profiling information
#ifdef PROFILING
    std::cout << PROF::prof;
#endif

    return 0;
}

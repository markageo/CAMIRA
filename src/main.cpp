/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "VTKWriter.h"
#include "Solver.h"
#include "utils.h"

#include "FiniteVolumeFunctions.h"


#include <iostream>
#include <type_traits>


int main(int argc, char const *argv[]) 
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    // Input file from first command line argument
    std::string inputFilename;
    if (argc == 1) {
        std::cout << "Please enter input file name: " << "\n";
        std::cin >> inputFilename;
        std::cout << "\n";
        std::cin.ignore();
    } else if (argc == 2) {
        inputFilename = argv[1];
    } else {
        throw std::invalid_argument( "Invalid command line options." );
    }

    // User input data
    std::string inputFileRetryChoice;
    CFD::InputData inputData;
    while (true){
        try {

            std::cout << "Reading input file: '" + inputFilename + "' ..." << "\n\n";
            inputData = CFD::ReadInputData(inputFilename);
            std::cout << "Success!" << "\n\n";
            break;

        // } catch (std::runtime_error &e) {
        } catch (int &e) {

            // std::cout << "Failure reading input file! \n" 
            //           << e.what() 
            //           << "\n\n";

            std::cout << "Would you like to try again? (y/n)" 
                      << "\n";
            std::cin  >> inputFileRetryChoice;
            if (inputFileRetryChoice != "y") 
                exit(-1);
            std::cout << "\n";

            std::cout << "Please enter input file name: " << "\n";
            std::cin  >> inputFilename;
            std::cout << "\n";
            std::cin.ignore();

        }
    }
    std::cout << "Press enter to begin.";
    std::cin.ignore();
    std::cout << std::endl;


    /*-------------------------------------------------------------------------------------*\
                                           Solve
    \*-------------------------------------------------------------------------------------*/

    using F = CFD::Fields::ENUMDATA;
    using BC = CFD::BoundaryConditions::ENUMDATA;
    using BP = CFD::BoundaryPatches::ENUMDATA;

    const CFD::Mesh mesh(inputData);
    CFD::ArrayAllocator<F> fields({F::U, F::V, F::W, F::P}, mesh.nCells);

    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    CFD::ArrayAllocator<F> faceVelocities( { {F::U, {mesh.nCells(0)+1, mesh.nCells(1)  , mesh.nCells(2)  }}, 
                                             {F::V, {mesh.nCells(0)  , mesh.nCells(1)+1, mesh.nCells(2)  }}, 
                                             {F::W, {mesh.nCells(0)  , mesh.nCells(1)  , mesh.nCells(2)+1}} } );



    fields[F::U].setRandom();
    fields[F::V].setRandom();
    fields[F::W].setRandom();
    CFD::UpdateFaceVelocities(faceVelocities, mesh, fields, inputData.boundaryConditions);
    
    UTIL::writeArray("U_cell_centers.txt", fields[F::U]);
    UTIL::writeArray("U_cell_faces.txt", faceVelocities[F::U]);

    UTIL::writeArray("V_cell_centers.txt", fields[F::V]);
    UTIL::writeArray("V_cell_faces.txt", faceVelocities[F::V]);

    UTIL::writeArray("W_cell_centers.txt", fields[F::W]);
    UTIL::writeArray("W_cell_faces.txt", faceVelocities[F::W]);

    

    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    // Data to pass to writer
    VTK::dataType VTKDataType = VTK::DOUBLE;
    if (std::is_same<CFD::floatType, float>::value) {
        VTKDataType = VTK::FLOAT;
    } 
    VTK::VTKWriterConfig config( mesh.cellCenters_x.size(), mesh.cellCenters_y.size(), mesh.cellCenters_z.size(), VTKDataType);
        config.SetWriteMode("ascii");
        config.SetASCIIPrecision(8);
    VTK::gridVectorType gridVector = {mesh.cellCenters_x.data(), mesh.cellCenters_y.data(), mesh.cellCenters_z.data()};
    VTK::scalarMapType scalarMap = { };
    VTK::vectorMapType vectorMap = { };

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    writer.WriteData("mesh.vtk", "3D Rectilinear Grid");

    /*-------------------------------------------------------------------------------------*\
                                           Testing
    \*-------------------------------------------------------------------------------------*/


    return 0;
}

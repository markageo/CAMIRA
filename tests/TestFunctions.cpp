#include "TestFunctions.h"

#include "Types.h"
#include "utils.h"
#include <filesystem>


// Write mesh data to files
void TEST::WriteMesh(const CFD::Mesh &mesh, const std::string &filedir) {

    using AX = CFD::Axis::ENUMDATA;

    AX axis;
    std::string ext = ".dat";
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"_x", "_y", "_z"} };
    std::string fname;
    auto Fname = [&] (const std::string &s) -> std::string { return filedir + s + axisSuffix[axis] + ext; };
    std::filesystem::create_directory(filedir);

    for (int a = 0; a != 3; a++) {
        axis = static_cast<AX>(a);

        // Cell Lengths
        UTIL::writeArray(Fname("cell_lengths"), mesh.cellLengths[axis]);   

        // Inverse cell lengths
        UTIL::writeArray(Fname("inv_cell_lengths"), mesh.cellLengthsInv[axis]);   

        // Cell centers
        UTIL::writeArray(Fname("cell_centers"), mesh.cellCenters[AX::X]);  

        // Inverse of cell center difference
        UTIL::writeArray(Fname("inv_cell_center_diff"), mesh.cellCenterDiffInv[AX::X]);

        // Cell faces
        UTIL::writeArray(Fname("cell_faces"), mesh.cellFaces[AX::X]);

        // Write interpolation factors to a file
        UTIL::writeArray(Fname("interp_factors"), mesh.interpFactors[AX::X]);

    }

}
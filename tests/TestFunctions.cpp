#include "TestFunctions.h"

#include "Types.h"
#include "utils.h"
#include <filesystem>
#include <iostream>


// ReadConfig helper functions
namespace 
{

    // Remove whitespace from string
    void RemoveWhiteSpace(std::string &str) {
        str.erase(std::remove_if(str.begin(), str.end(), isspace), str.end());
    }

    // Wrapper for std::getline that removes whitespace
    std::ifstream &ReadInputLine(std::ifstream &ifs, std::string &str) 
    {
        std::getline(ifs, str);
        RemoveWhiteSpace(str);
        return ifs;
    }

    // Returns substring before equals sign. Reads to end of string if no equals is present.
    std::string ReadBeforeEquals(const std::string &str) 
    {
        std::string ostr; // part of string before equals sign
        for (auto istr = str.begin(); istr != str.end() && *istr != '='; ++istr) { ostr += *istr; }
        return ostr;
    }


    // Returns substring after equals sign. Returns empty string if equals not found.
    std::string ReadAfterEquals(const std::string &str) 
    {
        std::string ostr; // part of string after equals sign
        auto istr = str.begin();
        while (*istr != '=' && istr != str.end()) { ++istr; } // Find start of equals
        istr++; // The character after the equals
        for (/* NULL */; istr != str.end() && *istr != '#'; ++istr) { ostr += *istr; }  // Get values after equals
        return ostr;
    }

    // Convert string to given numeric type T.
    template <typename T> 
    T String2Type(const std::string &str)
    {
        // NOTE: This does not work for ints in scientific notation.
        std::istringstream strstream(str);
        T num;
        strstream >> num;
        return num;
    }

    // Convert a string to a test state
    TEST::testState String2TestState(const std::string &str)
    {
        if        (str == "test") {
            return TEST::test;
        } else if (str == "write") {
            return TEST::write;
        }
        return TEST::none;
    }

}


TEST::TestConfig TEST::ReadConfig(const std::string &filename)
{

    TEST::TestConfig testConfig;
    std::ifstream fileStream(filename);

    std::string line;
    std::string key;
    std::string value;
    while (ReadInputLine(fileStream, line)) {

        if (line[0] == '#' || line.size() == 0) {
            continue;
        }

        key = ReadBeforeEquals(line);
        value = ReadAfterEquals(line);
        if         (key == "testInputDirectory") {
            testConfig.testInputDirectory = String2Type<std::string>(value);

        } else if (key == "meshTestInputFilename") {
            testConfig.meshTestInputFilename = String2Type<std::string>(value);

        } else if (key == "meshTest") {
            testConfig.meshTest = String2TestState(value);

        } else if (key == "meshTestOutputDirectory") {
            testConfig.meshTestOutputDirectory = String2Type<std::string>(value);

        } else if (key == "meshTestReferenceDirectory") {
            testConfig.meshTestReferenceDirectory = String2Type<std::string>(value);

        }

    }

    return testConfig;
}


void TEST::WriteMesh(const CFD::Mesh &mesh, const std::string &filedir) 
{

    using AX = CFD::Axis::ENUMDATA;

    AX axis;
    std::string ext = ".dat";
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"_x", "_y", "_z"} };
    std::string fname;
    auto Fname = [&] (const std::string &s) -> std::string { return filedir + s + axisSuffix[axis] + ext; };
    std::filesystem::create_directory(filedir);

    // Temporary array for extrapolation factors
    // rows  : extrapFactors.p, extrapFactors.a
    // cols  : positivePatch, negativepatch
    CFD::array2D extrapFactorsArray(2, 2); 

    for (int a = 0; a != 3; a++) {
        axis = static_cast<AX>(a);

        // Cell Lengths
        UTIL::WriteArray(Fname("cell_lengths"), mesh.cellLengths[axis]);   

        // Inverse cell lengths
        UTIL::WriteArray(Fname("inv_cell_lengths"), mesh.cellLengthsInv[axis]);   

        // Cell centers
        UTIL::WriteArray(Fname("cell_centers"), mesh.cellCenters[axis]);  

        // Inverse of cell center difference
        UTIL::WriteArray(Fname("inv_cell_center_diff"), mesh.cellCenterDiffInv[axis]);

        // Cell faces
        UTIL::WriteArray(Fname("cell_faces"), mesh.cellFaces[axis]);

        // Write interpolation factors to a file
        UTIL::WriteArray(Fname("interp_factors"), mesh.interpFactors[axis]);

        // Put the extrapolation factors in an array and write it out
        extrapFactorsArray(0, 0) = mesh.extrapFactors[CFD::positivePatches[axis]].p;
        extrapFactorsArray(1, 0) = mesh.extrapFactors[CFD::positivePatches[axis]].a;
        extrapFactorsArray(0, 1) = mesh.extrapFactors[CFD::negativePatches[axis]].p;
        extrapFactorsArray(1, 1) = mesh.extrapFactors[CFD::negativePatches[axis]].a;
        UTIL::WriteArray(Fname("extrap_factors"), extrapFactorsArray);

    }

}


bool TEST::CompareMesh(const CFD::Mesh &, const std::string &)
{


    return false;

}
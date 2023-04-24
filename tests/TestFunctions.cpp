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

}   // end anonymous namespace


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
            testConfig.testInputDirectory = value;

        } else if (key == "meshTestInputFilename") {
            testConfig.meshTestInputFilename = value;

        } else if (key == "meshTest") {
            testConfig.meshTest = String2TestState(value);

        } else if (key == "meshTestOutputDirectory") {
            testConfig.meshTestOutputDirectory = value;

        } else if (key == "meshTestReferenceDirectory") {
            testConfig.meshTestReferenceDirectory = value;

        } else if (key == "faceVelUniformTestInputFilename") {
            testConfig.faceVelUniformTestInputFilename = value;

        } else if (key == "faceVelZeroGradientTestInputFilename") {
            testConfig.faceVelZeroGradientTestInputFilename = value;

        } else if (key == "faceVelExtrapolatedTestInputFilename") {
            testConfig.faceVelExtrapolatedTestInputFilename = value;

        } else if (key == "faceVelTest") {
            testConfig.faceVelTest = String2TestState(value);

        } else if (key ==  "faceVelTestOutputDirectory") {
            testConfig.faceVelTestOutputDirectory = value;

        } else if (key == "faceVelTestReferenceDirectory") {
            testConfig.faceVelTestReferenceDirectory = value;
        }

    }

    return testConfig;
}


void TEST::WriteMesh(const CFD::Mesh &mesh, const std::string &filedir, const int precision) 
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
        UTIL::WriteArray(Fname("cell_lengths"), mesh.cellLengths[axis], precision);   

        // Inverse cell lengths
        UTIL::WriteArray(Fname("inv_cell_lengths"), mesh.cellLengthsInv[axis], precision);   

        // Cell centers
        UTIL::WriteArray(Fname("cell_centers"), mesh.cellCenters[axis], precision);  

        // Inverse of cell center difference
        UTIL::WriteArray(Fname("inv_cell_center_diff"), mesh.cellCenterDiffInv[axis], precision);

        // Cell faces
        UTIL::WriteArray(Fname("cell_faces"), mesh.cellFaces[axis], precision);

        // Write interpolation factors to a file
        UTIL::WriteArray(Fname("interp_factors"), mesh.interpFactors[axis], precision);

        // Put the extrapolation factors in an array and write it out
        extrapFactorsArray(0, 0) = mesh.extrapFactors[CFD::positivePatches[axis]].p;
        extrapFactorsArray(1, 0) = mesh.extrapFactors[CFD::positivePatches[axis]].a;
        extrapFactorsArray(0, 1) = mesh.extrapFactors[CFD::negativePatches[axis]].p;
        extrapFactorsArray(1, 1) = mesh.extrapFactors[CFD::negativePatches[axis]].a;
        UTIL::WriteArray(Fname("extrap_factors"), extrapFactorsArray, precision);

    }

}


// Helper functions Comparision tests
namespace 
{

    // Returns true if all elements of a tensor are exactly the same
    template<typename T>
    bool ArraysSame(const T &array1, const T &array2) 
    {
        using dimType = long int;

        // Check their dimensions are the same
        if (array1.size() != array2.size())
            return false;

        // Loop through and compare
        dimType n = array1.size();
        auto *array1Data = array1.data();
        auto *array2Data = array2.data();
        for (dimType i = 0; i != n; i++) {
            if ( array1Data[i] != array2Data[i] )
                return false;
        }
        return true;
    }


    // Read and compare arrays. Returns true if they are the same
    template<typename T>
    bool ReadAndCompare(const std::string &filename1, const std::string &filename2) {

        // Temporaries for storing read in data
        T array1, array2;
        array1 = UTIL::ReadArray<CFD::array1D>(filename1);
        array2 = UTIL::ReadArray<CFD::array1D>(filename2);
        if (!ArraysSame(array1, array2))
            return false;

        return true;

    }


}   // end anonymous namespace


bool TEST::CompareMesh(const std::string &outputDir, const std::string &referenceDir)
{

    using AX = CFD::Axis::ENUMDATA;

    AX axis;
    std::string ext = ".dat";
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"_x", "_y", "_z"} };
    std::string fname;
    auto Fname = [&] (const std::string &name, const std::string &dir) -> std::string { return dir + name + axisSuffix[axis] + ext; };

    // Temporary array for extrapolation factors
    // rows  : extrapFactors.p, extrapFactors.a
    // cols  : positivePatch, negativepatch
    CFD::array2D extrapFactorsArray(2, 2); 

    // Temporary variable for output and reference directories
    CFD::array1D outputArray1D, referenceArray1D;
    CFD::array2D outputArray2D, referenceArray2D;

    for (int a = 0; a != 3; a++) {
        axis = static_cast<AX>(a);

        // Cell Lengths
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_lengths", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_lengths", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;

        // Inverse cell lengths
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("inv_cell_lengths", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("inv_cell_lengths", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;  

        // Cell centers
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_centers", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_centers", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;

        // Inverse of cell center difference
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("inv_cell_center_diff", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("inv_cell_center_diff", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;

        // Cell faces
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_faces", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("cell_faces", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;

        // Write interpolation factors to a file
        outputArray1D = UTIL::ReadArray<CFD::array1D>(Fname("interp_factors", outputDir));
        referenceArray1D = UTIL::ReadArray<CFD::array1D>(Fname("interp_factors", referenceDir));
        if (!ArraysSame(outputArray1D, referenceArray1D))
            return false;

        // Put the extrapolation factors in an array and write it out
        outputArray2D = UTIL::ReadArray<CFD::array2D>(Fname("extrap_factors", outputDir));
        referenceArray2D = UTIL::ReadArray<CFD::array2D>(Fname("extrap_factors", referenceDir));
        if (!ArraysSame(outputArray2D, referenceArray2D))
            return false;


    }

    return true;

}
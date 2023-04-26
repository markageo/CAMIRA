#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "FiniteVolumeStructures.h"
#include <string>

namespace TEST {

// File extension for output files
constexpr char testFileExtension[] = ".dat";

enum testState {
    none,
    write,
    test
};

struct TestConfig
{
    std::string testInputDirectory;

    // Mesh testing
    std::string meshTestInputFilename;
    testState meshTest;
    std::string meshTestOutputDirectory;
    std::string meshTestReferenceDirectory;

    // General testing that involves face velocities
    std::string velUniformTestInputFilename;
    std::string velZeroGradientTestInputFilename;
    std::string velExtrapolatedTestInputFilename;

    // Face velocity testing
    testState faceVelTest;
    std::string faceVelTestOutputDirectory;
    std::string faceVelTestReferenceDirectory;

    // Finite volume coefficient testing
    testState fvCoeffTest;
    std::string fvCoeffTestOutputDirectory;
    std::string fvCoeffTestReferenceDirectory;

};

// Read test config files
TestConfig ReadConfig(const std::string &);


// Write sets of data to file
void WriteMesh(const CFD::Mesh &, const std::string &, const int = 6);
void WriteFields(const CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &, const std::string &, const int = 6);
void WriteFaceVels(const CFD::EnumVector<CFD::BoundaryConditions, CFD::ArrayAllocator<CFD::Fields, CFD::array3D> > &, const std::string &, const int = 6);
void WriteFVCoeffs(const CFD::EnumVector<CFD::BoundaryConditions, CFD::FVCoefficients > &, const std::string &, const int = 6);

// Compare data in two different directories
bool CompareMesh(const std::string &, const std::string &);
bool CompareFaceVels(const std::string &, const std::string &);
bool CompareFVCoeffs(const std::string &, const std::string &);


}   // end namespace TEST


#endif // TEST_FUNCTIONS    
#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "FiniteVolumeStructures.h"

namespace TEST {

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

    // Face velocity testing
    std::string faceVelUniformTestInputFilename;
    std::string faceVelZeroGradientTestInputFilename;
    std::string faceVelExtrapolatedTestInputFilename;
    testState faceVelTest;
    std::string faceVelTestOutputDirectory;
    std::string faceVelTestReferenceDirectory;
};

// Read test config files
TestConfig ReadConfig(const std::string &);


// Write sets of data to file
void WriteMesh(const CFD::Mesh &, const std::string &, const int = 6);
void WriteFields(const CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &, const std::string &, const int = 6);
void WriteFaceVels(const CFD::EnumVector<CFD::BoundaryConditions, CFD::ArrayAllocator<CFD::Fields, CFD::array3D> > &, const std::string &, const int = 6);

// Compare data in two different directories
bool CompareMesh(const std::string &, const std::string &);
bool CompareFaceVels(const std::string &, const std::string &);


}   // end namespace TEST


#endif // TEST_FUNCTIONS    
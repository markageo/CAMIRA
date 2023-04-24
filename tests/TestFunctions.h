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
};

// Read test config files
TestConfig ReadConfig(const std::string &);


// Write sets of data to file
void WriteMesh(const CFD::Mesh &, const std::string &);

// Compare generated mesh with stored one
bool CompareMesh(const CFD::Mesh &, const std::string &);


}   // end namespace TEST


#endif // TEST_FUNCTIONS    
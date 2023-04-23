#ifndef TEST_FUNCTIONS
#define TEST_FUNCTIONS

#include "FiniteVolumeStructures.h"

namespace TEST {

enum testState {
    null,
    write,
    test
};

struct TestConfig
{
    std::string testInputDirectory;
    std::string testInputFilename;

    // Whether to test or write data
    testState meshTest;

};

// Read test config files
TestConfig ReadConfig(const std::string &);

// Write sets of data to file
void WriteMesh(const CFD::Mesh &, const std::string &);

}   // end namespace TEST


#endif // TEST_FUNCTIONS    
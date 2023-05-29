#include "TestFunctions.h"

#include "Types.h"
#include "Utils.h"
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
    std::ifstream &ReadInputLine(std::ifstream &ifs, 
                                 std::string &str) 
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

        } else if (key == "velUniformTestInputFilename") {
            testConfig.velUniformTestInputFilename = value;

        } else if (key == "velZeroGradientTestInputFilename") {
            testConfig.velZeroGradientTestInputFilename = value;

        } else if (key == "velExtrapolatedTestInputFilename") {
            testConfig.velExtrapolatedTestInputFilename = value;

        } else if (key == "faceVelTest") {
            testConfig.faceVelTest = String2TestState(value);

        } else if (key ==  "faceVelTestOutputDirectory") {
            testConfig.faceVelTestOutputDirectory = value;

        } else if (key == "faceVelTestReferenceDirectory") {
            testConfig.faceVelTestReferenceDirectory = value;

        } else if (key == "fvCoeffTest") {
            testConfig.fvCoeffTest = String2TestState(value);

        } else if (key ==  "fvCoeffTestOutputDirectory") {
            testConfig.fvCoeffTestOutputDirectory = value;

        } else if (key == "fvCoeffTestReferenceDirectory") {
            testConfig.fvCoeffTestReferenceDirectory = value;
        }

    }

    return testConfig;
}


void TEST::WriteMesh(const CFD::Mesh &mesh, 
                     const std::string &filedir, 
                     const int precision) 
{

    using AX = CFD::Axis::ENUMDATA;

    AX axis;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"_x", "_y", "_z"} };
    std::filesystem::create_directory(filedir);
    auto Fname = [&] (const std::string &s) -> std::string { return filedir + s + axisSuffix[axis] + ext; };

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
        extrapFactorsArray(0, 0) = mesh.extrapFactors[CFD::PositivePatch[axis]].p;
        extrapFactorsArray(1, 0) = mesh.extrapFactors[CFD::PositivePatch[axis]].a;
        extrapFactorsArray(0, 1) = mesh.extrapFactors[CFD::NegativePatch[axis]].p;
        extrapFactorsArray(1, 1) = mesh.extrapFactors[CFD::NegativePatch[axis]].a;
        UTIL::WriteArray(Fname("extrap_factors"), extrapFactorsArray, precision);

    }

}


void TEST::WriteFields(const CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &fields, 
                       const std::string &filedir, 
                       const int precision) 
{
    using F = CFD::Fields::ENUMDATA;

    std::string ext = TEST::testFileExtension;
    std::filesystem::create_directory(filedir);
    auto Fname = [&] (const std::string &s) -> std::string { return filedir + "fields_" + s + ext; };
    
    UTIL::WriteArray(Fname("u"), fields[F::U], precision);
    UTIL::WriteArray(Fname("v"), fields[F::V], precision);
    UTIL::WriteArray(Fname("w"), fields[F::W], precision);
    UTIL::WriteArray(Fname("p"), fields[F::P], precision);
}


void TEST::WriteFaceVels(const CFD::EnumVector<CFD::BoundaryConditions, CFD::ArrayAllocator<CFD::Fields, CFD::array3D> > &faceVelocities, 
                         const std::string &filedir, 
                         const int precision) 
{
    using F = CFD::Fields::ENUMDATA;
    using BC = CFD::BoundaryConditions::ENUMDATA;

    BC boundaryCondition;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::BoundaryConditions, std::string> boundaryConditionSuffix{ {"zeroGradient", "uniform", "extrapolated"} };
    std::filesystem::create_directory(filedir);
    auto Fname = [&] (const std::string &s) -> std::string { return filedir + "face_velocities_" + boundaryConditionSuffix[boundaryCondition] + "_" + s + ext; };
    
    for (int bc = 0; bc != CFD::BoundaryConditions::count; bc++) {
        boundaryCondition = static_cast<BC>(bc);
        UTIL::WriteArray(Fname("u"), faceVelocities[boundaryCondition][F::U], precision);
        UTIL::WriteArray(Fname("v"), faceVelocities[boundaryCondition][F::V], precision);
        UTIL::WriteArray(Fname("w"), faceVelocities[boundaryCondition][F::W], precision);
    }
}


void TEST::WriteFVCoeffs(const CFD::EnumVector<CFD::BoundaryConditions, CFD::FVCoefficients > &fvCoeffs, 
                         const std::string &filedir, 
                         const int precision) 
{
    using BC = CFD::BoundaryConditions::ENUMDATA;
    using TC = CFD::TransportCoefficients::ENUMDATA;
    using AX =  CFD::Axis::ENUMDATA;

    BC boundaryCondition;
    TC transportCoeff;
    AX axis;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"x", "y", "z"} };
    CFD::EnumVector<CFD::BoundaryConditions, std::string> boundaryConditionSuffix{ {"zeroGradient", "uniform", "extrapolated"} };
    CFD::EnumVector<CFD::TransportCoefficients, std::string> transportCoeffSuffix{ { "tt", "nn", "ee", "t", "n", "e", "p", "w", "s", "b", "ww", "ss", "bb" } };
    std::filesystem::create_directory(filedir);
    auto Fname = [&] (const std::string &eq, const std::string &var) -> std::string { return filedir + "fvCoeffs_" + boundaryConditionSuffix[boundaryCondition] 
                                                                                                     + "_" + eq + "_" + var + "_" + transportCoeffSuffix[transportCoeff] 
                                                                                                     + ext; };

    // Boundary conditions
    for (int bc = 0; bc != CFD::BoundaryConditions::count; bc++) {
        boundaryCondition = static_cast<BC>(bc);

        // Iterate coefficients
        for (int tc = 0; tc != CFD::TransportCoefficients::count; tc++) {
            transportCoeff = static_cast<TC>(tc);

            // Diffusion coefficients
            for (int ax = 0; ax != CFD::Axis::count; ax++) {
                axis = static_cast<AX>(ax);
               
                if ( fvCoeffs[boundaryCondition].Umom.diff[axis].get(transportCoeff) ) {
                    UTIL::WriteArray(Fname("Umom", "diff_" + axisSuffix[axis]), fvCoeffs[boundaryCondition].Umom.diff[axis][transportCoeff], precision);
                }
                if ( fvCoeffs[boundaryCondition].Vmom.diff[axis].get(transportCoeff) ) {
                    UTIL::WriteArray(Fname("Vmom", "diff_" + axisSuffix[axis]), fvCoeffs[boundaryCondition].Vmom.diff[axis][transportCoeff], precision);
                }
                if ( fvCoeffs[boundaryCondition].Wmom.diff[axis].get(transportCoeff) ) {
                    UTIL::WriteArray(Fname("Wmom", "diff_" + axisSuffix[axis]), fvCoeffs[boundaryCondition].Wmom.diff[axis][transportCoeff], precision);
                }

            }


            // Momentum velocity coefficients
            if ( fvCoeffs[boundaryCondition].Umom.AU.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Umom", "AU"), fvCoeffs[boundaryCondition].Umom.AU[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Vmom.AV.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Vmom", "AV"), fvCoeffs[boundaryCondition].Vmom.AV[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Wmom.AW.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Wmom", "AW"), fvCoeffs[boundaryCondition].Wmom.AW[transportCoeff], precision);
            }

            // Momentum pressure coefficients
            if ( fvCoeffs[boundaryCondition].Umom.AP.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Umom", "AP"), fvCoeffs[boundaryCondition].Umom.AP[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Vmom.AP.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Vmom", "AP"), fvCoeffs[boundaryCondition].Vmom.AP[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Wmom.AP.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Wmom", "AP"), fvCoeffs[boundaryCondition].Wmom.AP[transportCoeff], precision);
            }


            // Continuity velocity coefficients
            if ( fvCoeffs[boundaryCondition].Cont.AU.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Cont", "AU"), fvCoeffs[boundaryCondition].Cont.AU[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Cont.AV.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Cont", "AV"), fvCoeffs[boundaryCondition].Cont.AV[transportCoeff], precision);
            }
            if ( fvCoeffs[boundaryCondition].Cont.AW.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Cont", "AW"), fvCoeffs[boundaryCondition].Cont.AW[transportCoeff], precision);
            }

            // Continuity pressure coefficients
            if ( fvCoeffs[boundaryCondition].Cont.AP.get(transportCoeff) ) {
                UTIL::WriteArray(Fname("Cont", "AP"), fvCoeffs[boundaryCondition].Cont.AP[transportCoeff], precision);
            }

        }

    }
}


// Helper functions Comparision tests
namespace 
{

    // Returns true if all elements of a tensor are exactly the same
    template<typename T>
    bool ArraysSame(const T &array1, 
                    const T &array2) 
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
    bool ReadAndCompare(const std::string &filename1, 
                        const std::string &filename2) {

        // Temporaries for storing read in data
        T array1, array2;
        array1 = UTIL::ReadArray<CFD::array1D>(filename1);
        array2 = UTIL::ReadArray<CFD::array1D>(filename2);
        if (!ArraysSame(array1, array2))
            return false;

        return true;

    }


}   // end anonymous namespace


bool TEST::CompareMesh(const std::string &outputDir, 
                       const std::string &referenceDir)
{

    using AX = CFD::Axis::ENUMDATA;

    AX axis;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::Axis, std::string> axisSuffix{ {"x", "y", "z"} };
    bool comparisonResult = true;

    // Temporary array for extrapolation factors
    // rows  : extrapFactors.p, extrapFactors.a
    // cols  : positivePatch, negativepatch
    CFD::array2D extrapFactorsArray(2, 2); 

    // Temporary variable for output and reference directories
    CFD::array1D outputArray1D, referenceArray1D;
    CFD::array2D outputArray2D, referenceArray2D;

    // Lambda for getting to filename
    auto Fname = [&] ( const std::string &name, 
                       const std::string &dir) 
                    -> std::string { 
        return dir + name + "_" + axisSuffix[axis] + ext; 
    };

    // Lambda for performing the check
    auto PerformCheck1D = [&] ( const std::string &meshQuantityName ) {
        outputArray1D = UTIL::ReadArray<CFD::array1D>( Fname( meshQuantityName, outputDir ) );
        referenceArray1D = UTIL::ReadArray<CFD::array1D>( Fname( meshQuantityName, referenceDir ) );
        if ( !ArraysSame( outputArray1D, referenceArray1D ) )
            comparisonResult = false;
    };

    auto PerformCheck2D = [&] ( const std::string &meshQuantityName ) {
        outputArray2D = UTIL::ReadArray<CFD::array2D>( Fname( meshQuantityName, outputDir ) );
        referenceArray2D = UTIL::ReadArray<CFD::array2D>( Fname( meshQuantityName, referenceDir ) );
        if ( !ArraysSame( outputArray1D, referenceArray1D ) )
            comparisonResult = false;
    };
            
    for (int a = 0; a != CFD::Axis::count; a++) {
        axis = static_cast<AX>(a);

        // Check each part of the mesh
        PerformCheck1D( "cell_lengths" );
        PerformCheck1D( "inv_cell_lengths" );
        PerformCheck1D( "cell_centers" );
        PerformCheck1D( "inv_cell_center_diff" );
        PerformCheck1D( "cell_faces" );
        PerformCheck1D( "interp_factors" );
        PerformCheck2D( "extrap_factors" );

    }


    return comparisonResult;
}



bool TEST::CompareFaceVels(const std::string &outputDir, 
                           const std::string &referenceDir)
{
    using BC = CFD::BoundaryConditions::ENUMDATA;

    BC boundaryCondition;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::BoundaryConditions, std::string> boundaryConditionSuffix{ {"zeroGradient", "uniform", "extrapolated"} };
    bool comparisonResult = true;
    
    
    // Temporary variable for output and reference directories
    CFD::array3D outputArray, referenceArray;

    // Lambda for getting filename
    auto Fname = [&] ( const std::string &s, 
                       const std::string &dir) 
                    -> std::string { 
        return dir + "face_velocities" + "_" + boundaryConditionSuffix[boundaryCondition] + "_" + s + ext; 
    };

    // Lambda for performing the check
    auto PerformCheck = [&] ( const std::string &fieldName ) {
        outputArray = UTIL::ReadArray<CFD::array3D>( Fname( fieldName, outputDir ) );
        referenceArray = UTIL::ReadArray<CFD::array3D>( Fname( fieldName, referenceDir ) );
        if (!ArraysSame(outputArray, referenceArray))
            comparisonResult = false;;
    };

    for (int bc = 0; bc != CFD::BoundaryConditions::count; bc++) {
        boundaryCondition = static_cast<BC>(bc);

        // Check each field with each boundary condition
        PerformCheck( "u" );
        PerformCheck( "v" );
        PerformCheck( "w" );

    }

    return comparisonResult;
}




bool TEST::CompareFVCoeffs( const std::string &outputDir, 
                            const std::string &referenceDir )
{

    using BC = CFD::BoundaryConditions::ENUMDATA;
    using enum CFD::TransportCoefficients::ENUMDATA;

    BC boundaryCondition;
    std::string ext = TEST::testFileExtension;
    CFD::EnumVector<CFD::BoundaryConditions, std::string> boundaryConditionSuffix{ {"zeroGradient", "uniform", "extrapolated"} };
    CFD::EnumVector<CFD::TransportCoefficients, std::string> transportCoefficientSuffix{ {"tt", "nn", "ee", "t", "n", "e", "p", "w", "s", "b", "ww", "ss", "bb"} };
    bool comparisonResult = true;
    
    
    // Temporary variable for output and reference directories
    CFD::array3D outputArray3D, referenceArray3D;
    CFD::array1D outputArray1D, referenceArray1D;

    // Lambda for getting filename
    auto Fname = [&] (const std::string &s, 
                      const CFD::TransportCoefficients::ENUMDATA tc, 
                      const std::string &dir) 
                      -> std::string { 
        return dir + "fvCoeffs" + "_" + boundaryConditionSuffix[boundaryCondition] + "_" + s + "_" + transportCoefficientSuffix[tc] + ext; 
    };

    // Lambda for performing the check
    auto PerformCheck1D = [&] ( const std::string &coeffName, 
                              const CFD::TransportCoefficients::ENUMDATA tc) {
        outputArray1D = UTIL::ReadArray<CFD::array1D>( Fname( coeffName, tc, outputDir ) );
        referenceArray1D = UTIL::ReadArray<CFD::array1D>( Fname( coeffName, tc, referenceDir ) );
        if ( !ArraysSame(outputArray1D, referenceArray1D) )
            comparisonResult = false;
    };

    auto PerformCheck3D = [&] ( const std::string &coeffName, 
                              const CFD::TransportCoefficients::ENUMDATA tc) {
        outputArray3D = UTIL::ReadArray<CFD::array3D>( Fname( coeffName, tc, outputDir ) );
        referenceArray3D = UTIL::ReadArray<CFD::array3D>( Fname( coeffName, tc, referenceDir ) );
        if ( !ArraysSame(outputArray3D, referenceArray3D) )
            comparisonResult = false;
    };

    for (int bc = 0; bc != CFD::BoundaryConditions::count; bc++) {
        boundaryCondition = static_cast<BC>(bc);

        // Continuity velocity coefficients
        PerformCheck1D("Cont_AU", p);
        PerformCheck1D("Cont_AU", e);
        PerformCheck1D("Cont_AU", w);

        PerformCheck1D("Cont_AV", p);
        PerformCheck1D("Cont_AV", n);
        PerformCheck1D("Cont_AV", s);

        PerformCheck1D("Cont_AW", p);
        PerformCheck1D("Cont_AW", t);
        PerformCheck1D("Cont_AW", b);

        // Continuity pressure coefficients
        PerformCheck3D("Cont_AP", n);
        PerformCheck3D("Cont_AP", e);
        PerformCheck3D("Cont_AP", s);
        PerformCheck3D("Cont_AP", w);
        PerformCheck3D("Cont_AP", t);
        PerformCheck3D("Cont_AP", b);
        PerformCheck3D("Cont_AP", p);
        PerformCheck3D("Cont_AP", nn);
        PerformCheck3D("Cont_AP", ee);
        PerformCheck3D("Cont_AP", ss);
        PerformCheck3D("Cont_AP", ww);
        PerformCheck3D("Cont_AP", tt);
        PerformCheck3D("Cont_AP", bb);


        // U momentum velocity
        PerformCheck3D("Umom_AU", p);
        PerformCheck3D("Umom_AU", n);
        PerformCheck3D("Umom_AU", e);
        PerformCheck3D("Umom_AU", s);
        PerformCheck3D("Umom_AU", w);
        PerformCheck3D("Umom_AU", t);
        PerformCheck3D("Umom_AU", b);

        // U momentum presuure
        PerformCheck1D("Umom_AP", p);
        PerformCheck1D("Umom_AP", e);
        PerformCheck1D("Umom_AP", w);


        // V momentum velocity
        PerformCheck3D("Vmom_AV", p);
        PerformCheck3D("Vmom_AV", n);
        PerformCheck3D("Vmom_AV", e);
        PerformCheck3D("Vmom_AV", s);
        PerformCheck3D("Vmom_AV", w);
        PerformCheck3D("Vmom_AV", t);
        PerformCheck3D("Vmom_AV", b);

        // V momentum presuure
        PerformCheck1D("Vmom_AP", p);
        PerformCheck1D("Vmom_AP", n);
        PerformCheck1D("Vmom_AP", s);


        // W momentum velocity
        PerformCheck3D("Wmom_AW", p);
        PerformCheck3D("Wmom_AW", n);
        PerformCheck3D("Wmom_AW", e);
        PerformCheck3D("Wmom_AW", s);
        PerformCheck3D("Wmom_AW", w);
        PerformCheck3D("Wmom_AW", t);
        PerformCheck3D("Wmom_AW", b);

        // W momentum presuure
        PerformCheck1D("Wmom_AP", p);
        PerformCheck1D("Wmom_AP", t);
        PerformCheck1D("Wmom_AP", b);

    }

    return comparisonResult;

}
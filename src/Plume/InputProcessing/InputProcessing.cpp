#include "InputProcessing.h"

#include "Core/IO/InputParser.h"
#include "Core/IO/CSVReader.h"
#include "Core/IO/IOTools.h"
#include "Core/Types.h"
#include "Core/FVLookups.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "Core/IO/ptreeTranslators.h"

#include <Eigen/Geometry>

#include <utility>
#include <optional>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <filesystem>
#include <memory>


namespace pt = boost::property_tree;

namespace CAMIRA
{

using namespace CORE;
namespace PLUME
{

namespace
{

    /*-------------------------------------------------------------------------------------*\
                                               Model
    \*-------------------------------------------------------------------------------------*/

    void ReadModel( InputData &inputData, 
                    const pt::ptree &tree )
    {
        const pt::ptree &modelTree = tree.get_child("Model");

        inputData.diffusionCoeff               = modelTree.get<floatType>("D");
        inputData.numberOfTimeSteps            = modelTree.get<intType>("numberOfTimeSteps");
        inputData.initialParticlesPerUnitMass  = modelTree.get<floatType>("initialParticlesPerUnitMass");
        inputData.timeStepSize                 = modelTree.get<floatType>("timeStepSize");
    }


    /*-------------------------------------------------------------------------------------*\
                                            VelocityField
    \*-------------------------------------------------------------------------------------*/

    void ReadVelocityField( InputData &inputData, 
                            const pt::ptree &tree )
    {
        const pt::ptree &velocityFieldTree = tree.get_child("VelocityField");

        inputData.velocityFieldFilename  = velocityFieldTree.get<std::string>("velocityFieldFilename");
        IOTOOLS::PrependRelativePath( inputData.velocityFieldFilename, inputData.inputFileDirectory );
    }



    /*-------------------------------------------------------------------------------------*\
                                             Sources
    \*-------------------------------------------------------------------------------------*/

    void ReadSources( InputData &inputData, 
                      const pt::ptree &tree )
    {
        const pt::ptree &sourcesTree = tree.get_child("Sources");

        for (auto source : sourcesTree ) {

            if        ( source.first == "ContinuousReleasePoint" ) {

                InputData::ContinuousReleasePointData tempContinuousReleasePointData;
                std::vector<floatType> tempLocation = source.second.get< std::vector<floatType> >( "location" ); 
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    tempContinuousReleasePointData.location(axis) = tempLocation[axis];
                } );

                tempContinuousReleasePointData.massFlowRate = source.second.get<floatType>( "massFlowRate" );

                inputData.continuousReleasePoints.push_back( tempContinuousReleasePointData );

            } else if (source.first == "InstantaneousReleasePoint") {

                InputData::InstantaneousReleasePointData tempInstantaneousReleasePointData;
                std::vector<floatType> tempLocation = source.second.get< std::vector<floatType> >( "location" ); 
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    tempInstantaneousReleasePointData.location(axis) = tempLocation[axis];
                } );

                tempInstantaneousReleasePointData.totalMass = source.second.get<floatType>( "totalMass" );

                inputData.instantaneousReleasePoints.push_back( tempInstantaneousReleasePointData );

            } else {
                throw std::runtime_error(  "'" + source.first + "' is not a valid Sources child name" );
            }

        }
        
    }


    /*-------------------------------------------------------------------------------------*\
                                       Boundary Conditions
    \*-------------------------------------------------------------------------------------*/

    InputData::BoundaryConditionInputData ReadBoundaryValueString( const std::string &fieldString,
                                                                   const pt::ptree &boundaryPatchTree,
                                                                   [[maybe_unused]] const BoundaryPatches::ENUMDATA boundaryPatch,
                                                                   [[maybe_unused]] const std::string &inputFileDirectory )
    {
        using BC = BoundaryConditions;

        // Get the boundary condition string from the tree string
        std::string boundaryString = boundaryPatchTree.get<std::string>( fieldString );

        InputData::BoundaryConditionInputData bcStruct;
        std::string bcTypeString;
        std::string::const_iterator stringIterator = boundaryString.begin();

        // Read the boundary condition type
        while (stringIterator != boundaryString.end()) {

            if (*stringIterator == MULTI_DELIMITER_CHAR) 
                break;
            
            bcTypeString += *stringIterator;
            stringIterator++;

        }


        // Store the type, some don't need a value
        if        (bcTypeString == "outflow") {

            bcStruct.type = BC::outflow;

        } else if (bcTypeString == "reflection") {

            bcStruct.type = BC::reflection;

        } else if (bcTypeString == "periodic") {

            bcStruct.type = BC::periodic;

        } else {

            throw std::runtime_error( "'" + bcTypeString + "' is not a valid boundary condition type." );
        }

        return bcStruct;
    }


    void ReadBoundaryConditions(InputData &inputData, 
                                const pt::ptree &tree)
    {
        using BP = BoundaryPatches::ENUMDATA;
        using enum Axis::ENUMDATA;

        const pt::ptree &boundaryConditionsTree = tree.get_child( "BoundaryConditions" );

        // Boundary name strings
        std::map<BP, std::string> boundaryPatchMap{ {BP::xPositive, "+x"}, {BP::xNegative, "-x"},
                                                    {BP::yPositive, "+y"}, {BP::yNegative, "-y"},
                                                    {BP::zPositive, "+z"}, {BP::zNegative, "-z"} };

        // Iterate boundary patches
        std::string valueString;
        const pt::ptree *boundaryPatchTreePointer = nullptr;
        for (const auto& [patchEnum, patchString] : boundaryPatchMap) {
            boundaryPatchTreePointer = &( boundaryConditionsTree.get_child( patchString ) );

            inputData.boundaryConditions[patchEnum] = ReadBoundaryValueString( "phi", *boundaryPatchTreePointer, patchEnum, inputData.inputFileDirectory );

        }


        // Verify that periodic boundary conditions are set for both sides
        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

            if ( inputData.boundaryConditions[bp].type == BoundaryConditions::periodic ) {

                Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ bp ];
                BoundaryPatches::ENUMDATA bpOpposite = ( bp == LUT::PositivePatch[axis] ) ? LUT::NegativePatch[axis] : LUT::PositivePatch[axis];
                if ( inputData.boundaryConditions[bpOpposite].type != BoundaryConditions::periodic ) {
                        throw std::runtime_error( "Specification of periodic boundary conditions requires both sides of domain to be specified as periodic." );
                }

            }

        } );

    }


    /*-------------------------------------------------------------------------------------*\
                                               Parallel
    \*-------------------------------------------------------------------------------------*/

    void ReadParallel( InputData &inputData, 
                       const pt::ptree &tree)
    {
        using enum Axis::ENUMDATA;

        const pt::ptree &parallelTree = tree.get_child("Parallel");

        inputData.parallelSettings.numberOfThreads = parallelTree.get<intType>("numberOfThreads");
    }


    /*-------------------------------------------------------------------------------------*\
                                            Output
    \*-------------------------------------------------------------------------------------*/

    void ReadMonitors( InputData &inputData,
                       const pt::ptree &outputTree )
    {

        // It's ok if there is no monitors tree
        boost::optional<const pt::ptree &> monitorsTreeOptional = outputTree.get_child_optional( "Monitors" );
        if ( !monitorsTreeOptional )
            return;

        const pt::ptree &monitorsTree = monitorsTreeOptional.get();

        for (auto monitor : monitorsTree) {

            if ( monitor.first == "Probe" ) {
                InputData::ProbeData tempProbeData;
                std::vector<floatType> tempLocation = monitor.second.get< std::vector<floatType> >( "location" ); 
                tempProbeData.filename              = monitor.second.get<std::string>( "filename" );
                IOTOOLS::PrependRelativePath( tempProbeData.filename, inputData.inputFileDirectory );
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    tempProbeData.location(axis) = tempLocation[axis];
                } );

                inputData.probes.push_back( tempProbeData );

            } else {
                throw std::runtime_error(  "'" + monitor.first + "' is not a valid Monitors child name" );
            }
            
        }
    }


    void VerifyWriteDirectory( const std::string &fileDir ) 
    {
        std::filesystem::path filepath = std::string(fileDir);
        if ( filepath.parent_path().empty() ) { // The user can specify the exectuable path
            return;
        }

        bool filePathExists = std::filesystem::exists(filepath.parent_path());
        if ( !filePathExists ) {
            throw std::runtime_error("Cannot find or access directory '" + std::string(filepath.parent_path()) + "'. Please Make sure it exists.");
        }
    }


    void VerifyOutputFiles( InputData &inputData ) 
    {
        VerifyWriteDirectory( inputData.fieldOutputFilename );
        for ( const auto &probe : inputData.probes ) {
            VerifyWriteDirectory( probe.filename );
        }
    }


    void ReadOutput( InputData &inputData,
                    const pt::ptree &tree )
    {
        const pt::ptree &outputTree = tree.get_child( "Output" );

        // Specifying field output format is optional, default to Binary
        inputData.outputFormatType = InputData::OutputFormatType::BINARY;
        boost::optional<std::string> outputFormatTypeOptional = outputTree.get_optional<std::string>( "outputFormatType" );
        if ( outputFormatTypeOptional ) {
            std::string outputFormatType = outputFormatTypeOptional.get();
            std::transform( outputFormatType.begin(), outputFormatType.end(), outputFormatType.begin(), ::tolower );
            if ( outputFormatType == "binary" ) {
                inputData.outputFormatType = InputData::OutputFormatType::BINARY;
            } else if ( outputFormatType  == "ascii") {
                inputData.outputFormatType = InputData::OutputFormatType::ASCII;
            } else {
                throw std::runtime_error(  "'" + outputFormatType + "' is not a valid output format type!" );
            }
        }


        // Profiling information filename
        inputData.profilingFilename = outputTree.get<std::string>( "profilingFilename" );
        IOTOOLS::PrependRelativePath( inputData.profilingFilename, inputData.inputFileDirectory );

        // Field output filename
        inputData.fieldOutputFilename = outputTree.get<std::string>( "fieldOutputFilename" );
        IOTOOLS::PrependRelativePath( inputData.fieldOutputFilename, inputData.inputFileDirectory );
        
        // Write interval for fields
        inputData.fieldWriteInterval = outputTree.get<intType>( "fieldWriteInterval" );

        // Monitors
        ReadMonitors( inputData, outputTree );

        // Verify the output filenames
        VerifyOutputFiles( inputData ); 
    }

}   // end anonymous namepsace




// Parse input file and read into InputData structure
InputData ReadInputData(const std::string &inputFilename) 
{
    pt::ptree tree = INP::ParseFile(inputFilename);

    InputData inputData;
    inputData.inputFileDirectory = IOTOOLS::RelativePath( inputFilename );

    ReadModel(inputData, tree);
    ReadVelocityField(inputData, tree);
    ReadSources(inputData, tree);
    ReadBoundaryConditions(inputData, tree);
    ReadParallel(inputData, tree);
    ReadOutput(inputData, tree);

    return inputData;
}



// Read command line input for input file name
InputData InputDataFromCommandLine(int argc, char const *argv[])
{
    // InputData object to return
    InputData inputData;

    // Input file from first command line argument
    std::string inputFilename;
    if (argc == 1) {
        std::cout << "Please enter input file name: "
                  << "\n";
        std::cin >> inputFilename;
        std::cout << "\n";
        std::cin.ignore();

    }
    else if (argc == 2) {
        inputFilename = argv[1];

    }
    else {
        throw std::invalid_argument("Invalid command line options.");

    }

    // User input data
    std::string inputFileRetryChoice;
    while (true)
    {
        try
        {
            std::cout << "Reading input file '" + inputFilename + "' ... ";
            inputData = ReadInputData(inputFilename);
            std::cout << "Success."
                      << "\n\n";
            break;
        } 
        catch (std::runtime_error &e) 
        {
            std::cout << "Failed. \n"
                      << e.what()
                      << "\n\n";

            std::cout << "Would you like to try again? (y/n)"
                      << "\n";
            std::cin >> inputFileRetryChoice;
            if (inputFileRetryChoice != "y")
                exit(-1);
            std::cout << "\n";

            std::cout << "Please enter input file name: "
                      << "\n";
            std::cin >> inputFilename;
            std::cout << "\n";
            std::cin.ignore();
        }
    }
    
    return inputData;
}



}   // end namespace PLUME

}   // end namespace CAMIRA

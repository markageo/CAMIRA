#include "InputParser.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Types.h"
#include "CSVReader.h"

#include "Boost/boost/property_tree/ptree.hpp"
#include <Eigen/Geometry>

#include <utility>
#include <optional>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>

#define VECTOR_START_CHAR       '('
#define VECTOR_END_CHAR         ')'
#define VECTOR_DELIMITER_CHAR   ','
#define MULTI_DELIMITER_CHAR    ','

namespace pt = boost::property_tree;


// Read command line input for input file name
CFD::InputData CFD::InputDataFromCommandLine(int argc, char const *argv[])
{
    // InputData object to return
    CFD::InputData inputData;

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
            std::cout << "Reading input file: '" + inputFilename + "' ..."
                      << "\n\n";
            inputData = CFD::ReadInputData(inputFilename);
            std::cout << "Success!"
                      << "\n\n";
            break;
        } 
        catch (std::runtime_error &e) 
        {
            std::cout << "Failure reading input file! \n"
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
    std::cout << "Press enter to begin.";
    std::cin.ignore();
    std::cout << std::endl;

    return inputData;
}


// InputData constructor
CFD::InputData::InputData() :
    meshSegments(),
    boundaryConditions()
    {
        linearSolverSettings.planeSweepDirection = CFD::BoundaryPatches::zPositive;
        linearSolverSettings.lineSweepDirection = CFD::BoundaryPatches::yPositive;
    };


/*-------------------------------------------------------------------------------------*\
                                      Translators
\*-------------------------------------------------------------------------------------*/

// Convert string to given numeric type T.
template <typename T> T 
String2Type(const std::string &str)
{
    // NOTE: This does not work for ints in scientific notation.
    std::istringstream strstream(str);
    T num;
    strstream >> num;
    return num;
}


// Parse vector string into an std::vector
template<typename T>
std::vector<T> ParseVectorString( const std::string &vecString )
{
    std::vector<T> vec;
    std::string::const_iterator stringIterator = vecString.begin();
    std::string valueString;

    if (*stringIterator != VECTOR_START_CHAR) {
        throw std::runtime_error( "'" + vecString + "' is not a valid input. Expected a vector." );
    }
    ++stringIterator;
    
    while ( stringIterator != vecString.end() ) {
        if ( *stringIterator == VECTOR_END_CHAR ) {
            vec.push_back( String2Type<T>(valueString) );
            break;
        }

        if ( *stringIterator == VECTOR_DELIMITER_CHAR ) {
            vec.push_back( String2Type<T>(valueString) );
            valueString.clear();
        } else {
            valueString += *stringIterator;
        }
        stringIterator++;
    }

    if (*stringIterator != VECTOR_END_CHAR) {
        throw std::runtime_error( "'" + vecString + "' vector not closed." );
    }

    return vec;
}


template< typename T >
struct VectorTranslator
{
    typedef std::string     internal_type;
    typedef std::vector<T>  external_type;

    external_type get_value( internal_type const &s ) {
        return ParseVectorString<T>( s );
    }
};


 // Specialization allows the translator to be used with ptree internally
namespace boost { namespace property_tree {

    template< typename T > 
    struct translator_between< std::string, std::vector< T > > 
    {
        typedef VectorTranslator< T > type;
    };

}   // end namespace property_tree  
}   // end namespace boost



namespace
{

    using namespace CFD;

    /*-------------------------------------------------------------------------------------*\
                                               Model
    \*-------------------------------------------------------------------------------------*/

    void ReadModel(InputData &inputData, 
                   const pt::ptree &tree)
    {
        using enum Axis::ENUMDATA;

        const pt::ptree &modelTree = tree.get_child("Model");

        inputData.rho = modelTree.get<floatType>("rho");
        inputData.nu  = modelTree.get<floatType>("nu");        
    }



    /*-------------------------------------------------------------------------------------*\
                                                Mesh
    \*-------------------------------------------------------------------------------------*/

    void ValidateSegmentBounds(const std::vector<InputData::MeshSegment> &meshSegment, 
                               const floatType &domainSize)
    {

        // The first one must start at zero
        if (meshSegment.front().lowerBound != 0.0) {
            throw std::runtime_error( "Mesh segments must start at coordinate 0." );
        }

        // Must end at the given domain size
        if (meshSegment.back().upperBound != domainSize ) {
            throw std::runtime_error( "Mesh segments to not match given domain bounds." );
        }


        // Must be no gaps or overlaps in the segments
        for (size_t i = 1; i != meshSegment.size(); i++) {
            if ( meshSegment[i].lowerBound != meshSegment[i-1].upperBound ) {
                throw std::runtime_error( "Mesh segments must not overlap or have gaps." );
            }
        }

    }


    void ReadGrid( std::vector<InputData::MeshSegment> &meshSegments,
                   const pt::ptree &meshTree,   
                   const std::string &gridString)
    {
        const pt::ptree &gridTree = meshTree.get_child(gridString);
        std::string boundsString, nCellsString, biasFactorString;
        std::vector<floatType> tempBoundsVector;
        InputData::MeshSegment tempMeshSegment;
        for (auto segment : gridTree) {
            if (segment.first != "Segment") {
                throw std::runtime_error(  "'" + segment.first + "' is not a valid grid child name" );
            }
            
            tempMeshSegment.nCells     = segment.second.get<intType>( "nCells" );
            tempMeshSegment.biasFactor = segment.second.get<floatType>( "biasFactor" );
            tempBoundsVector           = segment.second.get< std::vector<floatType> >( "bounds" );
            tempMeshSegment.lowerBound = tempBoundsVector[0];
            tempMeshSegment.upperBound = tempBoundsVector[1];

            meshSegments.push_back(tempMeshSegment);
        }
    }


    void ReadMesh(InputData &inputData, 
                  const pt::ptree &tree)
    {
        using enum Axis::ENUMDATA;

        const pt::ptree &meshTree = tree.get_child( "Mesh" );

        // Domain
        std::vector<floatType> domainSizeTemp = meshTree.get< std::vector<floatType> >( "domain" );
        inputData.domainSize(0) = domainSizeTemp[0];
        inputData.domainSize(1) = domainSizeTemp[1];
        inputData.domainSize(2) = domainSizeTemp[2];

        // Grids
        ReadGrid(inputData.meshSegments[X], meshTree, "Gridx");
        ReadGrid(inputData.meshSegments[Y], meshTree, "Gridy");
        ReadGrid(inputData.meshSegments[Z], meshTree, "Gridz");

        // Sort in increasing order of lower bound
        auto sortComparison = [](const auto& i, const auto& j) { return i.lowerBound < j.lowerBound; };
        std::sort( inputData.meshSegments[X].begin(), inputData.meshSegments[X].end(), sortComparison);
        std::sort( inputData.meshSegments[Y].begin(), inputData.meshSegments[Y].end(), sortComparison );
        std::sort( inputData.meshSegments[Z].begin(), inputData.meshSegments[Z].end(), sortComparison );

        // Check that the bounds are valid
        ValidateSegmentBounds(inputData.meshSegments[X], inputData.domainSize[X]);
        ValidateSegmentBounds(inputData.meshSegments[Y], inputData.domainSize[Y]);
        ValidateSegmentBounds(inputData.meshSegments[Z], inputData.domainSize[Z]);

    }


    /*-------------------------------------------------------------------------------------*\
                                       Boundary Conditions
    \*-------------------------------------------------------------------------------------*/

    // Checks the line that a valid value for the boundary condition has been given
    void CheckForBCValue( std::string::const_iterator &stringIterator,
                               const std::string &bcTypeString )
    {
        if (*stringIterator != MULTI_DELIMITER_CHAR) {
            throw std::runtime_error( "Expected value for " + bcTypeString + " boundary condition."  );
        }
        stringIterator++;
    }



    // Read uniform BC value. Assumes that the string iterator is after the comma
    floatType GetUniformBCValue( std::string::const_iterator &stringIterator,
                                 const std::string &boundaryString )
    {

        // Read the value to end of line
        std::string bcValueString;
        while (stringIterator != boundaryString.end()) {            
            bcValueString += *stringIterator;
            stringIterator++;
        }

        return String2Type<floatType>(bcValueString);
    }


    // Read uniform BC value. Assumes that the string iterator is after the comma
    InputData::Profile1D GetProfile1DValue( std::string::const_iterator &stringIterator,
                                            const std::string &boundaryString )
    {
        InputData::Profile1D profile1D;

        // Read the filename, which is after the comma and up to the end of the line
        std::string filename;
        for ( /* NULL */; stringIterator != boundaryString.end(); stringIterator++ ) { 
            filename += *stringIterator;
        }

        // Read in profile data from csv file
        std::vector< std::vector< std::string > > profileData = ReadCSV( filename );

        // First column tells us the axis
        if        ( profileData[0][0] == "x" ) {
            profile1D.axis = Axis::X;
        } else if ( profileData[0][0] == "y" ) {
            profile1D.axis = Axis::Y;
        } else if ( profileData[0][0] == "z" ) {
            profile1D.axis = Axis::Z;
        } else {
            throw std::runtime_error( "Invalid axis name in profile file '" + filename + "'."  );
        }

        // First column is the coordinate points, second column is the actual data
        intType nRows = profileData.size();
        intType nHeaderRows = 1;
        profile1D.coordinates = array1D( nRows - nHeaderRows );
        profile1D.values      = array1D( nRows - nHeaderRows );

        for ( intType i = 0; i != nRows-nHeaderRows; i++ ) {
            profile1D.coordinates(i) = String2Type<floatType>( profileData[i+nHeaderRows][0] );
            profile1D.values(i)      = String2Type<floatType>( profileData[i+nHeaderRows][1] );
        }

        return profile1D;
    }



    InputData::BoundaryConditionData ReadBoundaryValueString( const std::string &fieldString,
                                                              const pt::ptree &boundaryPatchTree )
    {
        using BC = BoundaryConditions::ENUMDATA;

        // Get the boundary condition string from the tree string
        std::string boundaryString = boundaryPatchTree.get<std::string>( fieldString );

        InputData::BoundaryConditionData bcStruct;
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
        if        (bcTypeString == "uniform") {

            CheckForBCValue( stringIterator, bcTypeString );
            bcStruct.uniformValue = GetUniformBCValue( stringIterator, boundaryString );
            bcStruct.hasUniformValue = true;
            bcStruct.type = BC::fixed;

        } else if (bcTypeString == "profile1D") {

            CheckForBCValue( stringIterator, bcTypeString );
            bcStruct.profile1D = GetProfile1DValue( stringIterator, boundaryString );
            bcStruct.hasProfile1D = true;
            bcStruct.type = BC::fixed;

        } else if (bcTypeString == "zeroGradient") {

            bcStruct.type = BC::zeroGradient;
            bcStruct.uniformValue = 0.0;
            return bcStruct;

        } else if (bcTypeString == "extrapolated") {

            bcStruct.type = BC::extrapolated;
            bcStruct.uniformValue = 0.0;
            return bcStruct;

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

            // Momentum 
            inputData.boundaryConditions.U[X][patchEnum] = ReadBoundaryValueString( "u", *boundaryPatchTreePointer );
            inputData.boundaryConditions.U[Y][patchEnum] = ReadBoundaryValueString( "v", *boundaryPatchTreePointer );
            inputData.boundaryConditions.U[Z][patchEnum] = ReadBoundaryValueString( "w", *boundaryPatchTreePointer );

            // Pressure
            inputData.boundaryConditions.P[patchEnum] = ReadBoundaryValueString( "p", *boundaryPatchTreePointer );

        }

    }


    /*-------------------------------------------------------------------------------------*\
                                           Solver
    \*-------------------------------------------------------------------------------------*/


    CFD::BoundaryPatches::ENUMDATA ReadAxisDirection(const std::string &axisString)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        if        ( axisString == "+x" ) {
            return BP::xPositive;

        } else if ( axisString == "-x" ) {
            return BP::xNegative;

        } else if ( axisString == "+y" ) {
            return BP::yPositive;
        
        } else if ( axisString == "-y" ) {
            return BP::yNegative;

        } else if ( axisString == "+z" ) {
            return BP::zPositive;

        } else if ( axisString == "-z" ) {
            return BP::zNegative;

        } else {
            // throw ERROR - invalid sweeping direction
            return BP::xPositive;
        }
    }


    void ReadSchemes( InputData &inputData, 
                      const pt::ptree &solverTree) 
    {
        const pt::ptree &schemesTree = solverTree.get_child( "Schemes" );
        std::string valueString;

        // Linearisation
        valueString = schemesTree.get<std::string>( "linearisation" );
        if        ( valueString == "Picard" ) {
            inputData.schemes.linearisation = Linearisation::Picard;
        } else if ( valueString == "Newton" ) {
            inputData.schemes.linearisation = Linearisation::Newton;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a valid linearisation." );
        }

        // Momentum interpolation
        valueString = schemesTree.get<std::string>( "momentumInterpolation" );
        if        ( valueString == "implicit" ) {
            inputData.schemes.momentumInterpolation = MomentumInterpolation::Implicit;
        } else if ( valueString == "semiExplicit" ) {
            inputData.schemes.momentumInterpolation = MomentumInterpolation::SemiExplicit;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a valid momentum interpolation method." );
        }

        // Advection scheme
        valueString = schemesTree.get<std::string>( "advectionScheme" );
        if        ( valueString == "upwind" ) {
            inputData.schemes.advectionScheme = AdvectionSchemes::Upwind;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a valid advection scheme." );
        }

        // Face interpolation scheme
        valueString = schemesTree.get<std::string>( "faceInterpolationScheme" );
        if        ( valueString == "weightedLinear" ) {
            inputData.schemes.faceInterpolationScheme = FaceInterpolationSchemes::WeightedLinear;
        } else if ( valueString == "average" ) {
            inputData.schemes.faceInterpolationScheme = FaceInterpolationSchemes::Average;   
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a valid face interpolation scheme." );
        }


        // Momentum implicit relaxation
        std::vector<floatType> momentumRelaxation = schemesTree.get< std::vector<floatType> >( "implicitMomentumRelaxation" );
        inputData.schemes.implicitRelaxation.U[Axis::X] = momentumRelaxation[0];
        inputData.schemes.implicitRelaxation.U[Axis::Y] = momentumRelaxation[1];
        inputData.schemes.implicitRelaxation.U[Axis::Z] = momentumRelaxation[2];

        // Pressure implicit relaxation
        inputData.schemes.implicitRelaxation.P = schemesTree.get<floatType>( "implicitPressureRelaxation" );

        // Max outer iterations
        inputData.schemes.maxOuterIterations = schemesTree.get<intType>( "maxOuterIterations" );

        // Max outer residuals
        floatType maxOuterResiduals = schemesTree.get<floatType>( "maxOuterResiduals" );
        inputData.schemes.maxOuterResiduals.U = maxOuterResiduals;
        inputData.schemes.maxOuterResiduals.P = maxOuterResiduals;
    }



    void ReadLinearSolverSettings( InputData &inputData, 
                                   const pt::ptree & solverTree) 
    {
        const pt::ptree &linearSolverTree = solverTree.get_child( "LinearSolver" );
        std::string valueString;

        // Solver type
        valueString = linearSolverTree.get<std::string>( "type" );
        if        ( valueString == "SUGS" ) {
            inputData.linearSolverSettings.type = LinearSolvers::SUGS;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a linear solver type." );
        }

        // Max iterations
        inputData.linearSolverSettings.maxIterations = linearSolverTree.get<intType>( "maxIterations" );

        // Max residuals
        floatType maxResiduals = linearSolverTree.get<floatType>( "maxResiduals" );
        inputData.linearSolverSettings.maxResiduals.U = maxResiduals;
        inputData.linearSolverSettings.maxResiduals.P = maxResiduals;

        // Momentum relaxation
        std::vector<floatType> momentumRelaxation = linearSolverTree.get< std::vector<floatType> >( "momentumRelaxation" );
        inputData.linearSolverSettings.relaxation.U[Axis::X] = momentumRelaxation[0];
        inputData.linearSolverSettings.relaxation.U[Axis::Y] = momentumRelaxation[1];
        inputData.linearSolverSettings.relaxation.U[Axis::Z] = momentumRelaxation[2];

        // Pressure relaxation
        inputData.linearSolverSettings.relaxation.P = linearSolverTree.get<floatType>( "pressureRelaxation" );

        // Plane sweep direction
        valueString = linearSolverTree.get<std::string>( "planeSweepDirection" );
        inputData.linearSolverSettings.planeSweepDirection = ReadAxisDirection( valueString );

        // Line sweep direction
        valueString = linearSolverTree.get<std::string>( "lineSweepDirection" );
        inputData.linearSolverSettings.lineSweepDirection = ReadAxisDirection( valueString );

    }



    void ReadSolver(InputData &inputData, 
                    const pt::ptree &tree)
    {
        const pt::ptree &solverTree = tree.get_child("Solver");

        // Read discretisation schemes
        ReadSchemes(inputData, solverTree);

        // Plane sweep settings
        ReadLinearSolverSettings(inputData, solverTree);

    }


    /*-------------------------------------------------------------------------------------*\
                                     Initial Conditions
    \*-------------------------------------------------------------------------------------*/

    void ReadInitialConditions( InputData &inputData,
                                const pt::ptree &tree )
    {
        using enum Axis::ENUMDATA;
        const pt::ptree &initialConditionsTree = tree.get_child("InitialConditions");

        inputData.initialConditions.U[X] = initialConditionsTree.get<floatType>( "u" );
        inputData.initialConditions.U[Y] = initialConditionsTree.get<floatType>( "v" );
        inputData.initialConditions.U[Z] = initialConditionsTree.get<floatType>( "w" );

        inputData.initialConditions.P    = initialConditionsTree.get<floatType>( "p" );
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

        InputData::ProbeData tempProbeData;
        std::vector<floatType> tempLocation;
        for (auto probe : monitorsTree) {
            if (probe.first != "Probe") {
                throw std::runtime_error(  "'" + probe.first + "' is not a valid Monitors child name" );
            }

            tempProbeData.filename = probe.second.get<std::string>( "filename" );
            tempLocation           = probe.second.get< std::vector<floatType> >( "location" ); 
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                tempProbeData.location(axis) = tempLocation[axis];
            } );

            inputData.probes.push_back( tempProbeData );
        }
    }

    void VerifyOutputFiles( InputData &inputData ) 
    {
        // Fields
        std::ofstream fileStream( inputData.fieldOutputFilename );
        if ( !fileStream )
            throw std::runtime_error("File '" + inputData.fieldOutputFilename + "' cannot be written to or accessed.");

        // Residual history
        fileStream.close();
        fileStream.open( inputData.residualHistoryFilename );
        if ( !fileStream )
            throw std::runtime_error("File '" + inputData.residualHistoryFilename + "' cannot be written to or accessed.");

        // Probes
        for ( const auto &probe : inputData.probes ) {
            fileStream.close();
            fileStream.open( probe.filename );
            if ( !fileStream )
                throw std::runtime_error("File '" + probe.filename + "' cannot be written to or accessed.");
        }

    }


    void ReadOutput( InputData &inputData,
                    const pt::ptree &tree )
    {
        const pt::ptree &outputTree = tree.get_child( "Output" );

        // Residual history filename
        inputData.residualHistoryFilename = outputTree.get<std::string>( "residualHistoryFilename" );

        // Field output filename
        inputData.fieldOutputFilename = outputTree.get<std::string>( "fieldOutputFilename" );

        // Probes
        ReadMonitors( inputData, outputTree );

        // Verify the output filenames
        VerifyOutputFiles( inputData ); 
    }

}   // end anonymous namepsace


// Parse input file and read into InputData structure
CFD::InputData CFD::ReadInputData(const std::string &inputFileName) 
{
    pt::ptree tree = ParseFile(inputFileName);
    
    InputData inputData;
    ReadModel(inputData, tree);
    ReadMesh(inputData, tree);
    ReadBoundaryConditions(inputData, tree);
    ReadInitialConditions(inputData, tree);
    ReadSolver(inputData, tree);
    ReadOutput(inputData, tree);

    return inputData;
}

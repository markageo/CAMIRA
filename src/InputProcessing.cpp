#include "InputProcessing.h"

#include "InputParser/InputParser.h"
#include "FiniteVolume.h"
#include "Types.h"

#include "Boost/boost/property_tree/ptree.hpp"
#include "Eigen/Geometry"

#include <utility>
#include <optional>
#include <iostream>
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
        planeSolverSettings.sweepDirection = CFD::BoundaryPatches::zPositive;
        lineSolverSettings.sweepDirection = CFD::BoundaryPatches::yPositive;
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
        // throw ERROR - expecting vector
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
        // throw ERROR - vector not closed
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
            // throw ERROR - must have mesh segment starting at coordinate 0
        }

        // Must end at the given domain size
        if (meshSegment.back().upperBound != domainSize ) {
            // throw ERROR - segments dont match given domain bounds
        }


        // Must be no gaps or overlaps in the segments
        for (size_t i = 1; i != meshSegment.size(); i++) {
            if ( meshSegment[i].lowerBound != meshSegment[i-1].upperBound ) {
                // throw ERROR - segments must not overlap or have gaps
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
                // throw ERROR invalid child name
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

    InputData::BoundaryConditionStruct ReadBoundaryValueString(const std::string &boundaryString)
    {
        using BC = BoundaryConditions::ENUMDATA;

        InputData::BoundaryConditionStruct bcStruct;
        std::string bcTypeString;
        std::string bcValueString;
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

            bcStruct.type = BC::uniform;

        } else if (bcTypeString == "zeroGradient") {

            bcStruct.type = BC::zeroGradient;
            bcStruct.value = 0.0;
            return bcStruct;

        } else if (bcTypeString == "extrapolated") {
            bcStruct.type = BC::extrapolated;
            bcStruct.value = 0.0;
            return bcStruct;
        } else {
            // throw ERROR -  invalid BC type
        }


        // Check for expected value
        if (*stringIterator != MULTI_DELIMITER_CHAR) {
            // throw ERROR - expected value for boudnary condition
        }
        stringIterator++;

        // Read the value to end of line
        while (stringIterator != boundaryString.end()) {            
            bcValueString += *stringIterator;
            stringIterator++;
        }

        // Convert and store value
        if (bcStruct.type == BC::uniform) {
            bcStruct.value = std::stod(bcValueString);
        }

        return bcStruct;
    }


    void ReadBoundaryConditions(InputData &inputData, 
                                const pt::ptree &tree)
    {
        using BP = BoundaryPatches::ENUMDATA;
        using F  = Fields::ENUMDATA;

        const pt::ptree &boundaryConditionsTree = tree.get_child( "BoundaryConditions" );

        // Boundary name strings
        std::map<BP, std::string> boundaryPatchMap{ {BP::xPositive, "+x"}, {BP::xNegative, "-x"},
                                                    {BP::yPositive, "+y"}, {BP::yNegative, "-y"},
                                                    {BP::zPositive, "+z"}, {BP::zNegative, "-z"} };

        // Field name strings
        std::map<F, std::string> fieldMap { {F::U, "u"},
                                            {F::V, "v"},
                                            {F::W, "w"},
                                            {F::P, "p"} };

        // Iterate boundary patches
        std::string valueString;
        const pt::ptree *boundaryPatchTreePointer = nullptr;
        for (const auto& [patchEnum, patchString] : boundaryPatchMap) {
            boundaryPatchTreePointer = &( boundaryConditionsTree.get_child(patchString) );

            // Iterate fields
            for (const auto& [fieldEnum, fieldString] : fieldMap) {

                valueString = boundaryPatchTreePointer->get<std::string>(fieldString);
                inputData.boundaryConditions[fieldEnum][patchEnum] = ReadBoundaryValueString(valueString);

            }
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
        using F = Fields::ENUMDATA;
        const pt::ptree &schemesTree = solverTree.get_child( "Schemes" );
        std::string valueString;

        // Linearisation
        valueString = schemesTree.get<std::string>( "linearisation" );
        if        ( valueString == "Picard" ) {
            inputData.schemes.linearisation = Linearisation::Picard;
        } else if ( valueString == "Newton" ) {
            inputData.schemes.linearisation = Linearisation::Newton;
        } else {
            // THROW ERROR - invalid linearisation
        }


        // Advection scheme
        valueString = schemesTree.get<std::string>( "advectionScheme" );
        if        ( valueString == "upwind" ) {
            inputData.schemes.advectionScheme = AdvectionSchemes::Upwind;
        } else {
            // THROW ERROR - invalid advection scheme
        }

        // Face interpolation scheme
        valueString = schemesTree.get<std::string>( "faceInterpolationScheme" );
        if        ( valueString == "weightedLinear" ) {
            inputData.schemes.faceInterpolationScheme = FaceInterpolationSchemes::WeightedLinear;
        } else if ( valueString == "average" ) {
            inputData.schemes.faceInterpolationScheme = FaceInterpolationSchemes::Average;   
        } else {
            // THROW ERROR - invalid face interpolation scheme
        }


        // Momentum implicit relaxation
        std::vector<floatType> momentumRelaxation = schemesTree.get< std::vector<floatType> >( "implicitMomentumRelaxation" );
        inputData.schemes.implicitRelaxation[F::U] = momentumRelaxation[0];
        inputData.schemes.implicitRelaxation[F::V] = momentumRelaxation[1];
        inputData.schemes.implicitRelaxation[F::W] = momentumRelaxation[2];

        // Pressure implicit relaxation
        inputData.schemes.implicitRelaxation[F::P] = schemesTree.get<floatType>( "implicitPressureRelaxation" );

        // Max outer iterations
        inputData.schemes.maxOuterIterations = schemesTree.get<intType>( "maxOuterIterations" );

        // Max outer residuals
        inputData.schemes.maxOuterResiduals = schemesTree.get<floatType>( "maxOuterResiduals" );

    }



    void ReadLineSolverSettings( InputData &inputData, 
                                 const pt::ptree &planeSolverTree) 
    {
        const pt::ptree &lineSolverTree = planeSolverTree.get_child( "LineSolver" );
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = lineSolverTree.get<std::string>( "type" );
        if        ( valueString == "SUGS" ) {
            inputData.lineSolverSettings.type = LineSolvers::SUGS;
        } else {
            // THROW ERROR - invalid line solver type
        }

        // Max iterations
        inputData.lineSolverSettings.maxIterations = lineSolverTree.get<intType>( "maxIterations" );

        // Max residuals
        inputData.lineSolverSettings.maxResiduals = lineSolverTree.get<floatType>( "maxResiduals" );
        
        // Momentum relaxation
        std::vector<floatType> momentumRelaxation = lineSolverTree.get< std::vector<floatType> >( "momentumRelaxation" );
        inputData.lineSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.lineSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.lineSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        inputData.lineSolverSettings.relaxation[F::P] = lineSolverTree.get<floatType>( "pressureRelaxation" );

        // Line sweep direction
        valueString = lineSolverTree.get<std::string>( "lineSweepDirection" );
        inputData.lineSolverSettings.sweepDirection = ReadAxisDirection( valueString );

    }



    void ReadPlaneSolverSettings( InputData &inputData, 
                                 const pt::ptree &linearSolverTree) 
    {
        const pt::ptree &planeSolverTree = linearSolverTree.get_child( "PlaneSolver" );
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = planeSolverTree.get<std::string>( "type" );
        if        ( valueString == "SUGS" ) {
            inputData.planeSolverSettings.type = PlaneSolvers::SUGS;
        } else {
            // THROW ERROR - invalid plane solver type
        }

        // Max iterations
        inputData.planeSolverSettings.maxIterations = planeSolverTree.get<intType>( "maxIterations" );

        // Max residuals
        inputData.planeSolverSettings.maxResiduals = planeSolverTree.get<floatType>( "maxResiduals" );
        
        // Momentum relaxation
        std::vector<floatType> momentumRelaxation = planeSolverTree.get< std::vector<floatType> >( "momentumRelaxation" );
        inputData.planeSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.planeSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.planeSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        inputData.planeSolverSettings.relaxation[F::P] = planeSolverTree.get<floatType>( "pressureRelaxation" );

        // Plane sweep direction
        valueString = planeSolverTree.get<std::string>( "planeSweepDirection" );
        inputData.planeSolverSettings.sweepDirection = ReadAxisDirection( valueString );

        // Line solver
        ReadLineSolverSettings(inputData, planeSolverTree);

    }



    void ReadLinearSolverSettings( InputData &inputData, 
                                   const pt::ptree & solverTree) 
    {
        const pt::ptree &linearSolverTree = solverTree.get_child( "LinearSolver" );
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = linearSolverTree.get<std::string>( "type" );
        if        ( valueString == "SUGS" ) {
            inputData.linearSolverSettings.type = LinearSolvers::SUGS;
        } else {
            // THROW ERROR - invalid plane solver type
        }

        // Max iterations
        inputData.linearSolverSettings.maxIterations = linearSolverTree.get<intType>( "maxIterations" );

        // Max residuals
        inputData.linearSolverSettings.maxResiduals = linearSolverTree.get<floatType>( "maxResiduals" );

        // Momentum relaxation
        std::vector<floatType> momentumRelaxation = linearSolverTree.get< std::vector<floatType> >( "momentumRelaxation" );
        inputData.linearSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.linearSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.linearSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        inputData.linearSolverSettings.relaxation[F::P] = linearSolverTree.get<floatType>( "pressureRelaxation" );

        // Plane solver settings
        ReadPlaneSolverSettings(inputData, linearSolverTree);

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
        const pt::ptree &initialConditionsTree = tree.get_child("InitialConditions");
        std::string valueString;

        EnumVector<Fields, std::string> fieldKeys( {"u", "v", "w", "p"} );

        Fields::ENUMDATA field;
        for ( int f = 0; f != Fields::count; f++ ) {
            field = static_cast<Fields::ENUMDATA>( f );

            valueString = initialConditionsTree.get<std::string>( fieldKeys[field] );        
            inputData.initialConditions[field] = String2Type<CFD::floatType>( valueString );
        }

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

    return inputData;
}

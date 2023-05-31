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
    boundaryConditions(),
    axisTransformation()
    {};

CFD::InputData::AxisTransformationMap::AxisTransformationMap() :
    m_codeMap({ {CFD::BoundaryPatches::ENUMDATA::xPositive, CFD::BoundaryPatches::ENUMDATA::xPositive},
                {CFD::BoundaryPatches::ENUMDATA::xNegative, CFD::BoundaryPatches::ENUMDATA::xNegative},
                {CFD::BoundaryPatches::ENUMDATA::yPositive, CFD::BoundaryPatches::ENUMDATA::yPositive},
                {CFD::BoundaryPatches::ENUMDATA::yNegative, CFD::BoundaryPatches::ENUMDATA::yNegative},
                {CFD::BoundaryPatches::ENUMDATA::zPositive, CFD::BoundaryPatches::ENUMDATA::zPositive},
                {CFD::BoundaryPatches::ENUMDATA::zNegative, CFD::BoundaryPatches::ENUMDATA::zNegative} }),
    m_userMap( m_codeMap )
    {};


namespace
{

    using namespace CFD;

    /*-------------------------------------------------------------------------------------*\
                                         Helper Functions
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
    std::vector<T> ParseVectorString( const std::string &vecString, 
                                      const int &dim)
    {
        std::vector<T> vec;
        std::string::const_iterator stringIterator = vecString.begin();
        std::string valueString;
        int dimCount = 0;

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
                dimCount++;
            } else {
                valueString += *stringIterator;
            }
            stringIterator++;
        }

        if (*stringIterator != VECTOR_END_CHAR) {
            // throw ERROR - vector not closed
        }

        if (dimCount != dim) {
            // throw ERROR - expecting vector of size dim
        }

        return vec;

    }


    /*-------------------------------------------------------------------------------------*\
                                               Model
    \*-------------------------------------------------------------------------------------*/

    void ReadModel(InputData &inputData, 
                   const pt::ptree &tree)
    {
        using enum Axis::ENUMDATA;

        const pt::ptree &meshTree = tree.get_child("Model");

        // Domain
        const std::string &rhoString = meshTree.get<std::string>("rho");
        inputData.rho = String2Type<floatType>(rhoString);

        const std::string &nuString = meshTree.get<std::string>("nu");
        inputData.nu = String2Type<floatType>(nuString);

    }



    /*-------------------------------------------------------------------------------------*\
                                                Mesh
    \*-------------------------------------------------------------------------------------*/

    void ValidateSegmentBounds(const std::vector<InputData::MeshSegment> &meshSegment, 
                               const floatType &domainSize)
    {
        using v_size_type = std::vector<InputData::MeshSegment>::size_type;

        // The first one must start at zero
        if (meshSegment.front().lowerBound != 0.0) {
            // throw ERROR - must have mesh segment starting at coordinate 0
        }

        // Must end at the given domain size
        if (meshSegment.back().upperBound != domainSize ) {
            // throw ERROR - segments dont match given domain bounds
        }


        // Must be no gaps or overlaps in the segments
        for (v_size_type i = 1; i != meshSegment.size(); i++) {
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
        std::vector<CFD::floatType> tempBoundsVector;
        InputData::MeshSegment tempMeshSegment;
        for (auto segment : gridTree) {
            if (segment.first != "Segment") {
                // throw ERROR invalid child name
            }
            nCellsString     = segment.second.get<std::string>("nCells");
            boundsString     = segment.second.get<std::string>("bounds");
            biasFactorString = segment.second.get<std::string>("biasFactor");
            tempBoundsVector = ParseVectorString<floatType>(boundsString, 2);

            tempMeshSegment.nCells = String2Type<CFD::intType>(nCellsString);
            tempMeshSegment.biasFactor = String2Type<CFD::floatType>(biasFactorString);
            tempMeshSegment.lowerBound = tempBoundsVector[0];
            tempMeshSegment.upperBound = tempBoundsVector[1];

            meshSegments.push_back(tempMeshSegment);
        }
    }


    void ReadMesh(InputData &inputData, 
                  const pt::ptree &tree)
    {
        using enum Axis::ENUMDATA;

        const pt::ptree &meshTree = tree.get_child("Mesh");

        // Domain
        const std::string &domainSizeString = meshTree.get<std::string>("domain");
        std::vector<floatType> domainSizeTemp = ParseVectorString<floatType>(domainSizeString, 3);
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

        const pt::ptree &boundaryConditionsTree = tree.get_child("BoundaryConditions");

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


    Eigen::Matrix<CFD::intType, 3, 1> Axis2Vector(const CFD::BoundaryPatches::ENUMDATA axis)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        switch (axis)
        {
            case (BP::xPositive):
                return {1, 0, 0};
                
            case (BP::xNegative):
                return {-1, 0, 0};

            case (BP::yPositive):
                return {0, 1, 0};

            case (BP::yNegative):
                return {0, -1, 0};

            case (BP::zPositive):
                return {0, 0, 1};

            case (BP::zNegative):
                return {0, 0, -1};

            default:
                // throw ERROR - invalud axis enum
                return {0, 0, 0};
        }
    }


    CFD::BoundaryPatches::ENUMDATA Vector2Axis(const Eigen::Matrix<CFD::intType, 3, 1> &vector)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        if        ( ( vector.array() ==  Axis2Vector(BP::xPositive).array() ).all() ) {
            return BP::xPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::xNegative).array() ).all() ) {
            return BP::xNegative;

        } else if ( ( vector.array() ==  Axis2Vector(BP::yPositive).array() ).all() ) {
            return BP::yPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::yNegative).array() ).all() ) {
            return BP::yNegative;

        } else if ( ( vector.array() ==  Axis2Vector(BP::zPositive).array() ).all() ) {
            return BP::zPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::zNegative).array() ).all() ) {
            return BP::zNegative;

        } else {
            // throw ERROR - invalid unit vector
            return BP::xPositive;
        }
    }


    void ReadSweepDirections( InputData &inputData, 
                              const pt::ptree &solverTree) 
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;

        std::string valueString;
        BP planeSweepDirection, lineSweepDirection, pointSweepDirection;
        Eigen::Matrix<CFD::intType, 3, 1> planeSweepVector, lineSweepVector, pointSweepVector;

        // User input for plane sweep direction
        valueString = solverTree.get<std::string>("planeSweepDirection");
        planeSweepDirection = ReadAxisDirection( valueString );
        planeSweepVector = Axis2Vector( planeSweepDirection );

        // User input for line sweep direction
        valueString = solverTree.get<std::string>("lineSweepDirection");
        lineSweepDirection = ReadAxisDirection( valueString );
        lineSweepVector = Axis2Vector( lineSweepDirection );

        // Plane and line sweep direction must be orthogonal
        if ( abs( planeSweepVector.dot( lineSweepVector ) ) == 1 ) {
            // throw ERROR - plane and line sweep directions cannot be the same
        }

        // Point sweep direction, chosen to make right handed coordinate system
        pointSweepVector = lineSweepVector.cross( planeSweepVector );
        pointSweepDirection = Vector2Axis( pointSweepVector );

        // Update the axisTransformation map
        inputData.axisTransformation.Set( BP::xPositive, pointSweepDirection);
        inputData.axisTransformation.Set( BP::xNegative, Vector2Axis( -pointSweepVector ));

        inputData.axisTransformation.Set( BP::yPositive, lineSweepDirection);
        inputData.axisTransformation.Set( BP::yNegative, Vector2Axis( -lineSweepVector ));

        inputData.axisTransformation.Set( BP::zPositive, planeSweepDirection);
        inputData.axisTransformation.Set( BP::zNegative, Vector2Axis( -planeSweepVector ));
    }



    void ReadSchemes( InputData &inputData, 
                      const pt::ptree &solverTree) 
    {
        using F = Fields::ENUMDATA;
        const pt::ptree &schemesTree = solverTree.get_child("Schemes");
        std::string valueString;

        // Linearisation
        valueString = schemesTree.get<std::string>("linearisation");
        if        ( valueString == "Picard" ) {
            inputData.schemes.linearisation = Linearisation::Picard;
        } else if ( valueString == "Newton" ) {
            inputData.schemes.linearisation = Linearisation::Newton;
        } else {
            // THROW ERROR - invalid linearisation
        }

        // Momentum implicit relaxation
        valueString = schemesTree.get<std::string>("implicitMomentumRelaxation");
        std::vector<floatType> momentumRelaxation = ParseVectorString<floatType>(valueString, 3);
        inputData.schemes.implicitRelaxation[F::U] = momentumRelaxation[0];
        inputData.schemes.implicitRelaxation[F::V] = momentumRelaxation[1];
        inputData.schemes.implicitRelaxation[F::W] = momentumRelaxation[2];

        // Pressure implicit relaxation
        valueString = schemesTree.get<std::string>("implicitPressureRelaxation");
        inputData.schemes.implicitRelaxation[F::P] = String2Type<floatType>(valueString);

        // Max outer iterations
        valueString = schemesTree.get<std::string>("maxOuterIterations");
        inputData.schemes.maxOuterIterations = String2Type<intType>(valueString);

        // Max outer residuals
        valueString = schemesTree.get<std::string>("maxOuterResiduals");
        inputData.schemes.maxOuterResiduals = String2Type<floatType>(valueString);

    }



    void ReadPlaneSolverSettings( InputData &inputData, 
                                 const pt::ptree & solverTree) 
    {
        const pt::ptree &planeSolverTree = solverTree.get_child("PlaneSolver");
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = planeSolverTree.get<std::string>("type");
        if        ( valueString == "SUGS" ) {
            inputData.planeSolverSettings.type = PlaneSolvers::SUGS;
        } else {
            // THROW ERROR - invalid plane solver type
        }

        // Max iterations
        valueString = planeSolverTree.get<std::string>("maxIterations");
        inputData.planeSolverSettings.maxIterations = String2Type<intType>(valueString);

        // Max residuals
        valueString = planeSolverTree.get<std::string>("maxResiduals");
        inputData.planeSolverSettings.maxResiduals = String2Type<floatType>(valueString);
        
        // Momentum relaxation
        valueString = planeSolverTree.get<std::string>("momentumRelaxation");
        std::vector<floatType> momentumRelaxation = ParseVectorString<floatType>(valueString, 3);
        inputData.planeSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.planeSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.planeSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        valueString = planeSolverTree.get<std::string>("pressureRelaxation");
        inputData.planeSolverSettings.relaxation[F::P] = String2Type<floatType>(valueString);

    }



    void ReadLineSolverSettings( InputData &inputData, 
                                 const pt::ptree & solverTree) 
    {
        const pt::ptree &lineSolverTree = solverTree.get_child("LineSolver");
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = lineSolverTree.get<std::string>("type");
        if        ( valueString == "SUGS" ) {
            inputData.lineSolverSettings.type = LineSolvers::SUGS;
        } else {
            // THROW ERROR - invalid line solver type
        }

        // Max iterations
        valueString = lineSolverTree.get<std::string>("maxIterations");
        inputData.lineSolverSettings.maxIterations = String2Type<intType>(valueString);

        // Max residuals
        valueString = lineSolverTree.get<std::string>("maxResiduals");
        inputData.lineSolverSettings.maxResiduals = String2Type<floatType>(valueString);
        
        // Momentum relaxation
        valueString = lineSolverTree.get<std::string>("momentumRelaxation");
        std::vector<floatType> momentumRelaxation = ParseVectorString<floatType>(valueString, 3);
        inputData.lineSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.lineSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.lineSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        valueString = lineSolverTree.get<std::string>("pressureRelaxation");
        inputData.lineSolverSettings.relaxation[F::P] = String2Type<floatType>(valueString);

    }



    void ReadLinearSolverSettings( InputData &inputData, 
                                   const pt::ptree & solverTree) 
    {
        const pt::ptree &linearSolverTree = solverTree.get_child("LinearSolver");
        using F = Fields::ENUMDATA;
        std::string valueString;

        // Solver type
        valueString = linearSolverTree.get<std::string>("type");
        if        ( valueString == "SUGS" ) {
            inputData.linearSolverSettings.type = LinearSolvers::SUGS;
        } else {
            // THROW ERROR - invalid plane solver type
        }

        // Max iterations
        valueString = linearSolverTree.get<std::string>("maxIterations");
        inputData.linearSolverSettings.maxIterations = String2Type<intType>(valueString);

        // Max residuals
        valueString = linearSolverTree.get<std::string>("maxResiduals");
        inputData.linearSolverSettings.maxResiduals = String2Type<floatType>(valueString);

        // Momentum relaxation
        valueString = linearSolverTree.get<std::string>("momentumRelaxation");
        std::vector<floatType> momentumRelaxation = ParseVectorString<floatType>(valueString, 3);
        inputData.linearSolverSettings.relaxation[F::U] = momentumRelaxation[0];
        inputData.linearSolverSettings.relaxation[F::V] = momentumRelaxation[1];
        inputData.linearSolverSettings.relaxation[F::W] = momentumRelaxation[2];

        // Pressure relaxation
        valueString = linearSolverTree.get<std::string>("pressureRelaxation");
        inputData.linearSolverSettings.relaxation[F::P] = String2Type<floatType>(valueString);


        // Plane solver settings
        ReadPlaneSolverSettings(inputData, linearSolverTree);

        // Line solver settings
        ReadLineSolverSettings(inputData, linearSolverTree);

    }



    void ReadSolver(InputData &inputData, 
                    const pt::ptree &tree)
    {
        const pt::ptree &solverTree = tree.get_child("Solver");

        // Sweep direction transformation data
        ReadSweepDirections(inputData, solverTree);

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


    /*-------------------------------------------------------------------------------------*\
                                    Axis Transformations
    \*-------------------------------------------------------------------------------------*/

    // The problem is remapped so that the plane sweeping direction is always in the z direction and
    // the line sweeping direction is always in the y direction (in the code). This is more memory
    // efficient and is simpler to implement.

    // Remaps the users boundary conditions
    void TransformBoundaryConditions(InputData &inputData)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        using F = CFD::Fields::ENUMDATA;

        // Temporary for boundary conditions as user specifies them
        InputData::BoundaryConditionData boundaryConditionsUser = inputData.boundaryConditions;
        

        EnumFor<Fields>([&] (F field) {
            
            EnumFor<BoundaryPatches>([&] (BP patch) { 

                inputData.boundaryConditions[field][patch] = boundaryConditionsUser[field][ inputData.axisTransformation.UserPatch( patch ) ];

            } );

        } );

    }


    // Reverse a mesh in a given axis
    void ReverseMesh(std::vector< InputData::MeshSegment > &meshSegments)
    {
        // Reverse the order of the segments
        std::reverse(meshSegments.begin(), meshSegments.end());

        // Now flip each segment
        for (auto &segment : meshSegments) {
            std::swap(segment.upperBound, segment.lowerBound);
            segment.upperBound = - segment.upperBound;
            segment.lowerBound = - segment.lowerBound;
            segment.biasFactor = - segment.biasFactor;
        }
    }


    // Remaps the user mesh
    void TransformMesh(InputData &inputData)
    {
        using enum Axis::ENUMDATA;
        using enum BoundaryPatches::ENUMDATA;


        // Create temporary copy of mesh data to take data from
        EnumVector<Axis, std::vector<InputData::MeshSegment> > userMeshSegments = inputData.meshSegments;
        floatVector3 userDomainSize;

        BoundaryPatches::ENUMDATA userPatch;
        Axis::ENUMDATA userAxis;

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            userPatch = inputData.axisTransformation.UserPatch( PositivePatch[ axis ] );
            userAxis = BoundaryPatchAxis[ userPatch ];

            inputData.meshSegments[ axis ] = userMeshSegments[ userAxis ];
            inputData.domainSize( axis )   = userDomainSize[ userAxis ];

            if ( userPatch == NegativePatch[ userAxis ] ) {
                ReverseMesh(inputData.meshSegments[ axis ]);
            }

        } );

    }


    // Transform an EnumVector of fields data. Only transforms the momentum equations part.
    void TransformFieldVector( EnumVector<Fields, floatType> &fieldsVector,
                               const InputData::AxisTransformationMap& axisTransformation )
    {
        using F = Fields::ENUMDATA;
        EnumVector<Axis, F> axisField({ F::U, F::V, F::W });

        // Create temporary copy to move data from 
        EnumVector<Fields, floatType> userFieldsVector = fieldsVector;

        BoundaryPatches::ENUMDATA userPatch;
        Axis::ENUMDATA userAxis;

        EnumFor<Axis>([&] (Axis::ENUMDATA axis) { 

            userPatch = axisTransformation.UserPatch( PositivePatch[ axis ] );
            userAxis = BoundaryPatchAxis[ userPatch ];

            fieldsVector[ axisField[axis] ] = userFieldsVector[ axisField[userAxis] ];

        } );
    }


    // Remaps the initial conditions
    void TransformInitialConditions( InputData &inputData )
    {
        TransformFieldVector( inputData.initialConditions, inputData.axisTransformation );
    }


    // Remaps any solver settings that have direction dependence
    void TransformSolver( InputData &inputData )
    {
        TransformFieldVector( inputData.schemes.implicitRelaxation      , inputData.axisTransformation );
        TransformFieldVector( inputData.linearSolverSettings.relaxation , inputData.axisTransformation );
        TransformFieldVector( inputData.planeSolverSettings.relaxation  , inputData.axisTransformation );
        TransformFieldVector( inputData.lineSolverSettings.relaxation   , inputData.axisTransformation );
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
    
    TransformBoundaryConditions(inputData);
    TransformInitialConditions(inputData);
    TransformMesh(inputData);
    TransformSolver(inputData);

    return inputData;
}


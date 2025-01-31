#include "InputProcessing.h"
#include "InputParser.h"
#include "CSVReader.h"
#include "IOTools.h"

#include "../FiniteVolume/FiniteVolume.h"
#include "../Types.h"
#include "../Tools/FVLookups.h"

#include "Boost/boost/property_tree/ptree.hpp"
#include <Eigen/Geometry>

#include <utility>
#include <optional>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <filesystem>

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
            std::cout << "Reading input file '" + inputFilename + "' ... ";
            inputData = CFD::ReadInputData(inputFilename);
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
            vec.push_back( CFD::IOTOOLS::String2Type<T>(valueString) );
            break;
        }

        if ( *stringIterator == VECTOR_DELIMITER_CHAR ) {
            vec.push_back( CFD::IOTOOLS::String2Type<T>(valueString) );
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

        std::string valueString = modelTree.get<std::string>("transient");
        if ( valueString == "yes" ) {
            inputData.transient = true;
        } else if ( valueString == "no" ) {
            inputData.transient = false;
        } else {
            throw std::runtime_error(  "'" + valueString + "' is not a transient option. Must be either 'yes' or 'no'." );
        }  
    }



    /*-------------------------------------------------------------------------------------*\
                                                Mesh
    \*-------------------------------------------------------------------------------------*/

    void ValidateSegmentBounds(const std::vector<InputData::MeshSegment> &meshSegment, 
                               const floatType &domainSize)
    {

        // The first one must start at zero
        if (meshSegment.front().startCoordinate != 0.0) {
            throw std::runtime_error( "Mesh segments must start at coordinate 0." );
        }

        // Must end at the given domain size
        if (meshSegment.back().endCoordinate != domainSize ) {
            throw std::runtime_error( "Mesh segments to not match given domain bounds." );
        }


        // Must be no gaps or overlaps in the segments
        for (size_t i = 1; i != meshSegment.size(); i++) {
            if ( meshSegment[i].startCoordinate != meshSegment[i-1].endCoordinate ) {
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
            tempMeshSegment.startCoordinate = tempBoundsVector[0];
            tempMeshSegment.endCoordinate   = tempBoundsVector[1];

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
        auto sortComparison = [](const auto& i, const auto& j) { return i.startCoordinate < j.startCoordinate; };
        std::sort( inputData.meshSegments[X].begin(), inputData.meshSegments[X].end(), sortComparison);
        std::sort( inputData.meshSegments[Y].begin(), inputData.meshSegments[Y].end(), sortComparison );
        std::sort( inputData.meshSegments[Z].begin(), inputData.meshSegments[Z].end(), sortComparison );

        // Check that the bounds are valid
        ValidateSegmentBounds(inputData.meshSegments[X], inputData.domainSize[X]);
        ValidateSegmentBounds(inputData.meshSegments[Y], inputData.domainSize[Y]);
        ValidateSegmentBounds(inputData.meshSegments[Z], inputData.domainSize[Z]);

    }


    /*-------------------------------------------------------------------------------------*\
                                         Solid Geometry
    \*-------------------------------------------------------------------------------------*/

    void ReadBlockData( InputData &inputData,
                        const std::pair<const std::string, pt::ptree> &solidObject )
    {
        InputData::SolidBlockData tempBlockData;
        std::vector< floatType > tempCenterPosition, tempDimensions, tempRotation;
        tempCenterPosition = solidObject.second.get< std::vector<floatType> >( "centerPosition" ); 
        tempDimensions     = solidObject.second.get< std::vector<floatType> >( "dimensions" ); 
        tempRotation       = solidObject.second.get< std::vector<floatType> >( "rotation" ); 

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            tempBlockData.centerPosition(axis) = tempCenterPosition[axis];
            tempBlockData.dimensions(axis)     = tempDimensions[axis];
            tempBlockData.rotation(axis)       = tempRotation[axis];
        } );

        inputData.solidBlocks.push_back( tempBlockData );
    }



    void ReadSphereData( InputData &inputData,
                        const std::pair<const std::string, pt::ptree> &solidObject )
    {
        InputData::SolidSphereData tempSphereData;
        std::vector< floatType > tempCenterPosition;
        floatType tempDiameter;
        tempCenterPosition = solidObject.second.get< std::vector<floatType> >( "centerPosition" ); 
        tempDiameter       = solidObject.second.get< floatType >( "diameter" );

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            tempSphereData.centerPosition(axis) = tempCenterPosition[axis];
        } );
        tempSphereData.diameter = tempDiameter;

        inputData.solidSpheres.push_back( tempSphereData );
    }



    void ReadGeometryFromFileData( InputData &inputData,
                                   const std::pair<const std::string, pt::ptree> &solidObject )
    {
        std::string valueString = solidObject.second.get< std::string >( "type" );
        if ( valueString == "STL" ) {
            valueString = solidObject.second.get< std::string >( "filename" );
            inputData.geometrySTLFiles.push_back( valueString );
        } else {
            throw std::runtime_error(  "'" + valueString + "' is not a supported geometry file type." );
        }
    }



    void ReadSolidGeometry( InputData &inputData, 
                            const pt::ptree &tree)
    {
        // It's ok if no geometry is specified
        boost::optional<const pt::ptree &> solidGeometryTreeOptional = tree.get_child_optional( "SolidGeometry" );
        if ( !solidGeometryTreeOptional ) {
            inputData.hasIBGeometry = false;
            return;
        }
            
        const pt::ptree &solidGeometryTree = solidGeometryTreeOptional.get();

        // Geometry boundary treatement
        std::string valueString = solidGeometryTree.get<std::string>( "boundaryTreatement" );
        if        ( valueString == "directionalImmersedBoundary" ) {
            inputData.geoemtryBoundaryTreatement = GeometryBoundaryTreatement::DirectionalImmersedBoundary;
        } else if ( valueString == "staircase" ) {
            inputData.geoemtryBoundaryTreatement = GeometryBoundaryTreatement::Staircase;
        } else {
            throw std::runtime_error(  "'" + valueString + "' is not a valid solid geometry treatement method." );
        }

        // Read blocks
        for (auto solidObject : solidGeometryTree) {

            if ( solidObject.first == "Block" ) {
                ReadBlockData( inputData, solidObject );
                inputData.hasIBGeometry = true;
                continue;
            }

            if ( solidObject.first == "Sphere" ) {
                ReadSphereData( inputData, solidObject );
                inputData.hasIBGeometry = true;
                continue;
            }

            if ( solidObject.first ==  "FromFile") {
                ReadGeometryFromFileData( inputData, solidObject );
                inputData.hasIBGeometry = true;
                continue;
            }

        }        
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

        return IOTOOLS::String2Type<floatType>(bcValueString);
    }


    // Read uniform BC value. Assumes that the string iterator is after the comma
    InputData::Profile1D GetProfile1DValue( std::string::const_iterator &stringIterator,
                                            const std::string &boundaryString,
                                            const BoundaryPatches::ENUMDATA boundaryPatch,
                                            const std::string &inputFileDirectory )
    {
        InputData::Profile1D profile1D;

        // Read the filename, which is after the comma and up to the end of the line
        std::string filename;
        for ( /* NULL */; stringIterator != boundaryString.end(); stringIterator++ ) { 
            filename += *stringIterator;
        }

        // Add the directory to the input filename
        if ( !inputFileDirectory.empty() ){
            filename = inputFileDirectory + "/" + filename;
        }
            

        // Read in profile data from csv file
        std::vector< std::vector< std::string > > profileData = ReadCSV( filename );

        // First column tells us the axis
        std::string axisString = IOTOOLS::RemoveWhitespace( profileData[0][0] );
        if        ( axisString == "x" ) {
            profile1D.axis = Axis::X;
        } else if ( axisString == "y" ) {
            profile1D.axis = Axis::Y;
        } else if ( axisString == "z" ) {
            profile1D.axis = Axis::Z;
        } else {
            throw std::runtime_error( "Invalid axis name in profile file '" + filename + "'."  );
        }

        // The profile must be orthogonal to the normal of the given patch
        if ( profile1D.axis == LUT::BoundaryPatchAxis[ boundaryPatch ] ) {
            throw std::runtime_error( "Profile axis must be in the plane of the specified patch." );
        }

        // First column is the coordinate points, second column is the actual data
        intType nRows = static_cast<intType>( profileData.size() );
        intType nHeaderRows = 1;
        profile1D.coordinates = Tensor1D( nRows - nHeaderRows );
        profile1D.values      = Tensor1D( nRows - nHeaderRows );

        for ( intType i = 0; i != nRows-nHeaderRows; i++ ) {
            profile1D.coordinates(i) = IOTOOLS::String2Type<floatType>( profileData[ static_cast<size_t>( i+nHeaderRows) ][0]  );
            profile1D.values(i)      = IOTOOLS::String2Type<floatType>( profileData[ static_cast<size_t>( i+nHeaderRows) ][1]  );
        }

        return profile1D;
    }



    InputData::BoundaryConditionInputData ReadBoundaryValueString( const std::string &fieldString,
                                                                   const pt::ptree &boundaryPatchTree,
                                                                   const BoundaryPatches::ENUMDATA boundaryPatch,
                                                                   const std::string &inputFileDirectory )
    {
        using BC = BoundaryConditions::ENUMDATA;

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
        if        (bcTypeString == "uniform") {

            CheckForBCValue( stringIterator, bcTypeString );
            bcStruct.uniformValue = GetUniformBCValue( stringIterator, boundaryString );
            bcStruct.hasUniformValue = true;
            bcStruct.type = BC::fixed;

        } else if (bcTypeString == "profile1D") {

            CheckForBCValue( stringIterator, bcTypeString );
            bcStruct.profile1D = GetProfile1DValue( stringIterator, boundaryString, boundaryPatch, inputFileDirectory );
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

        } else if (bcTypeString == "periodic") {

            bcStruct.type = BC::periodic;
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
            inputData.boundaryConditions.U[X][patchEnum] = ReadBoundaryValueString( "u", *boundaryPatchTreePointer, patchEnum, inputData.inputFileDirectory );
            inputData.boundaryConditions.U[Y][patchEnum] = ReadBoundaryValueString( "v", *boundaryPatchTreePointer, patchEnum, inputData.inputFileDirectory );
            inputData.boundaryConditions.U[Z][patchEnum] = ReadBoundaryValueString( "w", *boundaryPatchTreePointer, patchEnum, inputData.inputFileDirectory );

            // Pressure
            inputData.boundaryConditions.P[patchEnum] = ReadBoundaryValueString( "p", *boundaryPatchTreePointer, patchEnum, inputData.inputFileDirectory );

        }


        // Verify that periodic boundary conditions are set for both sides
        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

            ForAllFieldData( [&] (intType f) {
                if ( inputData.boundaryConditions[f][bp].type == BoundaryConditions::periodic ) {

                    Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ bp ];
                    BoundaryPatches::ENUMDATA bpOpposite = ( bp == LUT::PositivePatch[axis] ) ? LUT::NegativePatch[axis] : LUT::PositivePatch[axis];
                    if ( inputData.boundaryConditions[f][bpOpposite].type != BoundaryConditions::periodic ) {
                         throw std::runtime_error( "Specification of periodic boundary conditions requires both sides of domain to be specified as periodic." );
                    }

                }
            } );

        } );

    }


    /*-------------------------------------------------------------------------------------*\
                                           Solver
    \*-------------------------------------------------------------------------------------*/


    CFD::BoundaryPatches::ENUMDATA String2BoundaryPatch(const std::string &bpString)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        if        ( bpString == "+x" ) {
            return BP::xPositive;

        } else if ( bpString == "-x" ) {
            return BP::xNegative;

        } else if ( bpString == "+y" ) {
            return BP::yPositive;
        
        } else if ( bpString == "-y" ) {
            return BP::yNegative;

        } else if ( bpString == "+z" ) {
            return BP::zPositive;

        } else if ( bpString == "-z" ) {
            return BP::zNegative;

        } else {
            throw std::runtime_error( "'" + bpString + "' is an invalid axis direction" );
            return BP::xPositive;
        }
    }


    void ReadSchemes( InputData &inputData, 
                      const pt::ptree &solverTree) 
    {
        const pt::ptree &schemesTree = solverTree.get_child( "Schemes" );
        std::string valueString;

        // Time scheme
        if ( !inputData.transient ) {
            inputData.schemes.timeScheme = TimeSchemes::Steady;
        } else {
            valueString = schemesTree.get<std::string>( "timeScheme" );
            if        ( valueString == "backwardsEuler" ) {
                inputData.schemes.timeScheme = TimeSchemes::BackwardsEuler;
            } else if ( valueString == "backwardsThreeLevel" ) {
                inputData.schemes.timeScheme = TimeSchemes::BackwardsThreeLevel;
            } else {
                throw std::runtime_error( "'" + valueString + "' is not a valid time discretisation scheme." );
            }
            inputData.schemes.timeStep = schemesTree.get<floatType>("timeStepSize");
            inputData.schemes.numberOfTimesteps = schemesTree.get<intType>("numberOfTimesteps");
        }

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
        } else if ( valueString == "central" ) {
            inputData.schemes.advectionScheme = AdvectionSchemes::Central;
        } else if ( valueString == "SOU" ) {
            inputData.schemes.advectionScheme = AdvectionSchemes::SOU;
        } else if ( valueString == "QUICK" ) {
            inputData.schemes.advectionScheme = AdvectionSchemes::QUICK;
        } else  {
            throw std::runtime_error( "'" + valueString + "' is not a valid advection scheme." );
        }


        // Advection blending factor. Default this to 1
        boost::optional<floatType> advectionBlendingFactorOptional = schemesTree.get_optional<floatType>( "advectionBlendingFactor" );
        if ( !advectionBlendingFactorOptional ) {
            inputData.schemes.advectionBlendingFactor = 1.0f;
        } else {
            inputData.schemes.advectionBlendingFactor = advectionBlendingFactorOptional.get();
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
        floatType maxContinuityOuterResidual = schemesTree.get<floatType>( "maxContinuityOuterResidual" );
        inputData.schemes.maxOuterResiduals.P = maxContinuityOuterResidual;

        std::vector<floatType> maxMomentumOuterResiduals = schemesTree.get< std::vector<floatType> >( "maxMomentumOuterResiduals" );
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            inputData.schemes.maxOuterResiduals.U[axis] = maxMomentumOuterResiduals[axis];
        } );
        
    }



    void ReadLinearSolverSettings( InputData &inputData, 
                                   const pt::ptree & solverTree) 
    {
        const pt::ptree &linearSolverTree = solverTree.get_child( "LinearSolver" );
        std::string valueString;

        // Solver type
        valueString = linearSolverTree.get<std::string>( "type" );
        if        ( valueString == "nestedLineSymmetric" ) {
            inputData.linearSolverSettings.type = LinearSolvers::nestedLineSymmetric;
        } else if ( valueString == "domainSymmetric" ) {
            inputData.linearSolverSettings.type = LinearSolvers::domainSymmetric;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a valid linear solver type." );
        }

        // Max iterations
        inputData.linearSolverSettings.maxIterations = linearSolverTree.get<intType>( "maxIterations" );

        // Max residuals
        inputData.linearSolverSettings.maxResiduals = linearSolverTree.get<floatType>( "maxResiduals" );

        // Momentum relaxation
        std::vector<floatType> momentumRelaxation = linearSolverTree.get< std::vector<floatType> >( "momentumRelaxation" );
        inputData.linearSolverSettings.relaxation.U[Axis::X] = momentumRelaxation[0];
        inputData.linearSolverSettings.relaxation.U[Axis::Y] = momentumRelaxation[1];
        inputData.linearSolverSettings.relaxation.U[Axis::Z] = momentumRelaxation[2];

        // Pressure relaxation
        inputData.linearSolverSettings.relaxation.P = linearSolverTree.get<floatType>( "pressureRelaxation" );

        // Plane sweep direction
        valueString = linearSolverTree.get<std::string>( "planeSweepDirection" );
        inputData.linearSolverSettings.planeSweepDirection = String2BoundaryPatch( valueString );

        // Line sweep direction
        valueString = linearSolverTree.get<std::string>( "lineSweepDirection" );
        inputData.linearSolverSettings.lineSweepDirection = String2BoundaryPatch( valueString );

    }



    void ReadMultigridSettings( InputData &inputData, 
                                const pt::ptree &solverTree) 
    {
        const pt::ptree &multigridTree = solverTree.get_child( "Multigrid" );
        std::string valueString;

        // Multigrid cycle type
        valueString = multigridTree.get<std::string>( "cycle" );
        if        ( valueString == "V" ) {
            inputData.multigridSettings.cycle = MultigridCycleType::V;
        } else if ( valueString == "F" ) {
            inputData.multigridSettings.cycle = MultigridCycleType::F;
        } else if ( valueString == "W" ) {
            inputData.multigridSettings.cycle = MultigridCycleType::W;
        } else {
            throw std::runtime_error( "'" + valueString + "' is not a linear solver type." );
        }

        // Max coarse levels
        inputData.multigridSettings.maxCoarseLevels = multigridTree.get<intType>( "maxCoarseLevels" );

        // Max iterations
        inputData.multigridSettings.preSmoothingIterations  = multigridTree.get<intType>( "preSmoothingIterations" );
        inputData.multigridSettings.postSmoothingIterations = multigridTree.get<intType>( "postSmoothingIterations" );
        inputData.multigridSettings.maxCoarseGridIterations    = multigridTree.get<intType>( "maxCoarseGridIterations" );
        inputData.multigridSettings.fineGridIterations      = multigridTree.get<intType>( "fineGridIterations" );

        // Max residuals
        inputData.multigridSettings.maxCoarseGridResiduals    = multigridTree.get<floatType>( "maxCoarseGridResiduals" );

    }



    void ReadSolver(InputData &inputData, 
                    const pt::ptree &tree)
    {
        const pt::ptree &solverTree = tree.get_child("Solver");

        // Read discretisation schemes
        ReadSchemes(inputData, solverTree);

        // Read linear solver (plane sweeping) settings
        ReadLinearSolverSettings(inputData, solverTree);

        // Read multigrid settings
        ReadMultigridSettings(inputData, solverTree);

    }


    /*-------------------------------------------------------------------------------------*\
                                     Initial Conditions
    \*-------------------------------------------------------------------------------------*/

    void ReadInitialConditions( InputData &inputData,
                                const pt::ptree &tree )
    {
        using enum Axis::ENUMDATA;
        const pt::ptree &initialConditionsTree = tree.get_child("InitialConditions");

        #if defined( CFD_HAS_VTK_LIB )
            std::string valueString = initialConditionsTree.get<std::string>( "type" );
            if        ( valueString == "uniform" ) {
                inputData.initialConditionType = InputData::InitialConditionTypes::uniform;
            } else if ( valueString == "vtkFile" ) {
                inputData.initialConditionType = InputData::InitialConditionTypes::vtkFile;
            } else {
                throw std::runtime_error( "'" + valueString + "' is not a valid initial condition specification type." );
            }

            switch ( inputData.initialConditionType ) {
                case InputData::InitialConditionTypes::uniform:
                    inputData.constantInitialConditions.U[X] = initialConditionsTree.get<floatType>( "u" );
                    inputData.constantInitialConditions.U[Y] = initialConditionsTree.get<floatType>( "v" );
                    inputData.constantInitialConditions.U[Z] = initialConditionsTree.get<floatType>( "w" );
                    inputData.constantInitialConditions.P    = initialConditionsTree.get<floatType>( "p" );
                    break;

                case InputData::InitialConditionTypes::vtkFile:
                    inputData.initialConditionsFieldFilename = initialConditionsTree.get<std::string>( "filename" );
                    break;
            }
        #else
            std::string valueString = initialConditionsTree.get<std::string>( "type" );
            if        ( valueString == "uniform" ) {
                inputData.initialConditionType = InputData::InitialConditionTypes::uniform;
            } else if ( valueString == "vtkFile" ) {
                throw std::runtime_error( "Must have and compile with VTK library to use '" + valueString + "' initial condition." );;
            } else {
                throw std::runtime_error( "'" + valueString + "' is not a valid initial condition specification type." );
            }
            inputData.constantInitialConditions.U[X] = initialConditionsTree.get<floatType>( "u" );
            inputData.constantInitialConditions.U[Y] = initialConditionsTree.get<floatType>( "v" );
            inputData.constantInitialConditions.U[Z] = initialConditionsTree.get<floatType>( "w" );
            inputData.constantInitialConditions.P    = initialConditionsTree.get<floatType>( "p" );
        #endif
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
        VerifyWriteDirectory( inputData.residualHistoryFilename );
        for ( const auto &probe : inputData.probes ) {
            VerifyWriteDirectory( probe.filename );
        }
    }


    void ReadOutput( InputData &inputData,
                    const pt::ptree &tree )
    {
        const pt::ptree &outputTree = tree.get_child( "Output" );

        // Residual history filename
        inputData.residualHistoryFilename = outputTree.get<std::string>( "residualHistoryFilename" );

        // Profiling information filename
        inputData.profilingFilename = outputTree.get<std::string>( "profilingFilename" );

        // Field output filename
        inputData.fieldOutputFilename = outputTree.get<std::string>( "fieldOutputFilename" );

        // Geometry out filename, the user does not have to output the geometry to file
        inputData.outputGeometry = false;
        if ( inputData.hasIBGeometry ) {

            boost::optional<std::string> geometryOutputFilenameOptional = outputTree.get_optional<std::string>( "geometryOutputFilename" );
            if ( geometryOutputFilenameOptional ) {
                inputData.outputGeometry = true;
                inputData.geometryOutputFilename = geometryOutputFilenameOptional.get();
                return;
            }

        }
        
        
        // Write interval for fields
        inputData.fieldWriteInterval = outputTree.get<intType>( "fieldWriteInterval" );

        // Probes
        ReadMonitors( inputData, outputTree );

        // Verify the output filenames
        VerifyOutputFiles( inputData ); 
    }

}   // end anonymous namepsace


// Parse input file and read into InputData structure
CFD::InputData CFD::ReadInputData(const std::string &inputFilename) 
{
    pt::ptree tree = INP::ParseFile(inputFilename);

    InputData inputData;
    inputData.inputFileDirectory = IOTOOLS::RelativePath( inputFilename );

    ReadModel(inputData, tree);
    ReadMesh(inputData, tree);
    ReadSolidGeometry(inputData, tree);
    ReadBoundaryConditions(inputData, tree);
    ReadInitialConditions(inputData, tree);
    ReadSolver(inputData, tree);
    ReadOutput(inputData, tree);

    return inputData;
}

// Parse input file and just read sweep directions
std::tuple< BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA > CFD::ReadSweepDirections( const std::string &inputFilename )
{
    pt::ptree tree = INP::ParseFile(inputFilename);
    const pt::ptree &linearSolverTree = tree.get_child("Solver").get_child("LinearSolver");

    // Plane sweep direction
    std::string valueString = linearSolverTree.get<std::string>( "planeSweepDirection" );
    BoundaryPatches::ENUMDATA planeSweepDirection = String2BoundaryPatch( valueString );

    // Line sweep direction
    valueString = linearSolverTree.get<std::string>( "lineSweepDirection" );
    BoundaryPatches::ENUMDATA lineSweepDirection = String2BoundaryPatch( valueString );

    return { planeSweepDirection, lineSweepDirection };
}

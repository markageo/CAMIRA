#include "FiniteVolume.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

#include <cmath>

namespace CFD
{

/*-------------------------------------------------------------------------------------*\
                                         Mesh
\*-------------------------------------------------------------------------------------*/

namespace
{

std::vector<floatType> CalculateGrowthRates(const std::vector<InputData::MeshSegment> &meshSegments)
{
    size_t nSegments = meshSegments.size();
    std::vector<floatType> growthRates(nSegments);

    // Negative growth rate means shrinking grid
    for (size_t i = 0; i != nSegments; i++) {
        growthRates[i] = std::pow( std::abs(meshSegments[i].biasFactor) , 1.0f / static_cast<floatType>( meshSegments[i].nCells - 1 ) );
        if (meshSegments[i].biasFactor < 0 )   
            growthRates[i] = 1.0f/growthRates[i];
    }

    return growthRates;
}


void CalculateCellLengths(Tensor1D &cellLengths, 
                          const std::vector<InputData::MeshSegment> &meshSegments, 
                          const std::vector<floatType> &growthRates)
{
    size_t nSegments = meshSegments.size();

    floatType segmentLength, firstCellLength, geometricFactor;
    int cellIndex = 0;
    for (size_t s = 0; s != nSegments; s++) {    // Segments

        if (growthRates[s] != 1.0f) { 
            geometricFactor = (1.0f - std::pow( growthRates[s], static_cast<floatType>(meshSegments[s].nCells) )) / (1.0f - growthRates[s]);   // geometric series formula
        } else {
            geometricFactor = static_cast<floatType>( meshSegments[s].nCells );
        }
        segmentLength = meshSegments[s].upperBound - meshSegments[s].lowerBound;
        firstCellLength = segmentLength / geometricFactor; 

        for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
            cellLengths( cellIndex ) = firstCellLength*std::pow( growthRates[s], static_cast<floatType>(i) );
            cellIndex++;
        }
    }

}


void CalculateCellCenters(Tensor1D &cellCenters, 
                          const Tensor1D &cellLengths,
                          const floatType startPosition)
{
    intType nCellsTotal = cellLengths.size();

    floatType previousCellPosition = startPosition, previousCellLength = 0.0f;
    for (intType i = 0; i != nCellsTotal; i++) {
        cellCenters(i) = previousCellPosition + previousCellLength/2.0f + cellLengths(i)/2.0f;
        previousCellPosition = cellCenters(i);
        previousCellLength = cellLengths(i);
    }

}


void CalculateCellCenterDiffInv(Tensor1D &cellCenterDiffInv, 
                                const Tensor1D &cellCenters)
{
    // First and last element dont correspond to valid values
    intType nFaces = cellCenters.size() + 1;

    for (intType i = 1; i != nFaces-1; i++) {
        cellCenterDiffInv(i) = 1.0f/( cellCenters(i) - cellCenters(i-1) );
    }
}


void CalculateCellFaces(Tensor1D &cellFaces, 
                        const Tensor1D &cellLengths,
                        const floatType startPosition)
{
    intType nFaces = cellFaces.size();

    cellFaces(0) = startPosition;
    for (intType i = 1; i != nFaces; i++) {
        cellFaces(i) = cellFaces(i-1) + cellLengths(i-1);
    }
}


void CalculateCellFaceAreas(Tensor2D &cellFaceAreas, 
                            const Tensor1D &cellLengths_x, 
                            const Tensor1D &cellLengths_y)
{
    for (int j = 0; j != cellLengths_y.dimension(0); j++) {
        for (int i = 0; i != cellLengths_x.dimension(0); i++) {
            cellFaceAreas(i, j) = cellLengths_x(i) * cellLengths_y(j);
        }
    }
}


void CalculateInterpolationFactors_WeightedLinear( Tensor1D &interpFactors, 
                                                   const Tensor1D &cellCenters, 
                                                   const Tensor1D &cellFaces) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) 
                         / ( cellCenters(i) - cellCenters(i-1) );
    }
}


void CalculateInterpolationFactors_Average( Tensor1D &interpFactors ) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = 0.5f;
    }
}


Mesh::ExtrapFactorsStruct GetExtrapolationFactors(const Tensor1D &cellLengths, 
                                                  const intType fieldIndex_p, 
                                                  const intType fieldIndex_a)
{
    floatType extrapFactor_p = ( 2.0f*cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );

    floatType extrapFactor_a = - ( cellLengths(fieldIndex_p) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );
    
    return Mesh::ExtrapFactorsStruct{ extrapFactor_p, extrapFactor_a };
}


void CalculateExtrapolationFactors(EnumVector<BoundaryPatches, Mesh::ExtrapFactorsStruct > &extrapFactors, 
                                   const EnumVector<Axis, Tensor1D> &cellLengths, 
                                   const Axis::ENUMDATA axis)
{  

    using enum BoundaryPatches::ENUMDATA;
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA patchPositive = LUT::PositivePatch[ axis ], 
                              patchNegative = LUT::NegativePatch[ axis ];

    // If mesh is only 1 cell think (such as in a 2D simulation), make the extrapolatino equal to the single cell
    if ( cellLengths[axis].size() == 1 ) {
        extrapFactors[patchPositive].a = 0.0f;
        extrapFactors[patchPositive].p = 1.0f;

        extrapFactors[patchNegative].a = 0.0f;
        extrapFactors[patchNegative].p = 1.0f;
        return;
    }


    
    intType fieldIndex_p, fieldIndex_a; // Boundary cell node and the adjacent one

    // Positive patch boundary
    fieldIndex_p = cellLengths[axis].size() - 1;
    fieldIndex_a = fieldIndex_p - 1;
    extrapFactors[patchPositive] = GetExtrapolationFactors(cellLengths[axis], fieldIndex_p, fieldIndex_a);

    // Negative patch boundary
    fieldIndex_p = 0;
    fieldIndex_a = fieldIndex_p + 1;
    extrapFactors[patchNegative] = GetExtrapolationFactors(cellLengths[axis], fieldIndex_p, fieldIndex_a);
}


intType TotalCells(const std::vector<InputData::MeshSegment> &meshSegments)
{
    intType totalCells = 0;
    for (auto segment : meshSegments) {
        totalCells += segment.nCells;
    }
    return totalCells;
}



iArray3 NumberOfFaces( const iArray3 &nCells,
                       Axis::ENUMDATA axis)
{
    iArray3 nFaces;
    EnumFor<Axis>( [&] (Axis::ENUMDATA a) {
        nFaces(a) = nCells(a);
    } );
    nFaces(axis) += 1;  // There is one more faces than cells in the normal direction
    return nFaces;
}

}   // end anonymous namespace




// Constructor, creates the mesh from user inputdata
Mesh::Mesh(const InputData &inputData) :
    nCells( { TotalCells(inputData.meshSegments[Axis::X]),  TotalCells(inputData.meshSegments[Axis::Y]), TotalCells(inputData.meshSegments[Axis::Z])} ),

    nFacesNormal( { NumberOfFaces( nCells, Axis::X ), NumberOfFaces( nCells, Axis::Y ), NumberOfFaces( nCells, Axis::Z ) } ),

    cellCenters( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellFaces( {{Axis::ENUMDATA::X, nCells(0) + 1},
                {Axis::ENUMDATA::Y, nCells(1) + 1},
                {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellLengths( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellLengthsInv( {{Axis::ENUMDATA::X, nCells(0)},
                     {Axis::ENUMDATA::Y, nCells(1)},
                     {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellCenterDiffInv( {{Axis::ENUMDATA::X, nCells(0) + 1},
                        {Axis::ENUMDATA::Y, nCells(1) + 1},
                        {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    interpFactors( {{Axis::ENUMDATA::X, nCells(0) + 1},
                    {Axis::ENUMDATA::Y, nCells(1) + 1},
                    {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellFaceAreas( {{Axis::ENUMDATA::X, {nCells(1), nCells(2)} },
                    {Axis::ENUMDATA::Y, {nCells(0), nCells(2)} },
                    {Axis::ENUMDATA::Z, {nCells(0), nCells(1)} }} ),

    extrapFactors()

    { 
        std::vector< std::vector<floatType> > growthRates(Axis::count);

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            growthRates[axis] = CalculateGrowthRates(inputData.meshSegments[axis]);

            CalculateCellLengths(cellLengths[axis], inputData.meshSegments[axis], growthRates[axis]);
            cellLengthsInv[axis] = cellLengths[axis].inverse();

            CalculateCellCenters(cellCenters[axis], cellLengths[axis], inputData.meshSegments[axis].front().lowerBound);
            CalculateCellCenterDiffInv(cellCenterDiffInv[axis], cellCenters[axis]);

            CalculateCellFaces(cellFaces[axis], cellLengths[axis], inputData.meshSegments[axis].front().lowerBound);

            switch ( inputData.schemes.faceInterpolationScheme ) {
                case FaceInterpolationSchemes::Average:
                    CalculateInterpolationFactors_Average(interpFactors[axis]);
                    break;

                case FaceInterpolationSchemes::WeightedLinear:
                    CalculateInterpolationFactors_WeightedLinear(interpFactors[axis], cellCenters[axis], cellFaces[axis]);
                    break;
            }
                        
            CalculateExtrapolationFactors(extrapFactors, cellLengths, axis);

        } );


        // Cell face areas should be calculated on their own since they depend on other axis
        EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

            // Axis are ordered by numbering
            Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
            Axis::ENUMDATA axis2 = LUT::HiOrthogonalAxis[ axis ];
            CalculateCellFaceAreas(cellFaceAreas[axis], cellLengths[axis1], cellLengths[axis2]);

        } );

    };



/*-------------------------------------------------------------------------------------*\
                                BoundaryConditionData
\*-------------------------------------------------------------------------------------*/


namespace
{

    // Interpolates user defined profile onto the mesh in 1D. Points that are outside the user defined region are set 
    // to the value of the nearest point.
    Tensor1D InterpProfile1D( const Tensor1D x,
                             const Tensor1D v, 
                             const Tensor1D xquery )
    {
        const intType nValuePoints = x.size(); 
        const intType nQueryPoints = xquery.size(); 
        Tensor1D vquery = Tensor1D( nQueryPoints );

        intType i = 0;
        for ( intType iq = 0; iq != nQueryPoints; iq++ ) {

            bool isBelowMinPoint = ( xquery(iq) <= x(0) );
            if ( isBelowMinPoint ) {
                vquery(iq) = v(0);
                continue;
            }

            bool isAboveMaxPoint = ( xquery(iq) >= x(nValuePoints-1) );
            if ( isAboveMaxPoint ) {
                vquery(iq) = v(nValuePoints-1);
                continue;
            }

            // Internal points
            while ( i != nValuePoints ) {

                bool isBetweenPoints = ( xquery(iq) >= x(i) )  &&  ( xquery(iq) <= x(i+1) );
                if ( isBetweenPoints ) {

                    // Interpolate
                    floatType v0 = v( i ),
                              v1 = v( i + 1 ),
                              x0 = x( i ),
                              x1 = x( i + 1);
                    vquery(iq) = v0 + ( xquery(iq) - x0 ) * ( v1 - v0 ) / ( x1 - x0 );
                    break;
                }
                i++;
            }

        }

        return vquery;
    }



    // Interpolates a 1D profile onto a meshed boundary face
    Tensor2D SetBoundaryProfile1D( const InputData::Profile1D &profile1D,   
                                  const Mesh& mesh,
                                  const BoundaryPatches::ENUMDATA boundaryPatch )
    {
        Axis::ENUMDATA profileAxis  = profile1D.axis,
                    normalAxis   = LUT::BoundaryPatchAxis[ boundaryPatch ],
                    constantAxis = static_cast< Axis::ENUMDATA >( 3 - profileAxis - normalAxis );

        // 1D profile on the mesh points
        Tensor1D interpolatedProfile1D = InterpProfile1D( profile1D.coordinates, profile1D.values, mesh.cellCenters[profileAxis] );

        // Copy out in 2D
        intType nCellsLo = mesh.nCells( LUT::LoOrthogonalAxis[normalAxis] ),
                nCellsHi = mesh.nCells( LUT::HiOrthogonalAxis[normalAxis] );
        Tensor2D boundaryPatchValues( nCellsLo, nCellsHi );
        int constantAxis2D = ( constantAxis == LUT::LoOrthogonalAxis[normalAxis] ) ? 0 : 1; // 3D axis enums cannot be used on the 2D plane
        for ( intType i = 0; i != mesh.nCells(constantAxis); i++ ) {
            boundaryPatchValues.chip(i, constantAxis2D) = interpolatedProfile1D;
        }

        return boundaryPatchValues;
    }



    // Sets a constant value for a boundary face patch
    Tensor2D SetBoundaryProfileConstant( const floatType value,   
                                        const Mesh& mesh,
                                        const BoundaryPatches::ENUMDATA boundaryPatch )
    {
        Axis::ENUMDATA normalAxis   = LUT::BoundaryPatchAxis[ boundaryPatch ];

        intType nCellsLo = mesh.nCells( LUT::LoOrthogonalAxis[normalAxis] ),
                nCellsHi = mesh.nCells( LUT::HiOrthogonalAxis[normalAxis] );

        return Tensor2D( nCellsLo, nCellsHi ).setConstant( value );
    }


}   // end anonymous namespace




FieldData< BoundaryConditionData > SetBoundaryConditionData( const InputData &inputData,
                                                             const Mesh &mesh )
{

    FieldData< EnumVector< BoundaryPatches, BoundaryConditionConfig > > bcData;


    ForAllFieldData( [&] (intType f) {

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

            const auto &bcConfig = inputData.boundaryConditions[f][bp];

            bcData[f][bp].type = bcConfig.type;

            if ( bcConfig.hasUniformValue ) 
                bcData[f][bp].value = SetBoundaryProfileConstant( bcConfig.uniformValue, mesh, bp );

            if ( bcConfig.hasProfile1D ) 
                bcData[f][bp].value = SetBoundaryProfile1D( bcConfig.profile1D, mesh, bp );

        } );

    } );

    return bcData;
}





/*-------------------------------------------------------------------------------------*\
                                    FVCoefficients
\*-------------------------------------------------------------------------------------*/

namespace
{

std::vector< TransportCoefficients::ENUMDATA > PicardEnums( const Axis::ENUMDATA momentumEquation, 
                                                                      const Axis::ENUMDATA velocity)
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Y:
                    return {};
                    break;
                case Z:
                    return {};
                    break;
            }
            break;

        case Y: 

            switch ( velocity ) {
                case X: 
                    return {};
                    break;
                case Y:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Z:
                    return {};
                    break;
            }
            break;

        case Z:

            switch ( velocity ) {
                case X: 
                    return {};
                    break;

                case Y:
                    return {};
                    break;

                case Z:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
            }
            break;
    }
    return {};
}


std::vector< TransportCoefficients::ENUMDATA > NewtonEnums( const Axis::ENUMDATA momentumEquation, 
                                                            const Axis::ENUMDATA velocity)
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Y:
                    return {C::p, C::n, C::s};
                    break;
                case Z:
                    return {C::p, C::t, C::b};
                    break;
            }
            break;

        case Y: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::e, C::w};
                    break;
                case Y:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Z:
                    return {C::p, C::t, C::b};;
                    break;
            }
            break;

        case Z:

            switch ( velocity ) {
                case X: 
                    return {C::p, C::e, C::w};
                    break;

                case Y:
                    return {C::p, C::n, C::s};
                    break;

                case Z:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
            }
            break;
    }
    return {};
}



std::vector< TransportCoefficients::ENUMDATA > MomentumVelocityEnums( const Axis::ENUMDATA momentumEquation, 
                                                                      const Axis::ENUMDATA velocity,
                                                                      Linearisation li )
{
    switch ( li ) {
        case Linearisation::Picard:
            return PicardEnums( momentumEquation, velocity );

        case Linearisation::Newton:
            return NewtonEnums( momentumEquation, velocity );
    }
    return {};
}                                                                      



std::vector< TransportCoefficients::ENUMDATA > MomentumPressureEnums( const Axis::ENUMDATA momentumEquation )
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 
            return { C::p, C::e, C::w };
            break;

        case Y: 
            return { C::p, C::n, C::s };
            break;

        case Z:
            return { C::p, C::t, C::b };
            break;
    }
    return {};
}



std::vector< TransportCoefficients::ENUMDATA > ContinuityPressureEnums( MomentumInterpolation mi ) 
{
    using C = TransportCoefficients::ENUMDATA;
    switch ( mi ) {
        case MomentumInterpolation::Implicit:
            return {C::p, C::n, C::e, C::s, C::w, C::t, C::b, C::nn, C::ee, C::ss, C::ww, C::tt, C::bb};  

        case MomentumInterpolation::SemiExplicit:
            return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
    }
    return {};
}




}   // end anonymous namespace


using C = TransportCoefficients::ENUMDATA;
using enum Axis::ENUMDATA;

// Momentum equations constructor
MomentumEquation::MomentumEquation( const Axis::ENUMDATA axis, 
                                    const iArray3 &dims,
                                    Linearisation li ) :
    AU( { EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, X, li), dims),
          EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, Y, li), dims),
          EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, Z, li), dims) }  ),
    AP( MomentumPressureEnums( axis ), dims( axis ) ),
    B( CFD::Tensor3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    diagCoeffInv( CFD::Tensor3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    diff({ EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims(X) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims(Y) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims(Z) ) }),
    diffBoundary( 0.0f ),
    BUBoundary(),
    BPBoundary(),   // These should be dimensioned only if needed
    relaxation( 1.0f ),
    component( axis ),
    linearisation( li )
{};


// Continuity equations constructor
ContinuityEquation::ContinuityEquation( const iArray3 &dims,
                                        MomentumInterpolation mi ) :
    AU( { EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims( X )),
          EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims( Y )),
          EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims( Z )) } ),
    AP( ContinuityPressureEnums( mi ), dims ),
    B( Tensor3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    mwiSparseCoeffs( { std::array<Tensor1D, 4>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    mwiCompactCoeffs( { std::array<Tensor1D, 2>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    BUBoundary(),
    BPBoundary(),   // These should be dimensioned only if needed
    relaxation( 1.0f ),
    momentumInterpolation( mi )
{};


// Coefficients class constructor
FVCoefficients::FVCoefficients( const iArray3 &dims, 
                                Linearisation li,
                                MomentumInterpolation mi ) :
    Mom( { MomentumEquation(X, dims, li),  MomentumEquation(Y, dims, li),  MomentumEquation(Z, dims, li) } ),
    Cont( dims, mi ),
    nCells( dims )
{};



/*-------------------------------------------------------------------------------------*\
                                      Free Functions
\*-------------------------------------------------------------------------------------*/


FieldData<Tensor3D> InitialiseFields( const Mesh &mesh, 
                                     const InputData &inputData )
{

    FieldData<Tensor3D> fields( Tensor3D( mesh.nCells(0) + 2*CFD::nGhost, mesh.nCells(1) + 2*CFD::nGhost, mesh.nCells(2) + 2*CFD::nGhost).setZero() );

    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                 extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    // Set initial values
    ForAllFieldData( [&] (intType i) { 
        fields[i].slice( offsets, extents ).setConstant( inputData.initialConditions[i] );  
    } );
    
    return fields;
}


}   // end namespace CFD
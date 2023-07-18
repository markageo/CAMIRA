#include "FiniteVolume.h"

#include <cmath>

#include <iostream>

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


void CalculateCellLengths(array1D &cellLengths, 
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


void CalculateCellCenters(array1D &cellCenters, 
                          const array1D &cellLengths,
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


void CalculateCellCenterDiffInv(array1D &cellCenterDiffInv, 
                                const array1D &cellCenters)
{
    // First and last element dont correspond to valid values
    intType nFaces = cellCenters.size() + 1;

    for (intType i = 1; i != nFaces-1; i++) {
        cellCenterDiffInv(i) = 1.0f/( cellCenters(i) - cellCenters(i-1) );
    }
}


void CalculateCellFaces(array1D &cellFaces, 
                        const array1D &cellLengths,
                        const floatType startPosition)
{
    intType nFaces = cellFaces.size();

    cellFaces(0) = startPosition;
    for (intType i = 1; i != nFaces; i++) {
        cellFaces(i) = cellFaces(i-1) + cellLengths(i-1);
    }
}


void CalculateCellFaceAreas(array2D &cellFaceAreas, 
                            const array1D &cellLengths_x, 
                            const array1D &cellLengths_y)
{
    for (int j = 0; j != cellLengths_y.dimension(0); j++) {
        for (int i = 0; i != cellLengths_x.dimension(0); i++) {
            cellFaceAreas(i, j) = cellLengths_x(i) * cellLengths_y(j);
        }
    }
}


void CalculateInterpolationFactors_WeightedLinear( array1D &interpFactors, 
                                                   const array1D &cellCenters, 
                                                   const array1D &cellFaces) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) 
                         / ( cellCenters(i) - cellCenters(i-1) );
    }
}


void CalculateInterpolationFactors_Average( array1D &interpFactors ) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = 0.5f;
    }
}


Mesh::ExtrapFactorsStruct GetExtrapolationFactors(const array1D &cellLengths, 
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
                                   const EnumVector<Axis, array1D> &cellLengths, 
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

}   // end anonymous namespace


// Constructor, creates the mesh
Mesh::Mesh(const InputData &inputData) :
    nCells( { TotalCells(inputData.meshSegments[Axis::X]),  TotalCells(inputData.meshSegments[Axis::Y]), TotalCells(inputData.meshSegments[Axis::Z])} ),

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
            Axis::ENUMDATA axis1 = ( axis == Axis::X ) ? Axis::Y : Axis::X;
            Axis::ENUMDATA axis2 = ( axis == Axis::Z ) ? Axis::Y : Axis::Z;
            CalculateCellFaceAreas(cellFaceAreas[axis], cellLengths[axis1], cellLengths[axis2]);

        } );

    };




/*-------------------------------------------------------------------------------------*\
                                    FVCoefficients
\*-------------------------------------------------------------------------------------*/

namespace
{

std::vector< TransportCoefficients::ENUMDATA > MomentumVelocityEnums( const Axis::ENUMDATA momentumEquation, 
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


// Pressure coefficients are different depending on momentum interpolation treatement
template< MomentumInterpolation MI >
std::vector< TransportCoefficients::ENUMDATA > ContinuityPressureEnums() = delete;

template<>
std::vector< TransportCoefficients::ENUMDATA > ContinuityPressureEnums< MomentumInterpolation::Implicit >() {
    using C = TransportCoefficients::ENUMDATA;
    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b, C::nn, C::ee, C::ss, C::ww, C::tt, C::bb};
}

template<>
std::vector< TransportCoefficients::ENUMDATA > ContinuityPressureEnums< MomentumInterpolation::SemiExplicit >() {
    using C = TransportCoefficients::ENUMDATA;
    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
}



}   // end anonymous namespace


using C = TransportCoefficients::ENUMDATA;
using enum Axis::ENUMDATA;

// Momentum equations constructor
MomentumEquation::MomentumEquation( const Axis::ENUMDATA axis, 
                                                    const iVector3 &dims) :
    AU( { EnumVector<TransportCoefficients, array3D>( MomentumVelocityEnums(axis, X), dims),
          EnumVector<TransportCoefficients, array3D>( MomentumVelocityEnums(axis, Y), dims),
          EnumVector<TransportCoefficients, array3D>( MomentumVelocityEnums(axis, Z), dims) }  ),
    AP( MomentumPressureEnums( axis ), dims( axis ) ),
    B( CFD::array3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    diagCoeffInv( CFD::array3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    diff({ EnumVector<TransportCoefficients, array1D>( {C::p, C::e, C::w}, dims(X) ),
           EnumVector<TransportCoefficients, array1D>( {C::p, C::n, C::s}, dims(Y) ),
           EnumVector<TransportCoefficients, array1D>( {C::p, C::t, C::b}, dims(Z) ) }),
    boundaryDiff( 0.0f ),
    boundaryP( 0.0f ),
    boundaryVel( {array2D( dims(Y), dims(Z) ).setZero(), array2D( dims(Y), dims(Z) ).setZero(),
                  array2D( dims(X), dims(Z) ).setZero(), array2D( dims(X), dims(Z) ).setZero(),
                  array2D( dims(X), dims(Y) ).setZero(), array2D( dims(X), dims(Y) ).setZero()} ),    // This doesn't follow right hand rule
    relaxation( 1.0f )
{};


// Continuity equations constructor
template< MomentumInterpolation MI >
ContinuityEquation<MI>::ContinuityEquation( const iVector3 &dims ) :
    AU( { EnumVector<TransportCoefficients, array1D>( {C::p, C::e, C::w}, dims( X )),
          EnumVector<TransportCoefficients, array1D>( {C::p, C::n, C::s}, dims( Y )),
          EnumVector<TransportCoefficients, array1D>( {C::p, C::t, C::b}, dims( Z )) } ),
    AP( ContinuityPressureEnums<MI>(), dims ),
    B( array3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    mwiSparseCoeffs( { std::array<array1D, 4>{ array1D(dims(X)+1).setZero(), array1D(dims(X)+1).setZero(), array1D(dims(X)+1).setZero(), array1D(dims(X)+1).setZero() } ,
                       std::array<array1D, 4>{ array1D(dims(Y)+1).setZero(), array1D(dims(Y)+1).setZero(), array1D(dims(Y)+1).setZero(), array1D(dims(Y)+1).setZero() } ,
                       std::array<array1D, 4>{ array1D(dims(Z)+1).setZero(), array1D(dims(Z)+1).setZero(), array1D(dims(Z)+1).setZero(), array1D(dims(Z)+1).setZero() } } ),
    mwiCompactCoeffs( { std::array<array1D, 2>{ array1D(dims(X)+1).setZero(), array1D(dims(X)+1).setZero() } ,
                        std::array<array1D, 2>{ array1D(dims(Y)+1).setZero(), array1D(dims(Y)+1).setZero() } ,
                        std::array<array1D, 2>{ array1D(dims(Z)+1).setZero(), array1D(dims(Z)+1).setZero() } } ),
    boundaryP( {array2D( dims(Y), dims(Z) ).setZero(), array2D( dims(Y), dims(Z) ).setZero(),
                array2D( dims(X), dims(Z) ).setZero(), array2D( dims(X), dims(Z) ).setZero(),
                array2D( dims(X), dims(Y) ).setZero(), array2D( dims(X), dims(Y) ).setZero()} ),     // This doesn't follow right hand rule
    boundaryVel( 0.0f ),
    relaxation( 1.0f )
{};
template struct ContinuityEquation< MomentumInterpolation::Implicit >;
template struct ContinuityEquation< MomentumInterpolation::SemiExplicit >;


// Coefficients class constructor
template< MomentumInterpolation MI >
FVCoefficients<MI>::FVCoefficients( const iVector3 &dims ) :
    Mom( { MomentumEquation(X, dims),  MomentumEquation(Y, dims),  MomentumEquation(Z, dims) } ),
    Cont( dims ),
    nCells( dims )
{};
template struct FVCoefficients< MomentumInterpolation::Implicit >;
template struct FVCoefficients< MomentumInterpolation::SemiExplicit >;


/*-------------------------------------------------------------------------------------*\
                                      Free Functions
\*-------------------------------------------------------------------------------------*/


void RemoveGhostCells( array3D &array, 
                       const intType nGhostCells)
{
    arrayIndex3D offsets = { nGhostCells, nGhostCells, nGhostCells },
                 extents = { array.dimension(0) - 2*nGhostCells, 
                             array.dimension(1) - 2*nGhostCells,
                             array.dimension(2) - 2*nGhostCells };
                                              
    array = array3D( array ).slice(offsets, extents);
}


FieldData<array3D> InitialiseFields( const Mesh &mesh, 
                                     const InputData &inputData )
{

    FieldData<array3D> fields( array3D( mesh.nCells(0) + 2*CFD::nGhost, mesh.nCells(1) + 2*CFD::nGhost, mesh.nCells(2) + 2*CFD::nGhost).setZero() );

    arrayIndex3D offsets = {nGhost, nGhost, nGhost},
                 extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    // Set initial values
    ForAllFieldData( [&] (intType i) { 
        fields[i].slice( offsets, extents ).setConstant( inputData.initialConditions[i] );  
    } );
    
    return fields;
}


}   // end namespace CFD
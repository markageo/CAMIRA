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


void CalculateInterpolationFactors(array1D &interpFactors, 
                                   const array1D &cellCenters, 
                                   const array1D &cellFaces) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) 
                         / ( cellCenters(i) - cellCenters(i-1) );
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
                                   const ArrayAllocator<Axis, array1D> &cellLengths, 
                                   const Axis::ENUMDATA axis)
{  

    using enum BoundaryPatches::ENUMDATA;
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA patchPositive = PositivePatch[ axis ], 
                              patchNegative = NegativePatch[ axis ];

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
            CalculateInterpolationFactors(interpFactors[axis], cellCenters[axis], cellFaces[axis]);
            CalculateExtrapolationFactors(extrapFactors, cellLengths, axis);

        } );

    };


/*-------------------------------------------------------------------------------------*\
                                    FVCoefficients
\*-------------------------------------------------------------------------------------*/

namespace
{

// Return the vector of enums corresponding to the coefficients required to multiply be a variable for a particular equation
std::vector< TransportCoefficients::ENUMDATA > EquationEnums(const Fields::ENUMDATA equation, 
                                                             const Fields::ENUMDATA variable)
{
    using F = Fields::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;
    switch (equation) {
        case F::U:  // U momentum equation
            if        (variable == F::U) {
                return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
            } else if (variable == F::V) {
                return {};
            } else if (variable == F::W) {
                return {};
            } else if (variable == F::P) {
                return {C::p, C::e, C::w};
            }
            break;


        case F::V:  // V momentum equation
            if        (variable == F::U) {
                return {};
            } else if (variable == F::V) {
                return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
            } else if (variable == F::W) {
                return {};
            } else if (variable == F::P) {
                return {C::p, C::n, C::s};
            }
            break;

        case F::W:  // W momentum equation
            if        (variable == F::U) {
                return {};
            } else if (variable == F::V) {
                return {};
            } else if (variable == F::W) {
                return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
            } else if (variable == F::P) {
                return {C::p, C::t, C::b};
            }
            break;

        case F::P:  // Continuity equation (not Poisson pressure equation)
            if        (variable == F::U) {
                return {C::p, C::e, C::w};
            } else if (variable == F::V) {
                return {C::p, C::n, C::s};
            } else if (variable == F::W) {
                return {C::p, C::t, C::b};
            } else if (variable == F::P) {
                return {C::p, C::n, C::e, C::s, C::w, C::t, C::b, 
                        C::nn, C::ee, C::ss, C::ww, C::tt, C::bb};
            }
            break;
    }
    return {};
}

// Returns the dimension axis corresponding to the U, V, or W direction
Axis::ENUMDATA EquationDim(const Fields::ENUMDATA field) {
    using F = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    if (field == F::U) {
        return X;
    } else if (field == F::V) {
        return Y;
    } else if (field == F::W) {
        return Z;
    }
    return X;
}

}   // end anonymous namespace


using C = TransportCoefficients::ENUMDATA;
using F = Fields::ENUMDATA;
using enum Axis::ENUMDATA;

// Momentum equations constructor
FVCoefficients::MomentumEquation::MomentumEquation( const Fields::ENUMDATA field, 
                                                    const indexVector3 &dims) :
    AU( EquationEnums(field, F::U), dims ),
    AV( EquationEnums(field, F::V), dims ),
    AW( EquationEnums(field, F::W), dims ),
    AP( EquationEnums(field, F::P), dims( EquationDim(field) ) ),
    B( CFD::array3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    diff({ ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::e, C::w}, dims(X) ),
           ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::n, C::s}, dims(Y) ),
           ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::t, C::b}, dims(Z) ) }),
    boundaryDiff( 0.0f ),
    boundaryP( 0.0f ),
    boundaryVel( {array2D( dims(Y), dims(Z) ).setZero(), array2D( dims(Y), dims(Z) ).setZero(),
                  array2D( dims(X), dims(Z) ).setZero(), array2D( dims(X), dims(Z) ).setZero(),
                  array2D( dims(X), dims(Y) ).setZero(), array2D( dims(X), dims(Y) ).setZero()} )    // This doesn't follow right hand rule
{};


// Continuity equations constructor
FVCoefficients::ContinuityEquation::ContinuityEquation( const indexVector3 &dims ) :
    AU( EquationEnums(F::P, F::U), dims( EquationDim(F::U) ) ),
    AV( EquationEnums(F::P, F::V), dims( EquationDim(F::V) ) ),
    AW( EquationEnums(F::P, F::W), dims( EquationDim(F::W) ) ),
    AP( EquationEnums(F::P, F::P), dims ),
    B( CFD::array3D( dims(X), dims(Y), dims(Z) ).setZero() ),
    boundaryP( {array2D( dims(Y), dims(Z) ).setZero(), array2D( dims(Y), dims(Z) ).setZero(),
                array2D( dims(X), dims(Z) ).setZero(), array2D( dims(X), dims(Z) ).setZero(),
                array2D( dims(X), dims(Y) ).setZero(), array2D( dims(X), dims(Y) ).setZero()} ),     // This doesn't follow right hand rule
    boundaryVel( 0.0f )
{};


// Coefficients class constructor
FVCoefficients::FVCoefficients(const indexVector3 &dims) :
    Umom(F::U, dims),
    Vmom(F::V, dims),
    Wmom(F::W, dims),
    Cont(dims),
    nCells( dims )
{};


/*-------------------------------------------------------------------------------------*\
                                      Free Functions
\*-------------------------------------------------------------------------------------*/



void RemoveGhostCells( array3D &array, 
                       const intType nGhostCells)
{
    Eigen::array<Eigen::Index, 3> offsets = { nGhostCells, nGhostCells, nGhostCells },
                                  extents = { array.dimension(0) - 2*nGhostCells, 
                                              array.dimension(1) - 2*nGhostCells,
                                              array.dimension(2) - 2*nGhostCells };
                                              
    array = array3D( array ).slice(offsets, extents);
}


ArrayAllocator<Fields, array3D> InitialiseFields( const Mesh &mesh, 
                                                  const InputData &inputData )
{
    // Create fields
    ArrayAllocator<Fields, array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells + 2*CFD::nGhost);

    Eigen::array<Eigen::Index, 3> offsets = {nGhost, nGhost, nGhost},
                                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    // Set initial values
    EnumFor<Fields>( [&] (Fields::ENUMDATA f) {
        fields[f].slice( offsets, extents ).setConstant( inputData.initialConditions[f] );  // Don't set the ghost cells
    } );

    return fields;
}


}   // end namespace CFD
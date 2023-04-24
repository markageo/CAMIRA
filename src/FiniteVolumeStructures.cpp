#include "FiniteVolumeStructures.h"

#include <cmath>

#include <iostream>

using namespace CFD;

/*-------------------------------------------------------------------------------------*\
                                         Mesh
\*-------------------------------------------------------------------------------------*/

namespace
{

std::vector<floatType> CalculateGrowthRates(const std::vector<InputData::MeshSegment> &meshSegments)
{
    int nSegments = meshSegments.size();
    std::vector<floatType> growthRates(nSegments);

    // Negative growth rate means shrinking grid
    for (int i = 0; i != nSegments; i++) {
        growthRates[i] = std::pow( std::abs(meshSegments[i].biasFactor) , 1.0/( meshSegments[i].nCells - 1 ) );
        if (meshSegments[i].biasFactor < 0 )   
            growthRates[i] = 1.0/growthRates[i];
    }

    return growthRates;
}


void CalculateCellLengths(array1D &cellLengths, 
                          const std::vector<InputData::MeshSegment> &meshSegments, 
                          const std::vector<floatType> &growthRates)
{
    int nSegments = meshSegments.size();

    floatType segmentLength, firstCellLength, geometricFactor;
    int cellIndex = 0;
    for (int s = 0; s != nSegments; s++) {    // Segments

        if (growthRates[s] != 1.0) {
            geometricFactor = (1.0 - std::pow( growthRates[s], meshSegments[s].nCells )) / (1.0 - growthRates[s]);   // geometric series formula
        } else {
            geometricFactor = meshSegments[s].nCells;
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
                          const array1D &cellLengths)
{
    int nCellsTotal = cellLengths.size();

    floatType previousCellPosition = 0.0, previousCellLength = 0.0;
    for (int i = 0; i != nCellsTotal; i++) {
        cellCenters(i) = previousCellPosition + previousCellLength/2.0 + cellLengths(i)/2.0;
        previousCellPosition = cellCenters(i);
        previousCellLength = cellLengths(i);
    }

}


void CalculateCellCenterDiffInv(array1D &cellCenterDiffInv, 
                                const array1D &cellCenters)
{
    // First and last element dont correspond to valid values
    int nFaces = cellCenters.size() + 1;

    for (int i = 1; i != nFaces-1; i++) {
        cellCenterDiffInv(i) = 1.0/( cellCenters(i) - cellCenters(i-1) );
    }
}


void CalculateCellFaces(array1D &cellFaces, 
                        const array1D &cellLengths)
{
    int nFaces = cellLengths.size() + 1;

    cellFaces(0) = 0;
    for (int i = 1; i != nFaces; i++) {
        cellFaces(i) = cellFaces(i-1) + cellLengths(i-1);
    }
}


void CalculateInterpolationFactors(array1D &interpFactors, 
                                   const array1D &cellCenters, 
                                   const array1D &cellFaces) 
{
    for (int i = 1; i != interpFactors.dimension(0)-1; i++) {
        interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) / ( cellCenters(i) - cellCenters(i-1) );
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


Mesh::ExtrapFactorsStruct GetExtrapolationFactors(const array1D &cellLengths, 
                                                  const int fieldIndex_p, 
                                                  const int fieldIndex_a)
{
    floatType extrapFactor_p = ( 2.0*cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );

    floatType extrapFactor_a = - ( cellLengths(fieldIndex_p) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );
    
    return Mesh::ExtrapFactorsStruct( extrapFactor_p, extrapFactor_a );
}


void CalculateExtrapolationFactors(EnumVector<BoundaryPatches, Mesh::ExtrapFactorsStruct > &extrapFactors, 
                                   const ArrayAllocator<Axis, array1D> &cellLengths, 
                                   const Axis::ENUMDATA axis)
{  

    using enum BoundaryPatches::ENUMDATA;
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA patchPositive, patchNegative;
    if         (axis == X) {
        patchPositive = xPositive;
        patchNegative = xNegative;
    } else if  (axis == Y) {
        patchPositive = yPositive;
        patchNegative = yNegative;
    } else if  (axis == Z) {
        patchPositive = zPositive;
        patchNegative = zNegative;
    }

    int fieldIndex_p, fieldIndex_a; // Boundary cell node and the adjacent one

    // Positive patch boundary
    fieldIndex_p = cellLengths[axis].dimension(0) - 1;
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
                    {Axis::ENUMDATA::Y, {nCells(2), nCells(0)} },
                    {Axis::ENUMDATA::Z, {nCells(0), nCells(1)} }} ),

    extrapFactors()

    { 
        using enum Axis::ENUMDATA;

        std::vector< std::vector<floatType> > growthRates(Axis::count);

        for (int i = 0; i != Axis::count; i++) {
            Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>(i);
            growthRates[axis] = CalculateGrowthRates(inputData.meshSegments[axis]);

            CalculateCellLengths(cellLengths[axis], inputData.meshSegments[axis], growthRates[axis]);
            cellLengthsInv[axis] = cellLengths[axis].inverse();

            CalculateCellCenters(cellCenters[axis], cellLengths[axis]);
            CalculateCellCenterDiffInv(cellCenterDiffInv[axis], cellCenters[axis]);

            CalculateCellFaces(cellFaces[axis], cellLengths[axis]);
            CalculateInterpolationFactors(interpFactors[axis], cellCenters[axis], cellFaces[axis]);
            CalculateExtrapolationFactors(extrapFactors, cellLengths, axis);
        }

        // Cell face areas should be calculated on their own since they depend o other axis
        for (int i = 0; i != Axis::count; i++) {
            Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>(i);

            // Cyclic permutations to get correct order of axis for area, which is indexed by right hand rule
            Axis::ENUMDATA axis1 = static_cast<Axis::ENUMDATA>( (i+1) % Axis::count );
            Axis::ENUMDATA axis2 = static_cast<Axis::ENUMDATA>( (i+2) % Axis::count );
            CalculateCellFaceAreas(cellFaceAreas[axis], cellLengths[axis1], cellLengths[axis2]);
        }

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
intType EquationDim(const Fields::ENUMDATA field) {
    using F = Fields::ENUMDATA;
    if (field == F::U) {
        return 0;
    } else if (field == F::V) {
        return 1;
    } else if (field == F::W) {
        return 2;
    }
    return -1;
}

}   // end anonymous namespace


using C = TransportCoefficients::ENUMDATA;
using F = Fields::ENUMDATA;
using enum Axis::ENUMDATA;

// Momentum equations constructor
FVCoefficients::MomentumEquation::MomentumEquation(const Fields::ENUMDATA field, 
                                                   const CFD::indexVector3 &dims) :
    AU( EquationEnums(field, F::U), dims ),
    AV( EquationEnums(field, F::V), dims ),
    AW( EquationEnums(field, F::W), dims ),
    AP( EquationEnums(field, F::P), dims( EquationDim(field) ) ),
    B( dims(X), dims(Y), dims(Z) ),
    diff({ ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::e, C::w}, dims(X) ),
           ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::n, C::s}, dims(Y) ),
           ArrayAllocator<TransportCoefficients, array1D>( {C::p, C::t, C::b}, dims(Z) ) }),
    boundaryDiff(),
    boundaryP(),
    boundaryVel( {array2D( dims(Y), dims(Z) ), array2D( dims(Y), dims(Z) ),
                  array2D( dims(X), dims(Z) ), array2D( dims(X), dims(Z) ),
                  array2D( dims(X), dims(Y) ), array2D( dims(X), dims(Y) )} )   // This doesnt follow right hand rule
{};


// Continuity equations constructor
FVCoefficients::ContinuityEquation::ContinuityEquation(const indexVector3 &dims) :
    AU( EquationEnums(F::P, F::U), dims( EquationDim(F::U) ) ),
    AV( EquationEnums(F::P, F::V), dims( EquationDim(F::V) ) ),
    AW( EquationEnums(F::P, F::W), dims( EquationDim(F::W) ) ),
    AP( EquationEnums(F::P, F::P), dims ),
    B( dims(X), dims(Y), dims(Z) ),
    boundaryVel(),
    boundaryP()
{};


// Coefficients class constructor
FVCoefficients::FVCoefficients(const indexVector3 &dims) :
    Umom(F::U, dims),
    Vmom(F::V, dims),
    Wmom(F::W, dims),
    Cont(dims)
{};


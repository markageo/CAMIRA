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
    intType nFaces = cellLengths.size() + 1;

    cellFaces(0) = startPosition;
    for (intType i = 1; i != nFaces; i++) {
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

    BoundaryPatches::ENUMDATA patchPositive = xPositive, 
                              patchNegative = xNegative;
    if         (axis == X) {
        patchPositive = xPositive;
        patchNegative = xNegative;
    } else if  (axis == Y) {
        patchPositive = yPositive;
        patchNegative = yNegative;
    } else if  (axis == Z) {
        patchPositive = zPositive;
        patchNegative = zNegative;
    } else {
        /* NULL */
    }

    intType fieldIndex_p, fieldIndex_a; // Boundary cell node and the adjacent one

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

    extrapFactors()

    { 
        using enum Axis::ENUMDATA;

        std::vector< std::vector<floatType> > growthRates(Axis::count);

        for (int i = 0; i != Axis::count; i++) {
            Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>(i);
            growthRates[axis] = CalculateGrowthRates(inputData.meshSegments[axis]);

            CalculateCellLengths(cellLengths[axis], inputData.meshSegments[axis], growthRates[axis]);
            cellLengthsInv[axis] = cellLengths[axis].inverse();

            CalculateCellCenters(cellCenters[axis], cellLengths[axis], inputData.meshSegments[axis].front().lowerBound);
            CalculateCellCenterDiffInv(cellCenterDiffInv[axis], cellCenters[axis]);

            CalculateCellFaces(cellFaces[axis], cellLengths[axis], inputData.meshSegments[axis].front().lowerBound);
            CalculateInterpolationFactors(interpFactors[axis], cellCenters[axis], cellFaces[axis]);
            CalculateExtrapolationFactors(extrapFactors, cellLengths, axis);
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
    Cont(dims)
{};


/*-------------------------------------------------------------------------------------*\
                                      Free Functions
\*-------------------------------------------------------------------------------------*/


namespace
{

    void SwapMesh(Mesh &mesh, const Axis::ENUMDATA axis1, const Axis::ENUMDATA axis2)
    {
        std::swap( mesh.nCells(axis1)           , mesh.nCells(axis2) );
        std::swap( mesh.cellCenters[axis1]      , mesh.cellCenters[axis2] );
        std::swap( mesh.cellFaces[axis1]        , mesh.cellFaces[axis2] );
        std::swap( mesh.cellLengths[axis1]      , mesh.cellLengths[axis2] );
        std::swap( mesh.cellLengthsInv[axis1]   , mesh.cellLengthsInv[axis2] );
        std::swap( mesh.cellCenterDiffInv[axis1], mesh.cellCenterDiffInv[axis2] );
        std::swap( mesh.interpFactors[axis1]    , mesh.interpFactors[axis2] );
        std::swap( mesh.interpFactors[axis1]    , mesh.interpFactors[axis2] );

        std::swap( mesh.extrapFactors[ PositivePatch[axis1] ] , mesh.extrapFactors[ PositivePatch[axis2] ] );
        std::swap( mesh.extrapFactors[ NegativePatch[axis1] ] , mesh.extrapFactors[ NegativePatch[axis2] ] );
    }


    // Returns a copy of a 1D array that has been reversed 
    array1D ReversedArray( array1D &array )
    {
        Eigen::array<bool, 1> rev({true});
        return array.reverse( rev );
    };


    void ReverseMeshAxis(Mesh &mesh, const Axis::ENUMDATA axis)
    {
        // Eigen .reverse is done in place, so need to made a copy!
        mesh.cellCenters[axis]       = - ReversedArray( mesh.cellCenters[axis] );
        mesh.cellFaces[axis]         = - ReversedArray( mesh.cellFaces[axis] );
        mesh.cellLengths[axis]       = ReversedArray( mesh.cellLengths[axis] );
        mesh.cellLengthsInv[axis]    = ReversedArray( mesh.cellLengthsInv[axis] );
        mesh.cellCenterDiffInv[axis] = ReversedArray( mesh.cellCenterDiffInv[axis] );
        mesh.interpFactors[axis]     = 1 - ReversedArray( mesh.interpFactors[axis] );
        std::swap( mesh.extrapFactors[ PositivePatch[axis] ], mesh.extrapFactors[ NegativePatch[axis] ] );
    }


}   // end anonymous namespace


void TransformToUserCoordinates( Mesh &mesh, 
                                 ArrayAllocator<Fields, array3D> &fields, 
                                 const InputData::AxisTransformationMap &axisTransformation)                                  
{
    Axis::ENUMDATA codeAxis, userAxis;
    BoundaryPatches::ENUMDATA codePositivePatch, userPatchFromPositive;
    std::array<Axis::ENUMDATA, Axis::count> shuffleOrder = {X, Y, Z}; // Tracks where the user axis have been shuffled to when transforming mesh
    Eigen::array<intType , Axis::count> shuffleArray;
    Eigen::array<bool, Axis::count> reverseArray;

    // Mesh
    for (int a = 0; a != Axis::count; a++) {
        codeAxis = static_cast<Axis::ENUMDATA>(a);
        codePositivePatch = PositivePatch[codeAxis];

        userPatchFromPositive = axisTransformation.UserPatch(codePositivePatch);
        userAxis = BoundaryPatchAxis[userPatchFromPositive];

        if ( shuffleOrder[userAxis] != codeAxis ) {     // Only if it needs to be swapped
            std::swap( shuffleOrder[codeAxis], shuffleOrder[userAxis]  );
            SwapMesh(mesh, codeAxis, userAxis);
        }

        if ( userPatchFromPositive == NegativePatch[userAxis] ) {
            ReverseMeshAxis(mesh, userAxis);
        }

        // Fill the shuffle and revsere arrays, used for 3D arrays
        shuffleArray[codeAxis] = BoundaryPatchAxis[ axisTransformation.CodePatch( codePositivePatch ) ];
        reverseArray[codeAxis] = axisTransformation.CodePatch( codePositivePatch ) == NegativePatch[ static_cast<size_t>( shuffleArray[codeAxis] ) ];    // This reverse is after the shuffling
    }


    // 3D arrays
    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<F>(f);

        // Cell center values
        fields[field] = array3D( fields[field] ).shuffle(shuffleArray).reverse(reverseArray);   // Have to make a copy

        // if (field == F::P) 
        //     continue;

        // // Face velocities
        // faceVelocities[field] = array3D( faceVelocities[field] ).shuffle(shuffleArray).reverse(reverseArray);    // Have to make a copy

    }

}





ArrayAllocator<Fields, array3D> InitialiseFields(const Mesh &mesh, const InputData &inputData)
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
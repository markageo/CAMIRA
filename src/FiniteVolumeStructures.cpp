#include "FiniteVolumeStructures.h"

#include <cmath>

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


    void CalculateCellLengths(array1D &cellLengths, const std::vector<InputData::MeshSegment> &meshSegments, const std::vector<floatType> &growthRates)
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
            firstCellLength = segmentLength/geometricFactor; 

            for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
                cellLengths( cellIndex ) = firstCellLength*std::pow( growthRates[s], static_cast<floatType>(i) );
                cellIndex++;
            }
        }

    }


    void CalculateCellCenters(array1D &cellCenters, const array1D &cellLengths)
    {
        int nCellsTotal = cellLengths.size();

        floatType previousCellPosition = 0.0, previousCellLength = 0.0;
        for (int i = 0; i != nCellsTotal; i++) {
            cellCenters(i) = previousCellPosition + previousCellLength/2.0 + cellLengths(i)/2.0;
            previousCellPosition = cellCenters(i);
            previousCellLength = cellLengths(i);
        }

    }


    void CalculateCellFaces(array1D &cellFaces, const array1D &cellLengths)
    {
        int nFaces = cellLengths.size() + 1;

        floatType previousFacePosition = 0.0, previousCellLength = 0.0;
        for (int i = 0; i != nFaces; i++) {
            cellFaces(i) = previousFacePosition + previousCellLength;
            previousFacePosition = cellFaces(i);
            previousCellLength = cellLengths(i);
        }

    }


    void CalculateInterpolationFactors(array1D &interpFactors, const array1D &cellCenters, const array1D &cellFaces) 
    {
        for (int i = 1; i != interpFactors.dimension(0)-1; i++) {
            interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) / ( cellCenters(i) - cellCenters(i-1) );
        }
    }


    void CalculateCellFaceAreas(array2D &cellFaceAreas, const array1D &cellLengths_x, const array1D &cellLengths_y)
    {
        for (int j = 0; j != cellLengths_y.dimension(0); j++) {
            for (int i = 0; i != cellLengths_x.dimension(0); i++) {
                cellFaceAreas(i, j) = cellLengths_x(i) * cellLengths_y(j);
            }
        }
    }


    Mesh::ExtrapFactorsStruct GetExtrapolationFactors(const array1D &cellLengths, const int fieldIndex_p, const int fieldIndex_a)
    {
        floatType extrapFactor_p = ( 2.0*cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) )
                                 / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );

        floatType extrapFactor_a = - ( cellLengths(fieldIndex_p) )
                                   / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );
        
        return Mesh::ExtrapFactorsStruct( extrapFactor_p, extrapFactor_a );
    }


    void CalculateExtrapolationFactors(std::vector< Mesh::ExtrapFactorsStruct > &extrapFactors, const ArrayAllocator<Axis, array1D> &cellLengths, const Axis::ENUMDATA axis)
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

}

// Constructor, creates the mesh
Mesh::Mesh(const InputData &inputData) :
    nCells( { TotalCells(inputData.meshSegments_x),  TotalCells(inputData.meshSegments_y), TotalCells(inputData.meshSegments_z)} ),

    cellCenters( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellFaces( {{Axis::ENUMDATA::X, nCells(0) + 1},
                {Axis::ENUMDATA::Y, nCells(1) + 1},
                {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellLengths( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    interpFactors( {{Axis::ENUMDATA::X, nCells(0) + 1},
                    {Axis::ENUMDATA::Y, nCells(1) + 1},
                    {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellFaceAreas( {{Axis::ENUMDATA::X, {nCells(1), nCells(2)} },
                    {Axis::ENUMDATA::Y, {nCells(2), nCells(0)} },
                    {Axis::ENUMDATA::Z, {nCells(0), nCells(1)} }} ),

    extrapFactors( BoundaryPatches::count )

    { 
        std::vector<floatType> growthRates_x = CalculateGrowthRates(inputData.meshSegments_x);
        std::vector<floatType> growthRates_y = CalculateGrowthRates(inputData.meshSegments_y);
        std::vector<floatType> growthRates_z = CalculateGrowthRates(inputData.meshSegments_z);

        using enum Axis::ENUMDATA;

        CalculateCellLengths(cellLengths[X], inputData.meshSegments_x, growthRates_x);
        CalculateCellLengths(cellLengths[Y], inputData.meshSegments_y, growthRates_y);
        CalculateCellLengths(cellLengths[Z], inputData.meshSegments_z, growthRates_z);

        CalculateCellCenters(cellCenters[X], cellLengths[X]);
        CalculateCellCenters(cellCenters[Y], cellLengths[Y]);
        CalculateCellCenters(cellCenters[Z], cellLengths[Z]);

        CalculateCellFaces(cellFaces[X], cellLengths[X]);
        CalculateCellFaces(cellFaces[Y], cellLengths[Y]);
        CalculateCellFaces(cellFaces[Z], cellLengths[Z]);

        CalculateCellFaceAreas(cellFaceAreas[X], cellLengths[Y], cellLengths[Z]);
        CalculateCellFaceAreas(cellFaceAreas[Y], cellLengths[Z], cellLengths[X]);
        CalculateCellFaceAreas(cellFaceAreas[Z], cellLengths[X], cellLengths[Y]);

        CalculateInterpolationFactors(interpFactors[X], cellCenters[X], cellFaces[X]);
        CalculateInterpolationFactors(interpFactors[Y], cellCenters[Y], cellFaces[Y]);
        CalculateInterpolationFactors(interpFactors[Z], cellCenters[Z], cellFaces[Z]);

        CalculateExtrapolationFactors(extrapFactors, cellLengths, X);
        CalculateExtrapolationFactors(extrapFactors, cellLengths, Y);
        CalculateExtrapolationFactors(extrapFactors, cellLengths, Z);
    };


/*-------------------------------------------------------------------------------------*\
                                    FVCoefficients
\*-------------------------------------------------------------------------------------*/

using C = TransportCoefficients::ENUMDATA;

FVCoefficients::FVCoefficients(const indexVector3 &dims) :
    auu({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
    avv({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
    aww({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),

    aup({C::p, C::e, C::w}, dims(0)),
    avp({C::p, C::n, C::s}, dims(1)),
    awp({C::p, C::t, C::b}, dims(2)),

    acu({C::p, C::e, C::w}, dims(0)),
    acv({C::p, C::n, C::s}, dims(1)),
    acw({C::p, C::t, C::b}, dims(2)),
    
    acp({C::p, C::n, C::e, C::s, C::w, C::t, C::b, C::nn, C::ee, C::ss, C::ww, C::tt, C::bb}, dims),

    bu(dims(0), dims(1), dims(2)),
    bv(dims(0), dims(1), dims(2)),
    bw(dims(0), dims(1), dims(2)),
    bc(dims(0), dims(1), dims(2))
    {};


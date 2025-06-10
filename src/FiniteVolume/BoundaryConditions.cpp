#include "FiniteVolumeFunctions.h"
#include "../IO/InputProcessing.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
#include "../IO/VTKReader.h"

#include <cmath>
#include <stdexcept>

namespace CAMIRA
{

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




BoundaryConditionData SetBoundaryConditionData( const InputData &inputData,
                                                const Mesh &mesh )
{

    BoundaryConditionData bcData;

    // Set boundary condition data for use in solver
    ForAllFieldData( [&] (intType f) {

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

            const auto &bcConfig = inputData.boundaryConditions[f][bp];

            bcData.fields[f][bp].type = bcConfig.type;

            if ( bcConfig.hasUniformValue ) 
                bcData.fields[f][bp].value = SetBoundaryProfileConstant( bcConfig.uniformValue, mesh, bp );

            if ( bcConfig.hasProfile1D ) 
                bcData.fields[f][bp].value = SetBoundaryProfile1D( bcConfig.profile1D, mesh, bp );

        } );

    } );


    // Pressure field will be floating if none of the pathces are fixed
    bcData.pressureFieldIsFloating = true;
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        if ( bcData.fields.P[bp].type == BoundaryConditions::fixed ) {
            bcData.pressureFieldIsFloating = false;
        }
    } );

    return bcData;
}



}   // end namespace CAMIRA
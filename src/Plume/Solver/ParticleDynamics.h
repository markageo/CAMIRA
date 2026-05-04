#ifndef CAMIRA_PLUME_LAGRANGIAN
#define CAMIRA_PLUME_LAGRANGIAN

#include "Core/Types.h"
#include "Core/FVLookups.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "Plume/Particles/Particles.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/ConfigEnums.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{


namespace 
{


inline floatType signum( floatType f) {
    if ( f > 0.0f ) return  1.0f;
    if ( f < 0.0f ) return -1.0f;
    return 0.0f;
}



class StochasticVariableEta {

    public:
        
        StochasticVariableEta() : rd(), gen(rd()), d(0.0f, 1.0f) {};

        floatType operator()() 
        { return d(gen); }

    private:
        std::random_device rd;
        std::mt19937 gen;
        std::normal_distribution<floatType> d;
};



inline fVector3 GetParticleStep( const fVector3 &velocity,
                                 const floatType turbulentDiffusivity,
                                 const fVector3 &gradTurbulentDiffusivity,
                                 const InputData &inputData,
                                 StochasticVariableEta &eta )
{

    const fVector3 advection = inputData.timeStepSize * velocity;

    fVector3 etaVec = { eta(), eta(), eta() };

    const fVector3 molecularDiffusion = sqrt( 2.0f * inputData.diffusionCoeff * inputData.timeStepSize ) * etaVec;

    const fVector3 turbulentDiffusion = gradTurbulentDiffusivity * inputData.timeStepSize
                                      + sqrt( 2.0f * turbulentDiffusivity * inputData.timeStepSize ) * etaVec;

    return advection + molecularDiffusion + turbulentDiffusion;
}



fVector3 RecursiveReflection( const Tree &tree,
                              const fVector3 oldPosition,
                              const fVector3 newPosition )
{
    const fVector3 currentDelta = newPosition - oldPosition;
    const auto intersectionData = SegmentIntersectionAndNormal( tree, oldPosition, currentDelta );
    const fVector3 intersectionPoint = intersectionData.point;
    const fVector3 normalVector = intersectionData.normal; 

    fVector3 reflectedPosition = newPosition - 2.0f * ( normalVector.dot( newPosition - intersectionPoint )  ) * normalVector;

    // Shift the intersection point along delta to prevent an unintended intersection from initial intersection point
    const floatType shiftScale = 1e-8;
    const fVector3 shiftedIntersectionPoint = intersectionPoint + shiftScale * (reflectedPosition - intersectionPoint);

    if ( SegmentIntersects( tree, shiftedIntersectionPoint, reflectedPosition ) ) {
        reflectedPosition = RecursiveReflection( tree, shiftedIntersectionPoint, reflectedPosition );
    }

    return reflectedPosition;
}



void StepParticle( Particles &particles,
                   const intType idx,
                   const fVector3 &delta,
                   const Mesh &mesh, 
                   const Tree &tree,
                   const EnumVector<BoundaryPatches, InputData::BoundaryConditionInputData> &boundaryConditions )
{
    using enum Axis::ENUMDATA;

    const fVector3 oldPosition = { particles.x[idx], particles.y[idx], particles.z[idx] };
    fVector3 newPosition = oldPosition + delta;

    // Domain boundary intersection
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        floatType hiBounds = mesh.cellFaces[axis]( mesh.nFacesNormal[axis][axis]-1 ),
                  loBounds = mesh.cellFaces[axis]( 0 );

        BoundaryPatches::ENUMDATA boundaryPatch;
        bool exitedBounds = false;
        
        if        ( newPosition(axis) >= hiBounds ) {
            boundaryPatch = LUT::PositivePatch[axis];
            exitedBounds = true;
        } else if ( newPosition(axis) <= loBounds ) {
            boundaryPatch = LUT::NegativePatch[axis];
            exitedBounds = true;
        }


        if ( exitedBounds ) {

            switch ( boundaryConditions[ boundaryPatch ].type ) {

                case BoundaryConditions::outflow:
                    particles.active[idx] = false;
                    break;

                case BoundaryConditions::reflection:    
                {
                    const floatType wallPosition = ( newPosition(axis) >= hiBounds ) ? hiBounds 
                                                                                     : loBounds;
                    newPosition(axis) = 2.0f * wallPosition - newPosition(axis);
                    break;
                }

                case BoundaryConditions::periodic:
                    newPosition(axis) = signum( delta(axis) ) * ( loBounds - hiBounds ) + oldPosition(axis) + delta(axis); 
                    break;
            }

        }

    } );


    // Solid geometry intersection and reflection 
    if ( !tree.empty() ) {
        if ( SegmentIntersects( tree, oldPosition, newPosition ) ) {
            newPosition = RecursiveReflection( tree, oldPosition, newPosition );
        }
    }
    
    particles.x[idx] = newPosition(0);
    particles.y[idx] = newPosition(1);
    particles.z[idx] = newPosition(2);
}



}   // end anonymous namespace



inline void UpdateParticles( Particles &particles,
                             const Mesh &mesh, 
                             const EnumVector<Axis, Tensor3D> &velocityField,
                             const Tensor3D &nuTurbField,
                             const Tree &tree,
                             const InputData &inputData )
{

    #pragma omp parallel 
    {

    StochasticVariableEta eta;

    #pragma omp for
    for ( intType idx = 0; idx != particles.N; idx++ ) {

        UpdateParticlePositionIndexLinearSearch( particles, idx, mesh );

        fVector3 localVelocity;
        EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
            localVelocity[axis] = GetFieldQuantityTrilinearInterp( particles, idx, mesh, velocityField[axis] );
        } );

        const floatType turbulentDiffusivity = GetFieldQuantityTrilinearInterp( particles, idx, mesh, nuTurbField ) / inputData.turbulentSchmidtNumber;

        const fVector3 gradTurbulentDiffisivity = GetFieldQuantityGradient( particles, idx, mesh, nuTurbField ) / inputData.turbulentSchmidtNumber;

        const fVector3 delta = GetParticleStep( localVelocity, turbulentDiffusivity, gradTurbulentDiffisivity, inputData, eta );

        StepParticle( particles, idx, delta, mesh, tree, inputData.boundaryConditions );

    }

    }   // end omp parallel region


    // Remove any inactive particles
    particles.RemoveInactiveParticles();

}


inline void SplitParticles( Particles &particles )
{

    const intType initialNumberOfParticles = particles.N;

    for ( intType idx = 0; idx != initialNumberOfParticles; idx++ ) {

        // Half the mass of the particle
        particles.mass[idx] /= 2.0f;

        // Copy the particle to the end
        particles.AddIdenticalParticleBack( idx );

    }

}




}   // end namespace PLUME

}   // end namespace CAMIRA    


#endif // CAMIRA_PLUME_LAGRANGIAN
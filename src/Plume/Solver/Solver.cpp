#include "Solver.h"

#include "Core/Types.h"
#include "Core/Macros.h"
#include "Core/IO/VTKReader.h"
#include "Core/Mesh/Mesh.h"
#include "Plume/Particle/Particle.h"
#include "Plume/Sources/Sources.h"
#include "Plume/Solver/ParticleDynamics.h"
#include "Plume/Concentration/Concentration.h"

#include <iostream>
#include <memory>
#include <cmath>

namespace CAMIRA
{
    
using namespace CORE;

namespace PLUME
{

namespace
{

    // Calculates the total number of particles needed throughout the simulation, accounting for different souces so that memory can be 
    // allocated for them at the start
    // TODO: should account for particle splitting when this is added
    intType CalculateNumberOfParticlesNeeded( const InputData &inputData )
    {
        intType nParticles = 0;

        // Due to instantanious releases
        for ( const auto &release : inputData.instantaneousReleasePoints ) {
           intType particlesForSource = static_cast<intType>( std::ceil( release.totalMass * inputData.initialParticlesPerUnitMass ) );
           nParticles += particlesForSource;
        }

        // Due to continuous releases
        for ( const auto &release : inputData.continuousReleasePoints ) {
            floatType totalTime = static_cast<floatType>( inputData.numberOfTimeSteps ) * inputData.timeStepSize;
            intType particlesForSource = static_cast<intType>(  std::ceil(   release.massFlowRate 
                                                                           * totalTime 
                                                                           * static_cast<floatType>(inputData.initialParticlesPerUnitMass) ) );
            nParticles += particlesForSource;
        }

        return nParticles;
    }

}


void SolvePlume( const InputData &inputData )
{

    // Read and store windfield
    VTK::FieldFileData fieldFileData = VTK::ReadVTKFields( inputData.velocityFieldFilename );
    Mesh mesh( fieldFileData.cellFaces );
    EnumVector<Axis, Tensor3D> velocityField = fieldFileData.vertexFields.U;

    // Read and store geometry


    // Place to store the concentration field
    Tensor3D concentrationField( mesh.nCells(0), mesh.nCells(1), mesh.nCells(2) );
    concentrationField.setZero();

    // Allocate particles
    intType particlesNeeded = CalculateNumberOfParticlesNeeded( inputData );
    std::vector< Particle > particles;
    particles.reserve( particlesNeeded );

    // Add instantaneous release particles
    AddInstantaneousReleasePointParticles( particles, mesh, inputData );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Step through time
    floatType currentTime = 0;
    for ( intType timeStep = 1; timeStep != inputData.numberOfTimeSteps; timeStep++ ) {
        currentTime += inputData.timeStepSize;

        // Add continuous release particles
        AddContinuousReleasePointParticles( particles, mesh, inputData );

        // Update particle positions
        UpdateParticles( particles, mesh, velocityField, inputData );

        // Output
        if ( writeFields && (timeStep % inputData.fieldWriteInterval) == 0 ) {
            for ( auto &particle : particles ) {
                UpdateParticlePositionIndexBinarySearch( particle, mesh );
            }
            UpdateConcentrationField( concentrationField, particles, mesh );
        } 
        


    }

}


}   // end namespace PLUME

}   // end namespace CAMIRA
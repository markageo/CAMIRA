#include "Solver.h"

#include "Core/Types.h"
#include "Core/Macros.h"
#include "Core/IO/VTKReader.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "Plume/Particle/Particle.h"
#include "Plume/Sources/Sources.h"
#include "Plume/Solver/ParticleDynamics.h"
#include "Plume/Concentration/Concentration.h"
#include "Plume/Solver/Logging.h"

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
    // This is an estimate as it doesn't account for particles potentially leaving the domain
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
    EnumVector<Axis, Tensor3D> &velocityField = fieldFileData.vertexFields.U;

    // Read and store geometry
    Polyhedron P = MakeGeometry( inputData.stlGeometryFilename );
    Tree tree = MakeAABBTree( P );

    // Place to store the concentration field
    Tensor3D concentrationField( mesh.nCells(0), mesh.nCells(1), mesh.nCells(2) );
    concentrationField.setZero();

    // Allocate particles
    intType particlesNeeded = CalculateNumberOfParticlesNeeded( inputData );
    std::vector< Particle > particles;
    particles.reserve( particlesNeeded );
    bool splitParticles = ( inputData.particleSplitTimeStepInterval > 0 );
    intType numberOfParticleSplits = 0;

    // Logging objects
    ConcentrationFieldWriter concentrationFieldWriter( concentrationField, mesh, inputData );

    // Add instantaneous release particles
    AddInstantaneousReleasePointParticles( particles, mesh, inputData );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Write the initial condition
    UpdateConcentrationField( concentrationField, particles, mesh );
    concentrationFieldWriter.WriteData( 0.0f, 0 );

    // Step through time
    std::cout << "Starting time loop " << std::endl;
    floatType currentTime = 0;
    for ( intType timeStep = 1; timeStep <= inputData.numberOfTimeSteps; timeStep++ ) {
        currentTime += inputData.timeStepSize;

        std::cout << "Timestep: " << timeStep << "\n"
                  << "Number of particles: " << particles.size() << "\n"
                  << std::flush;

        // Split particles 
        bool splitParticlesThisIteration = splitParticles
                                        && ( numberOfParticleSplits < inputData.maxNumberOfParticleSplits )
                                        && ( timeStep % inputData.particleSplitTimeStepInterval == 0 );
        if ( splitParticlesThisIteration ) {
            numberOfParticleSplits++;
            SplitParticles( particles );
        }

        // Add continuous release particles
        AddContinuousReleasePointParticles( particles, mesh, inputData );
        
        // Update particle positions
        UpdateParticles( particles, mesh, velocityField, tree, inputData );

        // Output
        if ( writeFields && (timeStep % inputData.fieldWriteInterval) == 0 ) {
            for ( auto &particle : particles ) {
                UpdateParticlePositionIndexLinearSearch( particle, mesh );
            }
            UpdateConcentrationField( concentrationField, particles, mesh );
            concentrationFieldWriter.WriteData( currentTime, timeStep );
            std::cout << "(Concentration field written to file)" << "\n";
        } 
        
        std::cout << std::endl;
    }

}


}   // end namespace PLUME

}   // end namespace CAMIRA
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

    TIC("Pre processing")
    std::cout << "Setting up solver data and pre-processing ... \n\n" << std::flush;

    // Read and store windfield
    VTK::FieldFileData fieldFileData = VTK::ReadVTKFields( inputData.velocityFieldFilename );
    Mesh mesh( fieldFileData.cellFaces );
    EnumVector<Axis, Tensor3D> &velocityField = fieldFileData.vertexFields.U;
    Tensor3D &nuTurbField = fieldFileData.vertexNuTurb;

    // Read and store geometry
    Tree tree;
    if( inputData.hasSolidGeometry ) {
        Polyhedron P = MakeGeometry( inputData.stlGeometryFilename );
        tree = MakeAABBTree( P );
    }

    // Place to store the current concentration field
    Tensor3D concentrationField( mesh.nCells(0), mesh.nCells(1), mesh.nCells(2) );
    SetTensorZeroParallel( concentrationField );

    // Time averaged concentration fields
    std::vector<Tensor3D> timeAveragedConcentrationFields( inputData.timeAveragedConcentrationFieldData.size() );  
    for ( auto &field : timeAveragedConcentrationFields ) {
        field = Tensor3D( mesh.nCells(0), mesh.nCells(1), mesh.nCells(2) );
        SetTensorZeroParallel( field );
    }

    // Allocate particles
    intType particlesNeeded = CalculateNumberOfParticlesNeeded( inputData );
    std::vector< Particle > particles;
    particles.reserve( particlesNeeded );
    bool splitParticles = ( inputData.particleSplitTimeStepInterval > 0 );
    intType numberOfParticleSplits = 0;

    // Logging objects
    ConcentrationFieldWriter concentrationFieldWriter( mesh, inputData );

    // Add instantaneous release particles
    AddInstantaneousReleasePointParticles( particles, mesh, inputData );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Write the initial condition
    UpdateConcentrationField( concentrationField, particles, mesh );
    concentrationFieldWriter.WriteInstantaneousField( concentrationField, 0.0f, 0 );

    std::cout << "Pre-processing complete. Starting solve.\n\n" << std::flush;
    TOC()

    // Step through time
    TIC("Solver Loop")
    std::cout << "Starting time loop " << std::endl;
    floatType currentTime = 0;
    for ( intType timeStep = 1; timeStep <= inputData.numberOfTimeSteps; timeStep++ ) {
        currentTime += inputData.timeStepSize;

        std::cout << "Timestep: " << timeStep << "\n"
                  << "Number of particles: " << particles.size() << "\n"
                  << std::flush;

        // Split particles 
        TIC("Particle Splitting")
        bool splitParticlesThisIteration = splitParticles
                                        && ( numberOfParticleSplits < inputData.maxNumberOfParticleSplits )
                                        && ( timeStep % inputData.particleSplitTimeStepInterval == 0 );
        if ( splitParticlesThisIteration ) {
            numberOfParticleSplits++;
            SplitParticles( particles );
        }
        TOC()

        // Add continuous release particles
        TIC("Add Continuous Release")
        AddContinuousReleasePointParticles( particles, mesh, inputData );
        TOC()
        
        // Update particle positions
        TIC("Particle Update")
        UpdateParticles( particles, mesh, velocityField, nuTurbField, tree, inputData );
        TOC()

        // Concentration field, only update this iteration if it will be used
        bool outputThisIteration = ( writeFields && (timeStep % inputData.fieldWriteInterval) == 0 )
                                || ( timeStep == inputData.numberOfTimeSteps );
        bool updateConcentrationFieldThisIteration = outputThisIteration;
        if ( !updateConcentrationFieldThisIteration ) {
            for ( const auto &data : inputData.timeAveragedConcentrationFieldData ) {
                if ( timeStep > data.endTimeStep )
                    continue;

                if ( (timeStep - data.startTimeStep) % data.timeStepInterval == 0 )
                    updateConcentrationFieldThisIteration = true;
            }
        }

        if ( updateConcentrationFieldThisIteration ) {
            #pragma omp parallel for
            for ( auto &particle : particles ) {
                UpdateParticlePositionIndexLinearSearch( particle, mesh );
            }
            UpdateConcentrationField( concentrationField, particles, mesh );
        }

        // Instantanious concentration field output
        TIC("Field Output")
        if ( outputThisIteration ) {
            concentrationFieldWriter.WriteInstantaneousField( concentrationField, currentTime, timeStep );
            std::cout << "(Concentration field written to file)" << "\n";
        } 
        TOC()

        // Time averaged concentrations
        TIC("Time Averaged Concentration Fields")
        for ( size_t i = 0; i != timeAveragedConcentrationFields.size(); i++ ) {

            const auto data = inputData.timeAveragedConcentrationFieldData[i];

            if ( (timeStep - data.startTimeStep) % data.timeStepInterval != 0 )
                continue;

            timeAveragedConcentrationFields[i] += concentrationField * concentrationField.constant( static_cast<floatType>(data.timeStepInterval) * inputData.timeStepSize );

        }
        TOC()

        std::cout << std::endl;
    }

    // Finish average calculation and write to file
    for ( size_t i = 0; i != timeAveragedConcentrationFields.size(); i++ ) {
        
        const auto data = inputData.timeAveragedConcentrationFieldData[i];

        timeAveragedConcentrationFields[i] /= timeAveragedConcentrationFields[i].constant( static_cast<floatType>( data.endTimeStep - data.startTimeStep + 1 ) * inputData.timeStepSize );

        concentrationFieldWriter.WriteTimeAveragedField( timeAveragedConcentrationFields[i], data.filename, data.startTimeStep, data.endTimeStep, inputData.timeStepSize );

    }
    TOC()

}


}   // end namespace PLUME

}   // end namespace CAMIRA
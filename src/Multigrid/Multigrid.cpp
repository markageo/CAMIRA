#include "Multigrid.h"

#include "../FiniteVolume/Mesh.h"

namespace CFD
{

namespace
{



} // end anonymous namespace



std::vector<GridLevelData> CreateMGLevels( const InputData &inputData )
{

    const InputData::MultigridSettings &mgSettings = inputData.multigridSettings;

    std::vector<GridLevelData> mgLevels;

    mgLevels.emplace_back();
    mgLevels[0].mesh = CreateMesh( inputData );

    for ( intType level = 1; level != mgSettings.maxCoarseLevels; level++ ) {

        if ( !MeshCanBeCoarsened( mgLevels[level-1].mesh ) ) {
            break;
        }

        mgLevels.emplace_back();
        mgLevels[level].mesh = CoarsenMesh( mgLevels[level-1].mesh );
        
    }


    return mgLevels;

}



}   // end namespace CFD
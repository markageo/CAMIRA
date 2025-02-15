#ifndef FV_TOOLS
#define FV_TOOLS

#include "Types.h"
#include "../FiniteVolume/Mesh.h"


namespace CFD 
{

// Number of ghost cells in solution field
constexpr intType nGhost = 2;

// Finite Volume Tools
namespace FVT
{

    // -------------------------------------------------- Ghost Cell Indexing ------------------------------------------------- //

    // Eigen::Tensor ghost cell indexing - array indexing
    template<std::size_t dim> __attribute__((always_inline)) 
    inline Eigen::array<Eigen::Index, dim> G( Eigen::array<Eigen::Index, dim> idx )
    {
        for (auto it = idx.begin(); it != idx.end(); it++) {
            *it += nGhost;                               
        }
        return idx;
    }


    // Eigen::Tensor ghost cell indexing - list indexing
    template<class ...Args> 
    requires( sizeof...(Args) > 1 ) __attribute__((always_inline)) 
    inline Eigen::array< Eigen::Index, sizeof...(Args) > G( Args... args )
    {  
        Eigen::array< Eigen::Index, sizeof...(Args) > idx( {args...} );
        return G( idx ); 
    }

    // Special case for 1D, should just return a number
    __attribute__((always_inline)) 
    inline Eigen::Index G( Eigen::Index idx )
    {
        return idx + nGhost;
    }





    // --------------------------------------------------- Cell Face Indexing -------------------------------------------------- //

    // Returns start and end indexing for looping over internal faces of a mesh for a given normal direction.
    // Intended to be used as returning a structured binding.
    inline std::tuple<iArray3, iArray3> FaceInternalIndices( const Mesh &mesh, 
                                                             const Axis::ENUMDATA axis )
    {
        iArray3 startIndex = {0, 0, 0};
        iArray3 nFaces = mesh.nFacesNormal[axis];
        startIndex(axis) += 1;
        nFaces(axis) -= 1;
        return { startIndex, nFaces };
    }





    // --------------------------------------------------- Remove Ghost Cells -------------------------------------------------- //


    // Return a copy with ghost cells removed. Has to make a copy anyway.
    template< typename arrayType >
    arrayType RemoveGhostCells( const arrayType &array, 
                                const intType nGhostCells)
    {
        static_assert(std::is_same<arrayType, CFD::Tensor1D   >::value ||
                      std::is_same<arrayType, CFD::Tensor2D   >::value ||
                      std::is_same<arrayType, CFD::Tensor3D   >::value,
                      "Template parameter must be a Tensor.");

        Eigen::array< Eigen::Index, arrayType::NumDimensions > offsets, extents;
        for ( size_t i = 0; i != array.NumDimensions; i++ ) {
            offsets[i] = nGhostCells;
            extents[i] = array.dimension(i) - 2*nGhostCells;
        }

        return arrayType( array ).slice(offsets, extents);
    }


}   // end namespace FVT

}   // end namespace CFD

#endif // FV_TOOLS
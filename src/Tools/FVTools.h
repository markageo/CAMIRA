#ifndef FV_TOOLS
#define FV_TOOLS

#include "../Types.h"


namespace CFD 
{

// Number of ghost cells in solution field
constexpr intType nGhost = 2;

// Finite Volume Tools
namespace FVT
{

    // -------------------------------------------------- Ghost Cell Indexing ------------------------------------------------- //

    // Eigen::Tensor ghost cell indexing - array indexing
    template<std::size_t dim>
    inline Eigen::array<Eigen::Index, dim> G( Eigen::array<Eigen::Index, dim> idx )
    {
        for (auto it = idx.begin(); it != idx.end(); it++) {
            *it += nGhost;                               
        }
        return idx;
    }


    // Eigen::Tensor ghost cell indexing - list indexing
    template<class ...Args> 
    requires( sizeof...(Args) > 1 )
    inline Eigen::array< Eigen::Index, sizeof...(Args) > G( Args... args )
    {  
        Eigen::array< Eigen::Index, sizeof...(Args) > idx( {args...} );
        return G( idx ); 
    }

    // Special case for 1D, should just return a number
    inline Eigen::Index G( Eigen::Index idx )
    {
        return idx + nGhost;
    }



    // --------------------------------------------------- Cell Face Indexing -------------------------------------------------- //




}   // end namespace FVU

}   // end namespace CFD

#endif // FV_TOOLS
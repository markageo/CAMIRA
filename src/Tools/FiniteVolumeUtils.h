#ifndef CFD_UTILS
#define CFD_UTILS

#include "../Types.h"


namespace CFD 
{


// Finite Volume Utilities
namespace FVU
{

    // -------------------------------------------------- Ghost Cell Indexing ------------------------------------------------- //

    // Number of ghost cells in solution field
    constexpr intType nGhost = 2;


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





// Lookup Tables
namespace LUT
{

    // ----------------------------------------------------- Enum Lookups ----------------------------------------------------- //

    // Get index offset from TransportCoeffienct
    constexpr std::array<intType, TransportCoefficients::count> CoeffIndex = { 2,     // tt
                                                                            2,     // nn    
                                                                            2,     // ee
                                                                            1,     // t
                                                                            1,     // n
                                                                            1,     // e
                                                                            0,     // p
                                                                            -1,    // w
                                                                            -1,    // s
                                                                            -1,    // b
                                                                            -2,    // ww
                                                                            -2,    // ss
                                                                            -2};   // bb


    // Get BoundaryPatches from Axis
    constexpr std::array<BoundaryPatches::ENUMDATA, 3> PositivePatch = { BoundaryPatches::xPositive,
                                                                        BoundaryPatches::yPositive,
                                                                        BoundaryPatches::zPositive};

    constexpr std::array<BoundaryPatches::ENUMDATA, 3> NegativePatch{ BoundaryPatches::xNegative,
                                                                    BoundaryPatches::yNegative,
                                                                    BoundaryPatches::zNegative};

    // Get TransportCoefficient from Axis
    constexpr std::array<TransportCoefficients::ENUMDATA, 3> HiCoeff{ TransportCoefficients::e,
                                                                    TransportCoefficients::n,
                                                                    TransportCoefficients::t};
                                                                    
    constexpr std::array<TransportCoefficients::ENUMDATA, 3> HiHiCoeff{ TransportCoefficients::ee,
                                                                        TransportCoefficients::nn,
                                                                        TransportCoefficients::tt};


    constexpr std::array<TransportCoefficients::ENUMDATA, 3> LoCoeff{ TransportCoefficients::w,
                                                                    TransportCoefficients::s,
                                                                    TransportCoefficients::b};

    constexpr std::array<TransportCoefficients::ENUMDATA, 3> LoLoCoeff{ TransportCoefficients::ww,
                                                                        TransportCoefficients::ss,
                                                                        TransportCoefficients::bb};                                                                 


    // Get Axis from BoundaryPatches
    constexpr std::array<Axis::ENUMDATA, 6> BoundaryPatchAxis{ Axis::X,    // xPositive
                                                            Axis::X,    // xNegative
                                                            Axis::Y,    // yPositive
                                                            Axis::Y,    // yNegative
                                                            Axis::Z,    // zPositive
                                                            Axis::Z};   // zNegative


    // Axis orthogonal to a given one
    constexpr std::array<Axis::ENUMDATA, 3> LoOrthogonalAxis{ Axis::Y,
                                                            Axis::X,
                                                            Axis::X };

    constexpr std::array<Axis::ENUMDATA, 3> HiOrthogonalAxis{ Axis::Z,
                                                            Axis::Z,
                                                            Axis::Y };
 
}   // end namespace LUT


}   // end namespace CFD

#endif // CFD_UTILS
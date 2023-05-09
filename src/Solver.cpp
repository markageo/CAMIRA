#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"



namespace
{

using namespace CFD;

// Update the residual of the fields. Calculated as the L1 norm between current and previous iterations
void UpdateResiduals( EnumVector<Fields, floatType> &residuals,
                     const ArrayAllocator<Fields, array3D> &fields,
                     const ArrayAllocator<Fields, array3D> &fieldsOld,
                     const EnumVector<Fields, floatType> &residualsInitial = EnumVector<Fields, floatType>({1.0, 1.0, 1.0, 1.0}))
{
    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<Fields::ENUMDATA>(f);

        auto fieldDiff = fields[field] - fieldsOld[field];  // auto delays calculation
        Eigen::Tensor<floatType, 0> temp = fieldDiff.abs().mean();
        residuals[field] = temp(0) / residualsInitial[field];
    }
}


// Check it residual tolerence is met
bool MetResidualTolerence( const EnumVector<Fields, floatType> &residuals,
                           const EnumVector<Fields, floatType> &residualsTarget)
{
    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<Fields::ENUMDATA>(f);

        if ( residuals[field] > residualsTarget[field] )
            return false;

    }
    return true;
}

}   // end anonymous namespace




// Performs a single local update of block coupled equations
class BlockSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        BlockSolver(ArrayAllocator<Fields, array3D> &fields, const FVCoefficients &fvCoeffs) :
            m_fields( fields ), m_fvCoeffs( fvCoeffs ) {};


        // Core function which updates the local coupled system. Templated by staggering direction.
        template<TC Ustag, TC Vstag, TC Wstag >
        void UpdateBlock(const intType i, const intType j, const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;

            // Staggering must be valid
            static_assert( (Ustag == e) || (Ustag == w) || (Ustag == p), "Invalid U momentum staggering" );
            static_assert( (Vstag == n) || (Vstag == s) || (Vstag == p), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Index offset for the coupled cell
            static constexpr intType iUoffset = UoffsetCoupled(Ustag), 
                                     jVoffset = VoffsetCoupled(Vstag), 
                                     kWoffset = WoffsetCoupled(Wstag);

            // Index offset for the uncoupled pressure coordinates
            static constexpr intType iUleft = OffsetLeft(iUoffset), iUright = OffsetRight(iUoffset),
                                     jVleft = OffsetLeft(jVoffset), jVright = OffsetRight(jVoffset),
                                     kWleft = OffsetLeft(kWoffset), kWright = OffsetRight(kWoffset);
                                    

            // Compass coordinate of the pressure that is coupled into the continuity equation for each equation
            static constexpr TC cUcoupled = UcompassCoupled(Ustag),
                                cVcoupled = VcompassCoupled(Vstag), 
                                cWcoupled = WcompassCoupled(Wstag);

            // Compass coordinate of the uncoupled pressure coordinates for each equation
            static constexpr TC cUleft = UcompassLeft(cUcoupled), cUright = UcompassRight(cUcoupled),
                                cVleft = VcompassLeft(cVcoupled), cVright = VcompassRight(cVcoupled),
                                cWleft = WcompassLeft(cWcoupled), cWright = WcompassRight(cWcoupled);

            // Forces compile time evaluation
            static_assert( (iUoffset == 1) || (iUoffset == -1) || (iUoffset == 0) );        
            static_assert( (jVoffset == 1) || (jVoffset == -1) || (jVoffset == 0) );  
            static_assert( (kWoffset == 1) || (kWoffset == -1) || (kWoffset == 0) );  

            static_assert( (iUleft == -1) || (iUleft == 0) ); 
            static_assert( (jVleft == -1) || (jVleft == 0) ); 
            static_assert( (kWleft == -1) || (kWleft == 0) );  

            static_assert( (iUright == 1) || (iUright == 0) ); 
            static_assert( (jVright == 1) || (jVright == 0) ); 
            static_assert( (kWright == 1) || (kWright == 0) ); 

            static_assert( (cUcoupled == w) || (cUcoupled == p) || (cUcoupled == e) );  
            static_assert( (cUleft == w) || (cUleft == p) );  
            static_assert( (cUright == e) || (cUright == p) ); 

            static_assert( (cVcoupled == s) || (cVcoupled == p) || (cVcoupled == n) );  
            static_assert( (cVleft == s) || (cVleft == p) );  
            static_assert( (cVright == n) || (cVright == p) );  

            static_assert( (cWcoupled == b) || (cWcoupled == p) || (cWcoupled == t) );  
            static_assert( (cWleft == b) || (cWleft == p) );  
            static_assert( (cWright == t) || (cWright == p) );     


            // For indexing the staggered cells
            intType iU(i + iUoffset), jU(j           ), kU(k           ); // U momentum
            intType iV(i           ), jV(j + jVoffset), kV(k           ); // V momentum
            intType iW(i           ), jW(j           ), kW(k + kWoffset); // W momentum


            // Precompute momentum RHS divided by AP coefficients
            // U momentum
            floatType bU = ( m_fvCoeffs.Umom.B(iU, jU, kU)

                           - m_fvCoeffs.Umom.AU[n](iU, jU, kU) * m_fields[U]( G(iU  , jU+1, kU  ) )
                           - m_fvCoeffs.Umom.AU[e](iU, jU, kU) * m_fields[U]( G(iU+1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[s](iU, jU, kU) * m_fields[U]( G(iU  , jU-1, kU  ) )
                           - m_fvCoeffs.Umom.AU[w](iU, jU, kU) * m_fields[U]( G(iU-1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[t](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU+1) )
                           - m_fvCoeffs.Umom.AU[b](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU-1) )

                           - m_fvCoeffs.Umom.AP[cUleft ](iU) * m_fields[P]( G(iU+iUleft , jU  , kU  ) )
                           - m_fvCoeffs.Umom.AP[cUright](iU) * m_fields[P]( G(iU+iUright, jU  , kU  ) )

                           ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // V momentum 
            floatType bV = ( m_fvCoeffs.Vmom.B(iV, jV, kV)

                           - m_fvCoeffs.Vmom.AV[n](iV, jV, kV) * m_fields[V]( G(iV  , jV+1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[e](iV, jV, kV) * m_fields[V]( G(iV+1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[s](iV, jV, kV) * m_fields[V]( G(iV  , jV-1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[w](iV, jV, kV) * m_fields[V]( G(iV-1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[t](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV+1) )
                           - m_fvCoeffs.Vmom.AV[b](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV-1) )

                           - m_fvCoeffs.Vmom.AP[cVleft ](jV) * m_fields[P]( G(iV  , jV+jVleft , kV  ) )
                           - m_fvCoeffs.Vmom.AP[cVright](jV) * m_fields[P]( G(iV  , jV+jVright, kV  ) )

                           ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);


            // W momentum
            floatType bW = ( m_fvCoeffs.Wmom.B(iW, jW, kW)

                           - m_fvCoeffs.Wmom.AW[n](iW, jW, kW) * m_fields[W]( G(iW  , jW+1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[e](iW, jW, kW) * m_fields[W]( G(iW+1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[s](iW, jW, kW) * m_fields[W]( G(iW  , jW-1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[w](iW, jW, kW) * m_fields[W]( G(iW-1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[t](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW+1) )
                           - m_fvCoeffs.Wmom.AW[b](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW-1) )

                           - m_fvCoeffs.Wmom.AP[cWleft ](kW) * m_fields[P]( G(iW  , jW  , kW+kWleft ) )
                           - m_fvCoeffs.Wmom.AP[cWright](kW) * m_fields[P]( G(iW  , jW  , kW+kWright) )

                           ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);



            // Continuity for pressure
            floatType bP = (m_fvCoeffs.Cont.B(i, j, k)
                            
                          //- m_fvCoeffs.Cont.AU[e](i) * m_fields[U]( G(i+1, j, k) )
                          - m_fvCoeffs.Cont.AU[p](i) * m_fields[U]( G(i  , j, k) )
                          - m_fvCoeffs.Cont.AU[w](i) * m_fields[U]( G(i-1, j, k) )

                          - m_fvCoeffs.Cont.AU[e](i) * m_fields[U]( G(i+1, j, k) )
                          //- m_fvCoeffs.Cont.AU[p](i) * m_fields[U]( G(i  , j, k) )
                          - m_fvCoeffs.Cont.AU[w](i) * m_fields[U]( G(i-1, j, k) )

                          - m_fvCoeffs.Cont.AU[e](i) * m_fields[U]( G(i+1, j, k) )
                          - m_fvCoeffs.Cont.AU[p](i) * m_fields[U]( G(i  , j, k) )
                          //- m_fvCoeffs.Cont.AU[w](i) * m_fields[U]( G(i-1, j, k) )




                    
                          //- m_fvCoeffs.Cont.AV[n](j) * m_fields[V]( G(i, j+1, k) )
                          - m_fvCoeffs.Cont.AV[p](j) * m_fields[V]( G(i, j  , k) )
                          - m_fvCoeffs.Cont.AV[s](j) * m_fields[V]( G(i, j-1, k) )
                          
                          //- m_fvCoeffs.Cont.AW[t](k) * m_fields[W]( G(i, j, k+1) )
                          - m_fvCoeffs.Cont.AW[p](k) * m_fields[W]( G(i, j, k  ) )
                          - m_fvCoeffs.Cont.AW[b](k) * m_fields[W]( G(i, j, k-1) )
                          
                          - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields[P]( G(i  , j+1, k  ) )
                          - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields[P]( G(i+1, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields[P]( G(i  , j-1, k  ) )
                          - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields[P]( G(i-1, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields[P]( G(i  , j  , k+1) )
                          - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields[P]( G(i  , j  , k-1) )

                          - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields[P]( G(i  , j+2, k  ) )
                          - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields[P]( G(i+2, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields[P]( G(i  , j-2, k  ) )
                          - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields[P]( G(i-2, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields[P]( G(i  , j  , k+2) )
                          - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields[P]( G(i  , j  , k-2) )

                          ) / m_fvCoeffs.Cont.AP[p](i, j, k);

            // This only needs to be updated at linearisation
            floatType K = m_fvCoeffs.Cont.AP[p](i, j, k)
                        - m_fvCoeffs.Cont.AU[e](i) * m_fvCoeffs.Umom.AP[w](iU) / m_fvCoeffs.Umom.AU[p](iU, jU, kU)
                        - m_fvCoeffs.Cont.AV[n](j) * m_fvCoeffs.Vmom.AP[s](jV) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV)
                        - m_fvCoeffs.Cont.AW[t](k) * m_fvCoeffs.Wmom.AP[b](kW) / m_fvCoeffs.Umom.AW[p](iW, jW, kW);
            K = 1.0f / K;


            // Update P from continuity
            m_fields[P]( G(i, j, k) ) = ( bP
                                        - m_fvCoeffs.Cont.AU[e](i) * bU
                                        - m_fvCoeffs.Cont.AV[n](j) * bV
                                        - m_fvCoeffs.Cont.AW[t](k) * bW
                                        ) * K;


            // Update U from momentum 
            m_fields[U]( G(iU, jU, kU) ) = bU 
                                         - m_fvCoeffs.Umom.AP[w](iU) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // Update V from momentum
            m_fields[V]( G(iV, jV, kV) ) = bV 
                                         - m_fvCoeffs.Vmom.AP[s](jV) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);

            // Update W from momentum
            m_fields[W]( G(iW, jW, kW) ) = bW 
                                         - m_fvCoeffs.Wmom.AP[b](kW) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);

        }


    private:

        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;


        // Return the index by which the momentum is staggered
        static constexpr intType UoffsetCoupled(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return 0;
            } else if ( coeff == TC::e ) {
                return 1;
            } else if ( coeff == TC::w ) {
                return -1;
            }
        }

        static constexpr intType VoffsetCoupled(TC coeff)
        {
            if        ( coeff == TC::p ) {
                return 0;
            } else if ( coeff == TC::n ) {
                return 1;
            } else if ( coeff == TC::s ) {
                return -1;
            }
        }

        static constexpr intType WoffsetCoupled(TC coeff)
        {
            if        ( coeff == TC::p ) {
                return 0;
            } else if ( coeff == TC::t ) {
                return 1;
            } else if ( coeff == TC::b ) {
                return -1;
            }
        }


        // Returns staggering index for left and right uncoupled coordinates
        static constexpr intType OffsetLeft(intType offset)
        {
            if ( offset == 1 ) {
                return 0;
            } else if ( offset == 0 ) {
                return -1;
            } else if ( offset == -1 ) {
                return -1;
            }
        }

        static constexpr intType OffsetRight(intType offset)
        {
            if ( offset == 1 ) {
                return 1;
            } else if ( offset == 0 ) {
                return 1;
            } else if ( offset == -1 ) {
                return 0;
            }
        }


        // Return coupled, left, and right compass corrdinates if non staggered coordinates
        static constexpr TC UcompassCoupled(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::p;

            } else if ( coeff == TC::e ) {
                return TC::w;

            } else if ( coeff == TC::w ) {
                return TC::e;

            }
        }

        static constexpr TC UcompassLeft(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::w;

            } else if ( coeff == TC::e ) {
                return TC::w;

            } else if ( coeff == TC::w ) {
                return TC::p;

            }
        }

        static constexpr TC UcompassRight(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::e;

            } else if ( coeff == TC::e ) {
                return TC::p;

            } else if ( coeff == TC::w ) {
                return TC::e;
                
            }
        }


        static constexpr TC VcompassCoupled(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::p;

            } else if ( coeff == TC::n ) {
                return TC::s;

            } else if ( coeff == TC::s ) {
                return TC::n;

            }
        }

        static constexpr TC VcompassLeft(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::s;

            } else if ( coeff == TC::n ) {
                return TC::s;

            } else if ( coeff == TC::s ) {
                return TC::p;

            }
        }

        static constexpr TC VcompassRight(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::n;

            } else if ( coeff == TC::n ) {
                return TC::p;

            } else if ( coeff == TC::s ) {
                return TC::n;
                
            }
        }


        static constexpr TC WcompassCoupled(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::p;

            } else if ( coeff == TC::t ) {
                return TC::b;

            } else if ( coeff == TC::b ) {
                return TC::t;

            }
        }

        static constexpr TC WcompassLeft(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::b;

            } else if ( coeff == TC::t ) {
                return TC::b;

            } else if ( coeff == TC::b ) {
                return TC::p;

            }
        }

        static constexpr TC WcompassRight(TC coeff)
        {   
            if        ( coeff == TC::p ) {
                return TC::t;

            } else if ( coeff == TC::t ) {
                return TC::p;

            } else if ( coeff == TC::b ) {
                return TC::t;
                
            }
        }
};



// Solves line
class LineSolver
{
    public:

        LineSolver( ArrayAllocator<Fields, array3D> &fields,
                    FVCoefficients &fvCoeffs,
                    const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_fvCoeffs( fvCoeffs ),
        m_maxIterations( lineSolverSettings.maxIterations ),
        m_maxResiduals( lineSolverSettings.maxResiduals )
        {}

        void SolveLine(const intType j, const intType k)
        {
            BlockSolver blockSolver(m_fields, m_fvCoeffs);
            using enum TransportCoefficients::ENUMDATA;
            blockSolver.UpdateBlock<e, n, t>(1, 2, 3);

        }


    private:
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;

};


// Solves plane
class PlaneSolver
{
    public:

        PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
                     FVCoefficients &fvCoeffs,
                     const InputData::PlaneSolverSettings &planeSolverSettings,
                     const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_fvCoeffs( fvCoeffs ),
        m_maxIterations( planeSolverSettings.maxIterations ),
        m_maxResiduals( planeSolverSettings.maxResiduals ),
        m_lineSolver( fields, fvCoeffs, lineSolverSettings )
        {}

        void SolvePlane(const intType k)
        {

        }


    private:
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        LineSolver m_lineSolver;
};




void CFD::SweepSolve(ArrayAllocator<CFD::Fields, CFD::array3D>  &fields, 
                     const Mesh &mesh, 
                     const InputData &inputData) 
{
    // Extract from input data
    const InputData::PlaneSweepSettings  planeSweepSettings  = inputData.planeSweepSettings;
    const InputData::PlaneSolverSettings planeSolverSettings = inputData.planeSolverSettings;
    const InputData::LineSolverSettings  lineSolverSettings  = inputData.lineSolverSettings;

    // Initialise
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld( fields );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, inputData );
    
    // Counters and residuals
    intType nOuterIterations, nInnerIterations;
    EnumVector<Fields, floatType> residualsInner, residualsOuter, residualsOuterInitial, residualsInnerInitial;
    std::vector< EnumVector<Fields, floatType> > residualsHistory( planeSweepSettings.maxOuterIterations );


    // Outer iterations
    nOuterIterations = 0;
    while ( nOuterIterations < planeSweepSettings.maxOuterIterations ) 
    {
        
        // Inner iterations
        nInnerIterations = 0;
        while ( nInnerIterations < planeSweepSettings.maxInnerIterations ) {
        
            // Coupled sweep for RED nodes in +z direction



            // Coupled sweep for BLACK nodes in -z direction


            // Update residuals
            if ( nInnerIterations == 0 ) {
                UpdateResiduals(residualsInnerInitial, fields, fieldsOld);
            }
            UpdateResiduals(residualsInner, fields, fieldsOld);
            nInnerIterations++;


            // Check residual tolerence
            if ( MetResidualTolerence(residualsInner, planeSweepSettings.maxInnerResiduals) ) {
                break;
            }
        }
        

        // Update residuals
        if ( nOuterIterations == 0 ) {
            UpdateResiduals(residualsOuterInitial, fields, fieldsOld);
        }
        UpdateResiduals(residualsOuter, fields, fieldsOld, residualsOuterInitial);
        residualsHistory[nOuterIterations] = residualsOuter;
        nOuterIterations++;


        // Check residual tolerence
        if ( MetResidualTolerence(residualsOuter, planeSweepSettings.maxOuterResiduals) ) {
            break;
        }


        // Update nonlinear coefficients
        UpdateFaceVelocities( faceVelocities, mesh, fields, inputData );
        UpdateFVCoefficients( fvCoeffs, mesh, fields, inputData );

        // Update fields
        fieldsOld = fields;
    }


    // Strip unused data
    residualsHistory.resize( nOuterIterations );

        
}
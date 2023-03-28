#include "FiniteVolumeFunctions.h"

// Implementation file for finite volume coefficient structure and update functions


namespace CFD 
{

using C = TransportCoefficients::ENUMDATA;

FVCoefficients::FVCoefficients(const indexVector3 &dims) :
    auu({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
    avv({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
    aww({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),

    aup({C::p, C::e, C::w}, dims(0)),
    avp({C::p, C::n, C::s}, dims(1)),
    awp({C::p, C::t, C::b}, dims(2)),

    acu({C::p, C::e, C::w}, dims(0)),
    acv({C::p, C::n, C::s}, dims(1)),
    acw({C::p, C::t, C::b}, dims(2)),
    
    acp({C::p, C::n, C::e, C::s, C::w, C::t, C::b, C::nn, C::ee, C::ss, C::ww, C::tt, C::bb}, dims)
    {};

}  // end namespace CFD



namespace
{




}   // end anonymous namespace



namespace CFD 
{


void InitialiseFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields::ENUMDATA> &faceVelocities)
{

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields::ENUMDATA> &faceVelocities)
{

}


}   // end namespace CFD




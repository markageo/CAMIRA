#ifndef VTK_READER
#define VTK_READER

#include "../Types.h"

#include <string>


namespace VTK
{
    
using namespace CFD;

struct FieldFileData {
    EnumVector<Axis, Tensor1D> cellFaces;
    FieldData<Tensor3D> cellFields;
    FieldData<Tensor3D> vertexFields;
};


FieldFileData ReadVTKFields(std::string &);


} // end namespace VTK


#endif // VTK_READER
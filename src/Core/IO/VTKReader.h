#ifndef CAMIRA_VTK_READER
#define CAMIRA_VTK_READER

#include "Core/Types.h"

#include <string>


namespace VTK
{
    
using namespace CAMIRA;

struct FieldFileData {
    EnumVector<Axis, Tensor1D> cellFaces;
    FieldData<Tensor3D> cellFields;
    FieldData<Tensor3D> vertexFields;
};


FieldFileData ReadVTKFields(const std::string &);


} // end namespace VTK


#endif // CAMIRA_VTK_READER
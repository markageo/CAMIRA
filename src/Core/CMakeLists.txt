target_sources( libs 
    PRIVATE 
    PRIVATE 
        Macros.h
        Types.h
        ArrayIndexConversions.h
        FVLookups.h
        FVTools.h )

add_subdirectory(Geometry)
add_subdirectory(IO)
add_subdirectory(Mesh)
add_subdirectory(AxisTransformationMap)



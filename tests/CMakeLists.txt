option (CFD_TESTS "Create CFD debugging tests executable" off)
if (CFD_TESTS)
	message(STATUS "CFD debugging tests executable will be created")

    # Library used by tests
    add_library(test_libs TestFunctions.cpp)
    target_link_libraries(test_libs PUBLIC compiler_flags)

    # tests executable
    add_executable(test test.cpp)
    target_link_libraries(test PUBLIC model)
    target_link_libraries(test_libs PUBLIC model)
    target_link_libraries(test PUBLIC test_libs)
    
endif()

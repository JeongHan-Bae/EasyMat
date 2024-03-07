cmake_minimum_required(VERSION 3.24)
project(EasyMat)

# Add source files
set(SOURCE_FILES src/jh_matrix.cpp)

# Create a dynamic library from the source files
add_library(EasyMat SHARED ${SOURCE_FILES})

# Specify the include directories
target_include_directories(EasyMat PUBLIC /include/EasyMat)

set(DLL_DIR ${CMAKE_SOURCE_DIR}/lib)

set(DLL_FILES libEasyMat.dll)

set(DEST_DIR ${CMAKE_BINARY_DIR})

foreach (DLL_FILE ${DLL_FILES})
    add_custom_command(TARGET EasyMat POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${DLL_DIR}/${DLL_FILE}"
            "${DEST_DIR}")
endforeach ()
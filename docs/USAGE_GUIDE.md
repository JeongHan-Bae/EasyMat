# EasyMat Usage Guide

EasyMat is a simple and efficient C++ template library for working with matrices. This guide will walk you through the steps to use EasyMat in your C++ project with CMake.

## Prerequisites

Before getting started, ensure you have the following installed on your system:

- CMake (version 3.27 or higher)
- C++ compiler supporting C++17 standard

## Integration with CMake

1. **Download EasyMat**: Download the EasyMat library and extract it to a directory on your system.

2. **Create a CMake Project**: Create a new directory for your CMake project and navigate to it in your terminal.

3. **CMakeLists.txt Configuration**:

   Create a `CMakeLists.txt` file in your project directory and add the following lines:

   ```cmake
   cmake_minimum_required(VERSION 3.27)
   project(YourProject)

   set(CMAKE_CXX_STANDARD 17)

   # Add the source files of your project
   add_executable(YourProject main.cpp)

   # Set the path to the EasyMat header files
   set(EASY_MAT_DIR "path_to_easy_mat_include_directory/include/EasyMat")

   # Include EasyMat headers
   target_include_directories(YourProject PRIVATE "${EASY_MAT_DIR}")

   # Link against the EasyMat library
   set(EASY_MAT_LIB_DIR "${EASY_MAT_DIR}/../../lib")
   target_link_libraries(YourProject PRIVATE "${EASY_MAT_LIB_DIR}/libEasyMat.dll")

   # Set the working directory for the executable
   set_target_properties(YourProject PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/YourProject")

   # Copy the EasyMat DLL file to the output directory
   add_custom_command(TARGET YourProject POST_BUILD
           COMMAND ${CMAKE_COMMAND} -E copy_if_different
           "${EASY_MAT_LIB_DIR}/libEasyMat.dll"
           "$<TARGET_FILE_DIR:YourProject>"
           COMMENT "Copying libEasyMat.dll to YourProject output directory"
   )
   ```

   Replace `YourProject` and `path_to_easy_mat_include_directory` with appropriate values.

4. **Build your Project**: Run CMake to generate build files and then build your project using your preferred build tool.

## Using EasyMat in Your Code

Once you've integrated EasyMat into your project, you can start using it in your C++ code. Here's an example of how to use EasyMat to create a matrix with an initializer list:

```cpp
#include <iostream>
#include "jh_matrix.h" // Include EasyMat header file

using namespace jh; // Use the jh namespace

int main() {
    // Create a matrix using an initializer list
    matrix<int> mat = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    // Print the matrix size
    std::cout << "Matrix size: " << mat.get_row_count() << "x" << mat.get_col_count() << std::endl;

    return 0;
}
```

This example demonstrates how to create a matrix using EasyMat with an initializer list.

That's it! You've successfully integrated EasyMat into your C++ project and used it to create matrices.
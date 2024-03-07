# EasyMat API Documentation

## Overview

This API documentation provides details on the usage of the EasyMat library, a C++ template library for working with 2D matrices, developed by JeongHan Bae.

## Classes

### `matrix<T>`

A generic 2D matrix class.

The template parameter `T` can represent integer types such as `int`, `long`, or `long long int`, 
as well as floating-point types like `float`, `double`, or `long double`.



#### Constructors

- `matrix(size_t rows, size_t cols)`: 

  Constructs an initialized matrix with specified dimensions.

  This constructor initializes all elements of the matrix to a default value of the template type. 

  The matrix is created with the specified number of rows and columns.

- `matrix(const matrix<T>& other)`: Copy constructor.

  This constructor creates a new matrix by copying another matrix of the same type. 
  
  Additionally, you can create a matrix of a specific type by copying an instance of a corresponding type.

  **It's worth noting that** while creating a base matrix by copying a specific matrix is allowed, it's not recommended, as the specialized types inherit all functions from the base class.

- `matrix(std::initializer_list<std::initializer_list<T>> init)`: 

  Constructs a matrix from an initializer list of initializer lists.

  This constructor takes an initializer list of initializer lists to initialize the matrix.

  Any elements not explicitly specified in the initializer lists are filled with 0 values of the template type.

#### Methods

- `size_t get_row_count() const`: Gets the number of rows in the matrix.
- `size_t get_col_count() const`: Gets the number of columns in the matrix.
- `T& operator()(size_t row, size_t col)`: Accesses an element in the matrix for modification.
- `const T& operator()(size_t row, size_t col) const`: Accesses an element in the matrix for reading.
- `matrix<T> operator+(const matrix<T>& other) const`: Adds two matrices element-wise.
- `matrix<T>& operator+=(const matrix<T>& other)`: Adds another matrix to this matrix and assigns the result.
- `matrix<T> operator-(const matrix<T>& other) const`: Subtracts another matrix from this matrix element-wise.
- `matrix<T>& operator-=(const matrix<T>& other)`: Subtracts another matrix from this matrix and assigns the result.
- `matrix<T> operator*(const matrix<T>& other) const`: Multiplies two matrices.
- `matrix<T>& operator*=(const matrix<T>& other)`: Multiplies this matrix with another and assigns the result.


It's important to note that any result obtained from operations such as addition (`+`), subtraction (`-`), or multiplication (`*`) 
will be treated as the base class `matrix<T>`. When using compound assignment operators (`+=`, `-=`), the class won't change, 
but for `*=`, the specified matrix will be reset to the base (to avoid type mismatch).

---

#### `static Derived& static_cast_to(matrix<T>& mat)`: 

  This static method allows for casting the base matrix to a specified matrix of the same type.

  Example usage: 

  `matrix<T>::static_cast_to<Derived>(matrix<T>& mat);`

  It's important to note that if there's a size mismatch, an exception will be thrown.

---

#### `static void resize(matrix<T>& mat, size_t new_rows, size_t new_cols)`

Resizes the given matrix to new dimensions.

- `mat`: Matrix to resize.
- `new_rows`: New number of rows.
- `new_cols`: New number of columns.

Throws `invalid_argument` if the new dimensions don't have the same total number of elements.

#### `matrix<T> resize(size_t new_rows, size_t new_cols) const`

Create a matrix of resized dimensions.

- `new_rows`: New number of rows.
- `new_cols`: New number of columns.

Returns a resized matrix.

---

#### `static void transpose(matrix<T>& mat)`

Transposes the given matrix.

- `mat`: Matrix to transpose.

#### `matrix<T> transpose() const`

Gives a transpose of the matrix.

---

**Notice that** the dynamic methods of transpose and resize will not change the original matrix, but return a new matrix, while the class static function will modify the aimed matrix directly and change its type to the base case forcefully.

---


Returns a transposed matrix.
### `horizontal_matrix<T>`

A specialized horizontal matrix class inheriting from `matrix<T>`.

#### Constructors

- `horizontal_matrix(size_t cols)`: Constructs a horizontal matrix with specified number of columns.

#### Methods

- `horizontal_matrix<T>& operator=(const matrix<T>& other)`: Assignment operator for horizontal matrices.

### `vertical_matrix<T>`

A specialized vertical matrix class inheriting from `matrix<T>`.

#### Constructors

- `vertical_matrix(size_t rows)`: Constructs a vertical matrix with specified number of rows.

#### Methods

- `vertical_matrix<T>& operator=(const matrix<T>& other)`: Assignment operator for vertical matrices.

### `square_matrix<T>`

A specialized square matrix class inheriting from `matrix<T>`.

#### Constructors

- `square_matrix(size_t size)`: Constructs a square matrix with specified size.

#### Methods

- `static square_matrix<T>& identity(size_t size)`: Generates an identity square matrix of specified size.
- `static square_matrix<T>& diagonal(T value, size_t size)`: Generates a square matrix with diagonal elements set to a specified value.
- `square_matrix<T>& operator=(const matrix<T>& other)`: Assignment operator for square matrices.
- `square_matrix<T> inverse() const`: Performs matrix inversion.
- `square_matrix<T> operator/(const square_matrix<T>& other) const`: Performs matrix division.
- `square_matrix<T>& operator/=(const square_matrix<T>& other)`: Performs matrix division and assigns the result.

## Usage

To use the EasyMat library in your C++ project, include the header file `jh_matrix.h` and use the namespace `jh`. Here's an example demonstrating the basic usage:

```cpp
#include "jh_matrix.h"

using namespace jh;

int main() {
    // Create a 3x3 matrix
    matrix<int> mat1(3, 3);

    // Initialize the matrix with values using an initializer list
    matrix<int> mat2 = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

    // Perform matrix addition
    matrix<int> result = mat1 + mat2;

    return 0;
}
```

## Fibonacci Calculation Example

The `fibonacci` function demonstrates how to efficiently calculate Fibonacci numbers using matrix exponentiation:

```cpp
#include "jh_matrix.h"
#include <iostream>

using namespace jh;

template <typename T>
T fibonacci(size_t n, T init1 = 0, T init2 = 1) {
    // Define the matrix M for Fibonacci calculation
    square_matrix<T> M(2);
    M(0, 0) = M(0, 1) = M(1, 0) = 1;

    vertical_matrix<T> V(2);
    V(0, 0) = init2;
    V(1, 0) = init1;

    // Compute M^n
    auto result_matrix = M.pow(n) * V;

    // Return the element at (0, 1) which represents the nth Fibonacci number
    return result_matrix(1, 0);
}

int main() {
    // Calculate the 20th Fibonacci number with initial values 3 and 5
    long fib_20 = fibonacci<long>(20, 3, 5);
    std::cout << "20th Fibonacci number: " << fib_20

 << std::endl;
    return 0;
}
```

Ensure that you have integrated EasyMat into your project using CMake as described in the Integration with CMake section.

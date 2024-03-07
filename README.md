# EasyMat

EasyMat is a C++ template library developed by JeongHan Bae for working with 2D matrices. It provides a simple and efficient way to handle matrix operations in C++ programs.

## Project Structure

```markdown
EasyMatLib
│
├── include
│   └── EasyMat
│       └── jh_matrix.h
│
├── src
│   └── jh_matrix.cpp
│
├── lib
│   └── (DLL files)
│
├── CMakeLists.txt
│
├── README.md
│
└── docs
    ├── API_DOCUMENTATION.md
    └── USAGE_GUIDE.md
```
The project structure organizes source code, libraries, and documentation under a cohesive hierarchy.

## Documentation

  Detailed library documentation is available in two parts:

- [API Documentation](docs/API_DOCUMENTATION.md)
- [Usage Guide](docs/USAGE_GUIDE.md)

## Key Features

- Generic 2D matrix class supporting various numeric types (integer and floating-point).
- Constructors for initializing matrices with specified dimensions or using initializer lists.
- Basic matrix operations including addition, subtraction, multiplication, and more.
- Specialized matrix classes for horizontal, vertical, and square matrices.
- Functions for resizing and transposing matrices.
- Efficient implementation for calculating Fibonacci numbers using matrix exponentiation.

## About the Author

JeongHan Bae is a software engineer with a passion for developing libraries and tools to simplify complex programming tasks. With years of experience in C++ development, JeongHan created EasyMat to provide a user-friendly solution for matrix operations in C++ programs.

## About the Name "EasyMat"

The name "EasyMat" reflects the simplicity and ease of use of this matrix library. "Mat" is an abbreviation for both "matrix" and "mathematics," highlighting its focus on mathematical operations involving matrices. Designed with simplicity in mind, EasyMat aims to make matrix operations straightforward and accessible for C++ developers.

## License

This project is licensed under the [MIT License](LICENSE).
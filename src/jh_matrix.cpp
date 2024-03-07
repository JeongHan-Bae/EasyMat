#include "../include/EasyMat/jh_matrix.h"

#include <cstring>
#include <numeric>
#include <stdexcept>

namespace jh {
    template <>
    [[nodiscard]] square_matrix<int> square_matrix<int>::inverse() const = delete;

    template <>
    [[nodiscard]] square_matrix<long> square_matrix<long>::inverse() const = delete;

    template <>
    [[nodiscard]] square_matrix<long long int>
    square_matrix<long long int>::inverse() const = delete;

    template <>
    [[nodiscard]] square_matrix<int>
    square_matrix<int>::operator/(const square_matrix<int>& other) const = delete;

    template <>
    [[nodiscard]] square_matrix<long>
    square_matrix<long>::operator/(const square_matrix<long>& other) const = delete;

    template <>
    [[nodiscard]] square_matrix<long long int>
    square_matrix<long long int>::operator/(
            const square_matrix<long long int>& other) const = delete;

    template <>
    [[nodiscard]] square_matrix<int>&
    square_matrix<int>::operator/=(const square_matrix<int>& other) = delete;

    template <>
    [[nodiscard]] square_matrix<long>&
    square_matrix<long>::operator/=(const square_matrix<long>& other) = delete;

    template <>
    [[nodiscard]] square_matrix<long long int>&
    square_matrix<long long int>::operator/=(
            const square_matrix<long long int>& other) = delete;

    template class matrix<int>;
    template class horizontal_matrix<int>;
    template class vertical_matrix<int>;
    template class square_matrix<int>;

    template class matrix<long>;
    template class horizontal_matrix<long>;
    template class vertical_matrix<long>;
    template class square_matrix<long>;

    template class matrix<long long int>;
    template class horizontal_matrix<long long int>;
    template class vertical_matrix<long long int>;
    template class square_matrix<long long int>;

    template class matrix<float>;
    template class horizontal_matrix<float>;
    template class vertical_matrix<float>;
    template class square_matrix<float>;

    template class matrix<double>;
    template class horizontal_matrix<double>;
    template class vertical_matrix<double>;
    template class square_matrix<double>;

    template class matrix<long double>;
    template class horizontal_matrix<long double>;
    template class vertical_matrix<long double>;
    template class square_matrix<long double>;

    template <typename T>
    matrix<T>::matrix(size_t rows, size_t cols)
            : row_count(rows), col_count(cols), mat_array(new T[rows * cols]) {
        std::memset(mat_array.get(), 0, row_count * col_count * sizeof(T));
    }

    template <typename T>
    matrix<T>::matrix(const matrix<T>& other)
            : row_count(other.row_count), col_count(other.col_count),
              mat_array(new T[row_count * col_count]) {
        std::memcpy(mat_array.get(), other.mat_array.get(),
                    row_count * col_count * sizeof(T));
    }

    template <typename T>
    [[maybe_unused]] matrix<T>::matrix(
            std::initializer_list<std::initializer_list<T>> init) {
        // Get the size of the initializer list
        row_count = init.size();
        col_count = 0;

// Find the largest size of inner initializer lists
#pragma omp parallel for reduction(max : col_count)
        for (auto outer_it = init.begin(); outer_it != init.end(); ++outer_it) {
            col_count = std::max(col_count, outer_it->size());
        }

        // Allocate memory for the matrix and fill with zeros
        mat_array = std::make_unique<T[]>(row_count * col_count);
        std::fill(mat_array.get(), mat_array.get() + row_count * col_count, T{});

        // Fill the matrix with the provided values
        size_t i = 0;
        for (auto outer_it = init.begin(); outer_it != init.end(); ++outer_it) {
            size_t j = 0;
            for (auto inner_it = outer_it->begin(); inner_it != outer_it->end();
                 ++inner_it) {
                mat_array[i * col_count + j] = *inner_it;
                ++j;
            }
            ++i;
        }
    }

    template matrix<int>::matrix(
            std::initializer_list<std::initializer_list<int>> init);

    template matrix<long>::matrix(
            std::initializer_list<std::initializer_list<long>> init);

    template matrix<long long int>::matrix(
            std::initializer_list<std::initializer_list<long long int>> init);

    template matrix<float>::matrix(
            std::initializer_list<std::initializer_list<float>> init);

    template matrix<double>::matrix(
            std::initializer_list<std::initializer_list<double>> init);

    template matrix<long double>::matrix(
            std::initializer_list<std::initializer_list<long double>> init);

    template <typename T>
    template <typename O>
    [[maybe_unused]] matrix<T> matrix<T>::force_convert(const matrix<O>& mat) {
        size_t rows = mat.get_row_count();
        size_t cols = mat.get_col_count();
        matrix<T> new_mat(rows, cols);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                new_mat(i, j) = static_cast<T>(mat(i, j));
            }
        }

        return new_mat;
    }

    template matrix<double> matrix<double>::force_convert(const matrix<int>& mat);
    template matrix<double> matrix<double>::force_convert(const matrix<long>& mat);
    template matrix<double>
    matrix<double>::force_convert(const matrix<long long int>& mat);
    template matrix<double> matrix<double>::force_convert(const matrix<float>& mat);
    template matrix<double>
    matrix<double>::force_convert(const matrix<long double>& mat);

    template matrix<float> matrix<float>::force_convert(const matrix<int>& mat);
    template matrix<float> matrix<float>::force_convert(const matrix<long>& mat);
    template matrix<float>
    matrix<float>::force_convert(const matrix<long long int>& mat);
    template matrix<float> matrix<float>::force_convert(const matrix<double>& mat);
    template matrix<float>
    matrix<float>::force_convert(const matrix<long double>& mat);

    template matrix<long double>
    matrix<long double>::force_convert(const matrix<int>& mat);
    template matrix<long double>
    matrix<long double>::force_convert(const matrix<long>& mat);
    template matrix<long double>
    matrix<long double>::force_convert(const matrix<long long int>& mat);
    template matrix<long double>
    matrix<long double>::force_convert(const matrix<float>& mat);
    template matrix<long double>
    matrix<long double>::force_convert(const matrix<double>& mat);

    template matrix<int> matrix<int>::force_convert(const matrix<long>& mat);
    template matrix<int>
    matrix<int>::force_convert(const matrix<long long int>& mat);
    template matrix<int> matrix<int>::force_convert(const matrix<float>& mat);
    template matrix<int> matrix<int>::force_convert(const matrix<double>& mat);
    template matrix<int> matrix<int>::force_convert(const matrix<long double>& mat);

    template matrix<long> matrix<long>::force_convert(const matrix<int>& mat);
    template matrix<long>
    matrix<long>::force_convert(const matrix<long long int>& mat);
    template matrix<long> matrix<long>::force_convert(const matrix<float>& mat);
    template matrix<long> matrix<long>::force_convert(const matrix<double>& mat);
    template matrix<long>
    matrix<long>::force_convert(const matrix<long double>& mat);

    template matrix<long long int>
    matrix<long long int>::force_convert(const matrix<int>& mat);
    template matrix<long long int>
    matrix<long long int>::force_convert(const matrix<long>& mat);
    template matrix<long long int>
    matrix<long long int>::force_convert(const matrix<float>& mat);
    template matrix<long long int>
    matrix<long long int>::force_convert(const matrix<double>& mat);
    template matrix<long long int>
    matrix<long long int>::force_convert(const matrix<long double>& mat);

    template <typename T> size_t matrix<T>::get_row_count() const {
        return row_count;
    }

    template <typename T> size_t matrix<T>::get_col_count() const {
        return col_count;
    }

    template <typename T> T& matrix<T>::operator()(size_t row, size_t col) {
        return mat_array[row * col_count + col];
    }

    template <typename T>
    const T& matrix<T>::operator()(size_t row, size_t col) const {
        return mat_array[row * col_count + col];
    }

    template <typename T>
    matrix<T> matrix<T>::operator+(const matrix<T>& other) const {
        if (row_count != other.row_count || col_count != other.col_count) {
            throw std::invalid_argument(
                    "Matrix dimensions are not compatible for addition.");
        }

        matrix<T> result(row_count, col_count);
#pragma omp parallel for
        for (size_t i = 0; i < row_count * col_count; ++i) {
            result.mat_array[i] = mat_array[i] + other.mat_array[i];
        }
        return result;
    }

    template <typename T> matrix<T>& matrix<T>::operator+=(const matrix<T>& other) {
#pragma omp parallel for
        for (size_t i = 0; i < row_count * col_count; ++i) {
            mat_array[i] += other.mat_array[i];
        }
        return *this;
    }

    template <typename T>
    matrix<T> matrix<T>::operator-(const matrix<T>& other) const {
        if (row_count != other.row_count || col_count != other.col_count) {
            throw std::invalid_argument(
                    "Matrix dimensions are not compatible for subtraction.");
        }

        matrix<T> result(row_count, col_count);

#pragma omp parallel for
        for (size_t i = 0; i < row_count * col_count; ++i) {
            result.mat_array[i] = mat_array[i] - other.mat_array[i];
        }

        return result;
    }

    template <typename T> matrix<T>& matrix<T>::operator-=(const matrix<T>& other) {
        if (this == &other) {
            // If matrices are the same, set elements to zero using memset
            memset(mat_array.get(), 0, row_count * col_count * sizeof(T));
        } else {
            // If matrices are different, subtract elements as usual
#pragma omp parallel for
            for (size_t i = 0; i < row_count * col_count; ++i) {
                mat_array[i] -= other.mat_array[i];
            }
        }
        return *this;
    }

    template <typename T>
    matrix<T> matrix<T>::operator*(const matrix<T>& other) const {
        if (col_count != other.row_count) {
            throw std::invalid_argument(
                    "Matrix dimensions are not compatible for multiplication.");
        }

        matrix<T> result(row_count, other.col_count);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < row_count; ++i) {
            for (size_t j = 0; j < other.col_count; ++j) {
                T sum = 0;
                for (size_t k = 0; k < col_count; ++k) {
                    sum += (*this)(i, k) * other(k, j);
                }
                result(i, j) = sum;
            }
        }
        return result;
    }

    template <typename T> matrix<T>& matrix<T>::operator*=(const matrix<T>& other) {
        auto temp = (*this) * other;
        this->row_count = other.row_count;
        this->col_count = other.col_count;
        this->mat_array = std::move(temp.mat_array);
        return static_cast<matrix<T>&>(*this);
    }

    template <typename T>
    [[maybe_unused]] void matrix<T>::resize(matrix<T>& mat, size_t new_rows,
                                            size_t new_cols) {
        if (new_rows * new_cols != mat.row_count * mat.col_count) {
            throw std::invalid_argument(
                    "New dimensions must have the same total number of elements.");
        }

        // Update the matrix attributes with the new dimensions and data
        // As mat_array is 1D, no need to modify
        mat.row_count = new_rows;
        mat.col_count = new_cols;

        // Check if the input matrix is of any derived class type and cast it to
        // matrix if necessary
        if constexpr (!std::is_same_v<decltype(mat), matrix<T>&>) {
            mat = static_cast<matrix<T>&>(mat);
        }
    }

    template <typename T>
    [[maybe_unused]] matrix<T> matrix<T>::resize(size_t new_rows,
                                                 size_t new_cols) const {

        if (new_rows * new_cols != this->row_count * this->col_count) {
            throw std::invalid_argument(
                    "New dimensions must have the same total number of elements.");
        }
        matrix<T> new_mat(new_rows, new_cols);

        std::memcpy(new_mat.mat_array.get(), this->mat_array.get(),
                    new_rows * new_cols * sizeof(T));

        return std::move(new_mat);
    }

    template <typename T>
    [[maybe_unused]] void matrix<T>::transpose(matrix<T>& mat) {

        std::unique_ptr<T[]> new_mat_array(new T[mat.row_count * mat.col_count]);
#pragma omp parallel for collapse(2)
        for (size_t j = 0; j < mat.col_count; ++j) {
            for (size_t i = 0; i < mat.row_count; ++i) {
                new_mat_array[j * mat.row_count + i] =
                        mat.mat_array[i * mat.col_count + j];
            }
        }
        std::swap(mat.row_count, mat.col_count);
        mat.mat_array = std::move(new_mat_array);

        // Check if the resulting matrix is of any derived class type other than
        // square_matrix and cast it to matrix if necessary
        if constexpr (!std::is_same_v<decltype(mat), matrix<T>&> &&
                      !std::is_same_v<decltype(mat), square_matrix<T>&>) {
            mat = static_cast<matrix<T>&>(mat);
        }
    }

    template <typename T> [[maybe_unused]] matrix<T> matrix<T>::transpose() const {
        // If the matrix is square, create a square_matrix for transposition
        if (row_count == col_count) {
            square_matrix<T> transposed_square_mat(row_count);
#pragma omp parallel for collapse(2)
            for (size_t i = 0; i < row_count; ++i) {
                for (size_t j = i + 1; j < col_count; ++j) {
                    transposed_square_mat.mat_array[i * row_count + j] =
                            mat_array[j * row_count + i];
                    transposed_square_mat.mat_array[j * row_count + i] =
                            mat_array[i * row_count + j];
                }
                // Diagonal elements remain the same
                transposed_square_mat.mat_array[i * row_count + i] =
                        mat_array[i * row_count + i];
            }
            return transposed_square_mat;
        } else {
            // If it's not square, create a new matrix for transposition
            matrix<T> transposed_mat(col_count, row_count);
#pragma omp parallel for collapse(2)
            for (size_t i = 0; i < row_count; ++i) {
                for (size_t j = 0; j < col_count; ++j) {
                    transposed_mat.mat_array[j * row_count + i] =
                            mat_array[i * col_count + j];
                }
            }
            return transposed_mat;
        }
    }

    template <typename T>
    template <typename Derived>
    [[maybe_unused]] Derived& matrix<T>::static_cast_to(matrix<T>& mat) {
        if constexpr (std::is_same_v<Derived, horizontal_matrix<T>>) {
            if (mat.get_row_count() != 1) {
                throw std::invalid_argument(
                        "Invalid cast: row count must be 1 for horizontal_matrix.");
            }
        } else if constexpr (std::is_same_v<Derived, vertical_matrix<T>>) {
            if (mat.get_col_count() != 1) {
                throw std::invalid_argument(
                        "Invalid cast: column count must be 1 for vertical_matrix.");
            }
        } else if constexpr (std::is_same_v<Derived, square_matrix<T>>) {
            if (mat.get_row_count() != mat.get_col_count()) {
                throw std::invalid_argument("Invalid cast: row count must be equal "
                                            "to column count for square_matrix.");
            }
        }
        return static_cast<Derived&>(mat);
    }

    template horizontal_matrix<int>&
    matrix<int>::static_cast_to<horizontal_matrix<int>>(matrix<int>& mat);
    template horizontal_matrix<long>&
    matrix<long>::static_cast_to<horizontal_matrix<long>>(matrix<long>& mat);
    template horizontal_matrix<long long int>&
    matrix<long long int>::static_cast_to<horizontal_matrix<long long int>>(
            matrix<long long int>& mat);
    template horizontal_matrix<float>&
    matrix<float>::static_cast_to<horizontal_matrix<float>>(matrix<float>& mat);
    template horizontal_matrix<double>&
    matrix<double>::static_cast_to<horizontal_matrix<double>>(matrix<double>& mat);
    template horizontal_matrix<long double>&
    matrix<long double>::static_cast_to<horizontal_matrix<long double>>(
            matrix<long double>& mat);

    template vertical_matrix<int>&
    matrix<int>::static_cast_to<vertical_matrix<int>>(matrix<int>& mat);
    template vertical_matrix<long>&
    matrix<long>::static_cast_to<vertical_matrix<long>>(matrix<long>& mat);
    template vertical_matrix<long long int>&
    matrix<long long int>::static_cast_to<vertical_matrix<long long int>>(
            matrix<long long int>& mat);
    template vertical_matrix<float>&
    matrix<float>::static_cast_to<vertical_matrix<float>>(matrix<float>& mat);
    template vertical_matrix<double>&
    matrix<double>::static_cast_to<vertical_matrix<double>>(matrix<double>& mat);
    template vertical_matrix<long double>&
    matrix<long double>::static_cast_to<vertical_matrix<long double>>(
            matrix<long double>& mat);

    template square_matrix<int>&
    matrix<int>::static_cast_to<square_matrix<int>>(matrix<int>& mat);
    template square_matrix<long>&
    matrix<long>::static_cast_to<square_matrix<long>>(matrix<long>& mat);
    template square_matrix<long long int>&
    matrix<long long int>::static_cast_to<square_matrix<long long int>>(
            matrix<long long int>& mat);
    template square_matrix<float>&
    matrix<float>::static_cast_to<square_matrix<float>>(matrix<float>& mat);
    template square_matrix<double>&
    matrix<double>::static_cast_to<square_matrix<double>>(matrix<double>& mat);
    template square_matrix<long double>&
    matrix<long double>::static_cast_to<square_matrix<long double>>(
            matrix<long double>& mat);

// clang-format off
    template<typename T> // NOLINT
    matrix<T>& matrix<T>::operator=(const matrix<T>& other) // NOLINT
    {
        if (this != &other) {
            // Ensure dimensions match
            assert((row_count == other.row_count && col_count == other.col_count)
                   && "Matrix dimensions must match for assignment.");

            // Perform deep copy
            std::copy(other.mat_array.get(), other.mat_array.get()
                                             + row_count * col_count, mat_array.get());
        }
        return *this;
    }
// clang-format on

    template matrix<int>& matrix<int>::operator=(const matrix<int>& other);
    template matrix<long>& matrix<long>::operator=(const matrix<long>& other);
    template matrix<long long int>&
    matrix<long long int>::operator=(const matrix<long long int>& other);
    template matrix<float>& matrix<float>::operator=(const matrix<float>& other);
    template matrix<double>& matrix<double>::operator=(const matrix<double>& other);
    template matrix<long double>&
    matrix<long double>::operator=(const matrix<long double>& other);

    template <typename T>
    horizontal_matrix<T>::horizontal_matrix(size_t cols) : matrix<T>(1, cols) {}

    template <typename T>
    horizontal_matrix<T>& horizontal_matrix<T>::operator=(const matrix<T>& other) {
        // Call base class assignment operator
        matrix<T>::operator=(other);
        return *this;
    }

    template horizontal_matrix<int>&
    horizontal_matrix<int>::operator=(const matrix<int>& other);
    template horizontal_matrix<long>&
    horizontal_matrix<long>::operator=(const matrix<long>& other);
    template horizontal_matrix<long long int>&
    horizontal_matrix<long long int>::operator=(const matrix<long long int>& other);
    template horizontal_matrix<float>&
    horizontal_matrix<float>::operator=(const matrix<float>& other);
    template horizontal_matrix<double>&
    horizontal_matrix<double>::operator=(const matrix<double>& other);
    template horizontal_matrix<long double>&
    horizontal_matrix<long double>::operator=(const matrix<long double>& other);

    template <typename T>
    vertical_matrix<T>::vertical_matrix(size_t rows) : matrix<T>(rows, 1) {}

    template <typename T>
    vertical_matrix<T>& vertical_matrix<T>::operator=(const matrix<T>& other) {
        // Call base class assignment operator
        matrix<T>::operator=(other);
        return *this;
    }

    template vertical_matrix<int>&
    vertical_matrix<int>::operator=(const matrix<int>& other);
    template vertical_matrix<long>&
    vertical_matrix<long>::operator=(const matrix<long>& other);
    template vertical_matrix<long long int>&
    vertical_matrix<long long int>::operator=(const matrix<long long int>& other);
    template vertical_matrix<float>&
    vertical_matrix<float>::operator=(const matrix<float>& other);
    template vertical_matrix<double>&
    vertical_matrix<double>::operator=(const matrix<double>& other);
    template vertical_matrix<long double>&
    vertical_matrix<long double>::operator=(const matrix<long double>& other);

    template <typename T>
    square_matrix<T>::square_matrix(size_t size) : matrix<T>(size, size) {}

    template <typename T>
    square_matrix<T>& square_matrix<T>::identity(size_t size) {
        auto* result = new square_matrix<T>(size);
        for (size_t i = 0; i < size; ++i) {
            (*result)(i, i) = static_cast<T>(1);
        }
        return *result;
    }

    template <typename T>
    [[maybe_unused]] square_matrix<T>& square_matrix<T>::diagonal(T value,
                                                                  size_t size) {
        auto* result = new square_matrix<T>(size);
        for (size_t i = 0; i < size; ++i) {
            (*result)(i, i) = value;
        }
        return *result;
    }

    template <typename T>
    square_matrix<T>& square_matrix<T>::operator=(const matrix<T>& other) {
        // Call base class assignment operator
        matrix<T>::operator=(other);
        return *this;
    }

    template square_matrix<int>&
    square_matrix<int>::operator=(const matrix<int>& other);
    template square_matrix<long>&
    square_matrix<long>::operator=(const matrix<long>& other);
    template square_matrix<long long int>&
    square_matrix<long long int>::operator=(const matrix<long long int>& other);
    template square_matrix<float>&
    square_matrix<float>::operator=(const matrix<float>& other);
    template square_matrix<double>&
    square_matrix<double>::operator=(const matrix<double>& other);
    template square_matrix<long double>&
    square_matrix<long double>::operator=(const matrix<long double>& other);

    template <typename T>
    template <typename O>
    square_matrix<T> square_matrix<T>::force_convert(const square_matrix<O>& mat) {
        size_t size = mat.get_row_count();
        square_matrix<T> new_mat(size);

#pragma omp parallel for collapse(2)
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                new_mat(i, j) = static_cast<T>(mat(i, j));
            }
        }

        return new_mat;
    }

    template square_matrix<double>
    square_matrix<double>::force_convert(const square_matrix<int>& mat);
    template square_matrix<double>
    square_matrix<double>::force_convert(const square_matrix<long>& mat);
    template square_matrix<double>
    square_matrix<double>::force_convert(const square_matrix<long long int>& mat);
    template square_matrix<double>
    square_matrix<double>::force_convert(const square_matrix<float>& mat);
    template square_matrix<double>
    square_matrix<double>::force_convert(const square_matrix<long double>& mat);

    template square_matrix<float>
    square_matrix<float>::force_convert(const square_matrix<int>& mat);
    template square_matrix<float>
    square_matrix<float>::force_convert(const square_matrix<long>& mat);
    template square_matrix<float>
    square_matrix<float>::force_convert(const square_matrix<long long int>& mat);
    template square_matrix<float>
    square_matrix<float>::force_convert(const square_matrix<double>& mat);
    template square_matrix<float>
    square_matrix<float>::force_convert(const square_matrix<long double>& mat);

    template square_matrix<long double>
    square_matrix<long double>::force_convert(const square_matrix<int>& mat);
    template square_matrix<long double>
    square_matrix<long double>::force_convert(const square_matrix<long>& mat);
    template square_matrix<long double> square_matrix<long double>::force_convert(
            const square_matrix<long long int>& mat);
    template square_matrix<long double>
    square_matrix<long double>::force_convert(const square_matrix<float>& mat);
    template square_matrix<long double>
    square_matrix<long double>::force_convert(const square_matrix<double>& mat);

    template square_matrix<int>
    square_matrix<int>::force_convert(const square_matrix<long>& mat);
    template square_matrix<int>
    square_matrix<int>::force_convert(const square_matrix<long long int>& mat);
    template square_matrix<int>
    square_matrix<int>::force_convert(const square_matrix<float>& mat);
    template square_matrix<int>
    square_matrix<int>::force_convert(const square_matrix<double>& mat);
    template square_matrix<int>
    square_matrix<int>::force_convert(const square_matrix<long double>& mat);

    template square_matrix<long>
    square_matrix<long>::force_convert(const square_matrix<int>& mat);
    template square_matrix<long>
    square_matrix<long>::force_convert(const square_matrix<long long int>& mat);
    template square_matrix<long>
    square_matrix<long>::force_convert(const square_matrix<float>& mat);
    template square_matrix<long>
    square_matrix<long>::force_convert(const square_matrix<double>& mat);
    template square_matrix<long>
    square_matrix<long>::force_convert(const square_matrix<long double>& mat);

    template square_matrix<long long int>
    square_matrix<long long int>::force_convert(const square_matrix<int>& mat);
    template square_matrix<long long int>
    square_matrix<long long int>::force_convert(const square_matrix<long>& mat);
    template square_matrix<long long int>
    square_matrix<long long int>::force_convert(const square_matrix<float>& mat);
    template square_matrix<long long int>
    square_matrix<long long int>::force_convert(const square_matrix<double>& mat);
    template square_matrix<long long int>
    square_matrix<long long int>::force_convert(
            const square_matrix<long double>& mat);

    template <typename T> void square_matrix<T>::square() { *this *= *this; }

    template <typename T>
    square_matrix<T> square_matrix<T>::inverse_helper() const {

        // Check if the matrix is invertible
        if (this->row_count != this->col_count) {
            throw std::invalid_argument(
                    "Matrix is not square, and thus not invertible.");
        }

        size_t size = this->row_count;
        square_matrix<T> helper(*this);
        square_matrix<T> result = square_matrix<T>::identity(size);

        for (size_t i = 0; i < size; ++i) {
            // Find pivot index and value
            size_t pivot_index = i;
            T pivot = helper(i, i);

#pragma omp parallel for reduction(max : pivot) shared(pivot_index)
            for (size_t j = i + 1; j < size; ++j) {
                if (std::abs(helper(j, i)) > std::abs(pivot)) {
#pragma omp critical
                    {
                        pivot = helper(j, i);
                        pivot_index = j;
                    }
                }
            }

            // Check if the matrix is not invertible
            if (pivot == 0.0) {
                throw std::invalid_argument("Matrix is not invertible.");
            }

            // Swap rows of helper and result matrices
            std::swap_ranges(helper.mat_array.get() + i * size,
                             helper.mat_array.get() + (i + 1) * size,
                             helper.mat_array.get() + pivot_index * size);
            std::swap_ranges(result.mat_array.get() + i * size,
                             result.mat_array.get() + (i + 1) * size,
                             result.mat_array.get() + pivot_index * size);

// Normalize rows
#pragma omp parallel for
            for (size_t j = i + 1; j < size; ++j) {
                helper(i, j) /= pivot;
            }
#pragma omp parallel for
            for (size_t j = 0; j < size; ++j) {
                result(i, j) /= pivot;
            }

            // Eliminate elements below pivot
            for (size_t k = i + 1; k < size; ++k) {
                T factor = helper(k, i);
#pragma omp parallel for
                for (size_t j = k; j < size; ++j) {
                    helper(k, j) -= factor * helper(i, j);
                }
#pragma omp parallel for
                for (size_t j = 0; j < size; ++j) {
                    result(k, j) -= factor * result(i, j);
                }
            }
        }

// Back substitution to make diagonal elements 1 and left-below elements 0
#pragma omp parallel for
        for (size_t i = size - 1; i > 0; --i) {
            for (size_t k = i - 1; k < size; --k) {
                T factor = helper(k, i);
                for (size_t j = 0; j < size; ++j) {
                    result(k, j) -= result(i, j) * factor;
                }
            }
        }

        return result;
    }

// Template instantiation for square_matrix<float>
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<int>& mat);
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<long>& mat);
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<long long int>& mat);
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<float>& mat);
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<double>& mat);
    template square_matrix<float>
    square_matrix<float>::inverse_helper(const square_matrix<long double>& mat);

// Template instantiation for square_matrix<double>
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<int>& mat);
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<long>& mat);
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<long long int>& mat);
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<float>& mat);
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<double>& mat);
    template square_matrix<double>
    square_matrix<double>::inverse_helper(const square_matrix<long double>& mat);

// Template instantiation for square_matrix<long double>
    template square_matrix<long double>
    square_matrix<long double>::inverse_helper(const square_matrix<int>& mat);
    template square_matrix<long double>
    square_matrix<long double>::inverse_helper(const square_matrix<long>& mat);
    template square_matrix<long double> square_matrix<long double>::inverse_helper(
            const square_matrix<long long int>& mat);
    template square_matrix<long double>
    square_matrix<long double>::inverse_helper(const square_matrix<float>& mat);
    template square_matrix<long double>
    square_matrix<long double>::inverse_helper(const square_matrix<double>& mat);
    template square_matrix<long double> square_matrix<long double>::inverse_helper(
            const square_matrix<long double>& mat);

    template <typename T>
    template <typename O>
    square_matrix<T> square_matrix<T>::inverse_helper(const square_matrix<O>& mat) {
        try {
            // Convert the input matrix to double
            square_matrix<T> double_mat = square_matrix<T>::force_convert(mat);
            // Compute and return the inverse of the double matrix
            return double_mat.inverse();
        } catch (const std::invalid_argument& e) {
            // Rethrow the exception as an "invalid matrix" exception
            throw std::invalid_argument("Invalid matrix: " + std::string(e.what()));
        }
    }

    template <typename T>
    square_matrix<T>
    square_matrix<T>::operator/(const square_matrix<T>& other) const {
        try {
            square_matrix<T> inverse_other = other.inverse();
            square_matrix<T> result(*this);
            result *= inverse_other;
            return result;
        } catch (const std::invalid_argument& e) {
            throw std::runtime_error("Division by non-invertible matrix.");
        }
    }

// Template instantiation for square_matrix<float>
    template square_matrix<float>
    square_matrix<float>::operator/(const square_matrix<float>& other) const;

// Template instantiation for square_matrix<double>
    template square_matrix<double>
    square_matrix<double>::operator/(const square_matrix<double>& other) const;

// Template instantiation for square_matrix<long double>
    template square_matrix<long double> square_matrix<long double>::operator/(
            const square_matrix<long double>& other) const;

    template <typename T>
    square_matrix<T>& square_matrix<T>::operator/=(const square_matrix<T>& other) {
        auto temp = (*this) / other;
        this->mat_array = std::move(temp.mat_array);
        return *this;
    }

    template <typename T>
    [[maybe_unused]] square_matrix<T>
    square_matrix<T>::pow(const size_t exponent) const {
        size_t parse = (size_t)1 << (sizeof(size_t) * 8 - 1);
        auto result = square_matrix<T>::identity(this->row_count);
        while (parse > 1) {
            if ((parse & exponent) != 0) {
                result *= (*this);
            }
            result.square();
            parse >>= 1;
        }
        if ((exponent & 1) != 0) {
            result *= (*this);
        }

        return result;
    }

    template <typename T>
    template <typename O>
    square_matrix<T> square_matrix<T>::pow_helper(const square_matrix<O>& mat,
                                                  size_t exponent) {

        size_t parse = (size_t)1 << (sizeof(size_t) * 8 - 2);

        auto converted_mat = square_matrix<T>::force_convert(mat);
        auto result = square_matrix<T>::identity(converted_mat.row_count);
        while (parse > 1) {
            if ((parse & exponent) != 0) {
                result *= converted_mat;
            }
            result.square();
            parse >>= 1;
        }
        if ((exponent & 1) != 0) {
            result *= converted_mat;
        }

        return static_cast<square_matrix<T>>(result);
    }

// Template instantiation for square_matrix<float>
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<int>& mat, size_t exponent);
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<long>& mat, size_t exponent);
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<long long int>& mat,
                                   size_t exponent);
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<float>& mat,
                                   size_t exponent);
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<double>& mat,
                                   size_t exponent);
    template square_matrix<int>
    square_matrix<int>::pow_helper(const square_matrix<long double>& mat,
                                   size_t exponent);

// Template instantiation for square_matrix<double>
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<int>& mat, size_t exponent);
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<long>& mat,
                                    size_t exponent);
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<long long int>& mat,
                                    size_t exponent);
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<float>& mat,
                                    size_t exponent);
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<double>& mat,
                                    size_t exponent);
    template square_matrix<long>
    square_matrix<long>::pow_helper(const square_matrix<long double>& mat,
                                    size_t exponent);

// Template instantiation for square_matrix<long double>
    template square_matrix<long long int>
    square_matrix<long long int>::pow_helper(const square_matrix<int>& mat,
                                             size_t exponent);
    template square_matrix<long long int>
    square_matrix<long long int>::pow_helper(const square_matrix<long>& mat,
                                             size_t exponent);
    template square_matrix<long long int> square_matrix<long long int>::pow_helper(
            const square_matrix<long long int>& mat, size_t exponent);
    template square_matrix<long long int>
    square_matrix<long long int>::pow_helper(const square_matrix<float>& mat,
                                             size_t exponent);
    template square_matrix<long long int>
    square_matrix<long long int>::pow_helper(const square_matrix<double>& mat,
                                             size_t exponent);
    template square_matrix<long long int>
    square_matrix<long long int>::pow_helper(const square_matrix<long double>& mat,
                                             size_t exponent);

// Template instantiation for square_matrix<float>
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<int>& mat,
                                     size_t exponent);
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<long>& mat,
                                     size_t exponent);
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<long long int>& mat,
                                     size_t exponent);
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<float>& mat,
                                     size_t exponent);
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<double>& mat,
                                     size_t exponent);
    template square_matrix<float>
    square_matrix<float>::pow_helper(const square_matrix<long double>& mat,
                                     size_t exponent);

// Template instantiation for square_matrix<double>
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<int>& mat,
                                      size_t exponent);
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<long>& mat,
                                      size_t exponent);
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<long long int>& mat,
                                      size_t exponent);
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<float>& mat,
                                      size_t exponent);
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<double>& mat,
                                      size_t exponent);
    template square_matrix<double>
    square_matrix<double>::pow_helper(const square_matrix<long double>& mat,
                                      size_t exponent);

// Template instantiation for square_matrix<long double>
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<int>& mat,
                                           size_t exponent);
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<long>& mat,
                                           size_t exponent);
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<long long int>& mat,
                                           size_t exponent);
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<float>& mat,
                                           size_t exponent);
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<double>& mat,
                                           size_t exponent);
    template square_matrix<long double>
    square_matrix<long double>::pow_helper(const square_matrix<long double>& mat,
                                           size_t exponent);


    template<typename T>
    [[maybe_unused]]T square_matrix<T>::det() const {
        // Check if the matrix is square
        if (this->row_count != this->col_count) {
            throw std::invalid_argument(
                    "Matrix is not square, and thus not invertible.");
        }

        size_t size = this->row_count;

        square_matrix<T> helper(*this);

        T carry = (T) 1;

        bool sign = true;

        for (size_t i = 0; i < size - 1; ++i) {
            // Find pivot index and value
            size_t pivot_index = i;
            T pivot = helper(i, i);

#pragma omp parallel for reduction(max : pivot) shared(pivot_index)
            for (size_t j = i + 1; j < size; ++j) {
                if (std::abs(helper(j, i)) > std::abs(pivot)) {
#pragma omp critical
                    {
                        pivot = helper(j, i);
                        pivot_index = j;
                    }
                }
            }

            // pivot == 0
            if (pivot == T{}) {
                return T{};
            }

            // Swap rows of helper and result matrices
            std::swap_ranges(helper.mat_array.get() + i * size,
                             helper.mat_array.get() + (i + 1) * size,
                             helper.mat_array.get() + pivot_index * size);

            if ((pivot_index & 1) == (size & 1)) {
                sign = !sign;
            }

            // Eliminate elements below pivot
            for (size_t k = i + 1; k < size; ++k) {
                T factor = helper(k, i);
#pragma omp parallel for
                for (size_t j = i + 1; j < size; ++j) {
                    helper(k, j) *= pivot;
                    helper(k, j) -= factor * helper(i, j);
                    helper(k, j) /= carry;
                }
            }
            carry = pivot;
        }

        return sign ? helper.mat_array[size * size - 1] : -helper.mat_array[size * size - 1];
    }
} // namespace jh
#pragma once

#include <memory>
#include <cassert>
#include <initializer_list>

/**
 * @brief Namespace for JeongHan Bae's 2D matrix library.
 */
namespace jh{

    /**
     * @brief A generic 2D matrix class.
     *
     * @tparam T Type of numeric elements in the matrix.
     */
    template <typename T>
    class matrix {
    protected:
        size_t row_count;               /**< Number of rows in the matrix. */
        size_t col_count;               /**< Number of columns in the matrix. */
        std::unique_ptr<T[]> mat_array; /**< Array storing matrix elements. */

    public:
        /**
         * @brief Constructs an initialized matrix with specified dimensions.
         *
         * @param rows Number of rows.
         * @param cols Number of columns.
         */
        matrix(size_t rows, size_t cols);

        /**
         * @brief Copy constructor.
         *
         * @param other Another matrix of same element type to copy from.
         */
        matrix(const matrix<T>& other);

        /**
         * @brief Constructs a matrix from an initializer list of initializer lists.
         *
         * @param init Initializer list representing matrix elements.
         */
        [[maybe_unused]] matrix(std::initializer_list<std::initializer_list<T>> init);

        /**
         * @brief Performs a forced conversion of matrix type.
         *
         * @tparam O Type of the original matrix.
         * @tparam T Type of the new matrix.
         * @param mat Another matrix to convert from.
         * @return Converted matrix.
         */
        template<typename O>
        [[maybe_unused]] static matrix<T> force_convert(const matrix<O>& mat);

        /** Virtual destructor */
        virtual ~matrix() = default;

        /**
         * @brief Gets the number of rows in the matrix.
         *
         * @return Number of rows.
         */
        [[nodiscard]] size_t get_row_count() const;

        /**
         * @brief Gets the number of columns in the matrix.
         *
         * @return Number of columns.
         */
        [[nodiscard]] size_t get_col_count() const;

        /**
         * @brief Accesses an element in the matrix for modification.
         *
         * @param row Row index.
         * @param col Column index.
         * @return Reference to the accessed element.
         */
        T& operator()(size_t row, size_t col);

        /**
         * @brief Accesses an element in the matrix for reading.
         *
         * @param row Row index.
         * @param col Column index.
         * @return Const reference to the accessed element.
         */
        const T& operator()(size_t row, size_t col) const;

        /**
         * @brief Adds two matrices element-wise.
         *
         * @param other Another matrix to add.
         * @return Sum of this and other.
         * @throw std::invalid_argument When the dimensions of the matrices don't match.
         */
        matrix<T> operator+(const matrix<T>& other) const;

        /**
         * @brief Adds another matrix to this matrix and assigns the result.
         *
         * @param other Another matrix to add.
         * @return Reference to the modified matrix.
         * @throw std::invalid_argument When the dimensions of the matrices don't match.
         */
        matrix<T>& operator+=(const matrix<T>& other);

        /**
         * @brief Subtracts another matrix from this matrix element-wise.
         *
         * @param other Another matrix to subtract.
         * @return Difference of this and other.
         * @throw std::invalid_argument When the dimensions of the matrices don't match.
         */
        matrix<T> operator-(const matrix<T>& other) const;

        /**
         * @brief Subtracts another matrix from this matrix and assigns the result.
         *
         * @param other Another matrix to subtract.
         * @return Reference to the modified matrix.
         * @throw std::invalid_argument When the dimensions of the matrices don't match.
         */
        matrix<T>& operator-=(const matrix<T>& other);


        /**
         * @brief Multiplies two matrices.
         *
         * @param other Another matrix to multiply with.
         * @return Product of this and other.
         * @throw invalid_argument When this->col_count and other.row_count don't match.
         */
        matrix<T> operator*(const matrix<T>& other) const;

        /**
         * @brief Multiplies this matrix with another and assigns the result.
         *
         * @param other Another matrix to multiply with.
         * @return Reference to the product matrix.
         */
        matrix<T>& operator*=(const matrix<T>& other);

        /**
         * @brief Resizes the given matrix to new dimensions.
         *
         * @param mat Matrix to resize.
         * @param new_rows New number of rows.
         * @param new_cols New number of columns.
         * @throw invalid_argument New dimensions must have the same total number of elements.
         */
        [[maybe_unused]] static void resize(matrix<T>& mat, size_t new_rows, size_t new_cols);

        /**
         * @brief Resizes the matrix to new dimensions.
         *
         * @param new_rows New number of rows.
         * @param new_cols New number of columns.
         * @return A resized matrix.
         */
        [[maybe_unused]] matrix<T> resize(size_t new_rows, size_t new_cols) const;

        /**
         * @brief Transposes the given matrix.
         *
         * @param mat Matrix to transpose.
         */
        [[maybe_unused]] static void transpose(matrix<T>& mat);

        /**
         * @brief Transposes the matrix.
         *
         * @return A transposed matrix.
         */
        [[nodiscard]] [[maybe_unused]] matrix<T> transpose() const;

        /**
         * @brief Static cast to a derived matrix type. Only realized for the base class.
         *
         * @tparam Derived Derived matrix type.
         * @param mat Matrix to cast.
         * @return Reference to the casted matrix.
         */
        template<typename Derived>
        [[maybe_unused]] static Derived& static_cast_to(matrix<T>& mat);

        // clang-format off
        /**
         * @brief Assignment operator.
         *
         * @param other Another matrix to assign from.
         * @return Reference to the modified matrix.
         * @assert Matrix dimensions must match for assignment.
         */
        virtual matrix<T>& operator=(const matrix<T>& other); // NOLINT
        // clang-format on
    };

    /**
     * @brief A specialized horizontal matrix.
     *
     * @tparam T Type of elements in the matrix.
     */
    template <typename T>
    class horizontal_matrix : public matrix<T> {
    public:

        /**
         * @brief Constructs a horizontal matrix with specified number of columns.
         *
         * @param cols Number of columns.
         */
        explicit horizontal_matrix(size_t cols);

        /**
         * @brief Assignment operator for horizontal matrices.
         *
         * @param other Another matrix to assign from.
         * @return Reference to the modified matrix.
         * @assert Matrix dimensions must match for assignment.
         */
        horizontal_matrix<T>& operator=(const matrix<T>& other);
    };

    /**
     * @brief A specialized vertical matrix.
     *
     * @tparam T Type of elements in the matrix.
     */
    template <typename T>
    class vertical_matrix : public matrix<T> {
    public:

        /**
         * @brief Constructs a vertical matrix with specified number of rows.
         *
         * @param rows Number of rows.
         */
        explicit vertical_matrix(size_t rows);

        /**
         * @brief Assignment operator for vertical matrices.
         *
         * @param other Another matrix to assign from.
         * @return Reference to the modified matrix.
         * @assert Matrix dimensions must match for assignment.
         */
        vertical_matrix<T>& operator=(const matrix<T>& other);
    };

    /**
     * @brief A specialized square matrix.
     *
     * @tparam T Type of elements in the matrix.
     */
    template <typename T>
    class square_matrix : public matrix<T> {
    public:

        /**
         * @brief Constructs a square matrix with specified size.
         *
         * @param size Size of the matrix (number of rows/columns).
         */
        explicit square_matrix(size_t size);

        /**
         * @brief Generates an identity square matrix of specified size.
         *
         * @param size Size of the matrix (number of rows/columns).
         * @return Identity matrix.
         */
        static square_matrix<T>& identity(size_t size);

        /**
         * @brief Generates a square matrix with diagonal elements set to a specified value.
         *
         * @param value Value to set on the diagonal.
         * @param size Size of the square matrix (number of rows/columns).
         * @return Reference to the generated square matrix with diagonal elements set.
         */
        [[maybe_unused]] static square_matrix<T>& diagonal(T value, size_t size);

        /**
         * @brief Assignment operator for square matrices.
         *
         * @param other Another matrix to assign from.
         * @return Reference to the modified matrix.
         * @assert Matrix dimensions must match for assignment.
         */
        square_matrix<T>& operator=(const matrix<T>& other);

        /**
         * @brief Converts a square matrix of type O to type T by creating a new matrix of the same size.
         *
         * @tparam O Type of elements in the input matrix.
         * @param mat Matrix to convert.
         * @return New matrix of type T containing elements converted from the input matrix.
         */
        template<typename O>
        static square_matrix<T> force_convert(const square_matrix<O>& mat);

        /**
         * @brief Performs matrix inversion.
         *
         * @return Inverse of the matrix.
         */
        square_matrix<T> inverse() const {
            static_assert(
                    std::is_floating_point_v<T>,
                    "Inverse operation is only supported for floating-point types.");
            return this->inverse_helper();
        }

        /**
         * @brief Computes the inverse of a square matrix of type O and returns it as a square matrix of type T.
         *
         * @tparam O Type of elements in the input matrix.
         * @param mat Matrix to compute the inverse for.
         * @return Inverse of the input matrix as a square matrix of type T.
         * @note Only supported for floating-point types.
         */
        template <typename O>
        [[maybe_unused]] static square_matrix<T> inverse(const square_matrix<O>& mat) {
            static_assert(
                    std::is_floating_point_v<T>,
                    "Inverse operation is only supported for floating-point types.");
            return square_matrix<T>::inverse_helper(mat);
        }

        /**
         * @brief Performs matrix division.
         *
         * @param other Another matrix to divide with.
         * @return Result of matrix division.
         */
        square_matrix<T> operator/(const square_matrix<T>& other) const;

        /**
         * @brief Performs matrix division and assigns the result.
         *
         * @param other Another matrix to divide with.
         * @return Reference to the modified matrix.
         */
        square_matrix<T>& operator/=(const square_matrix<T>& other);

        /**
         * @brief Computes matrix raised to a power.
         *
         * @param exponent Exponent to raise the matrix to.
         * @return Result of matrix exponentiation.
         */
        [[maybe_unused]] square_matrix<T> pow(size_t exponent) const;

        /**
         * @brief Computes matrix raised to a power.
         *
         * @param mat Matrix to raise to a power.
         * @param exponent Exponent to raise the matrix to.
         * @return Result of matrix exponentiation.
         */
        template <typename O>
        [[maybe_unused]] static square_matrix<T> pow(const square_matrix<O>& mat,
                                                     long long int exponent) {
            static_assert(std::is_floating_point_v<T> || exponent >= 0,
                          "Power operation of negative exponent is only supported for "
                          "floating-point types.");
            if (exponent >= 0) {
                return square_matrix<T>::pow_helper(mat, exponent);
            }
            return square_matrix<T>::pow_helper(mat, -exponent).inverse();
        }

        [[maybe_unused]] T det() const;

    private:
        // helper methods

        /**
         * @brief let mat *= mat
         */
        void square();

        /**
         * @brief Helper function for matrix inversion.
         *
         * @tparam O Type of elements in the matrix to invert.
         * @param mat Matrix to invert.
         * @return Inverse of the matrix.
         */
        template<typename O>
        static square_matrix<T> inverse_helper(const square_matrix<O>& mat);

        /**
         * @brief Helper function for matrix inversion.
         *
         * @return Inverse of the matrix.
         */
        square_matrix<T> inverse_helper() const;

        /**
         * @brief Helper function for matrix exponentiation.
         *
         * @tparam O Type of elements in the matrix to raise to a power.
         * @param mat Matrix to raise to a power.
         * @param exponent Exponent to raise the matrix to.
         * @return Result of matrix exponentiation.
         */
        template<typename O>
        static square_matrix<T> pow_helper(const square_matrix<O>& mat, size_t exponent);
    };
}//namespace jh
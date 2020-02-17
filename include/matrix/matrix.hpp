#ifndef MATRIX_MATRIX_HPP_INCLUDED
#define MATRIX_MATRIX_HPP_INCLUDED


#include <optional>
#include <vector>
#include <type_traits>


namespace matrix
{


template <typename T>
struct ScalarTraits
{
    static bool isZero(const T& x) noexcept
    {
        return x == 0;
    }
};


template <>
struct ScalarTraits<float>
{
    static bool isZero(float x) noexcept
    {
        return abs(x) < 1e-7;
    }
};


template <>
struct ScalarTraits<double>
{
    static bool isZero(double x) noexcept
    {
        return abs(x) < 1e-8;
    }
};


template <>
struct ScalarTraits<long double>
{
    static bool isZero(long double x) noexcept
    {
        return abs(x) < 1e-9;
    }
};


template <typename T>
class Matrix
{
public:
    /// Constructor. Creates a matrix of size numRows × numCols and fills it with the given value
    Matrix(size_t numRows, size_t numCols, const T& fill = T{}):
        _numRows(numRows),
        _numCols(numCols),
        _data(numRows * numCols, fill)
    { }


    /// Indexing (const; boundary checks are not performed)
    const T& at(size_t row, size_t col) const noexcept
    {
        return _data[row * cols() + col];
    }


    /// Indexing (non-const; boundary checks are not performed)
    T& at(size_t row, size_t col) noexcept
    {
        return _data[row * cols() + col];
    }


    /// Return the number of rows
    size_t rows() const noexcept
    {
        return _numRows;
    }


    /// Return the number of columns
    size_t cols() const noexcept
    {
        return _numCols;
    }


    /// Add a multiple of one row to another (in place; boundary checks are not performed)
    void addMultipleOfRow(size_t destinationIndex, size_t sourceIndex, const T& k) noexcept
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            at(destinationIndex, j) += at(sourceIndex, j) * k;
        }
    }


    /// Multiply a row by a scalar (in place; boundary checks are not performed)
    void multiplyRow(size_t index, const T& k) noexcept
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            at(index, j) *= k;
        }
    }


    /// Swap two rows (in place; boundary checks are not performed)
    void swapRows(size_t index1, size_t index2) noexcept
    {
        for (size_t j = 0; j < cols(); ++j)
        {
            std::swap(at(index1, j), at(index2, j));
        }
    }


    /// Gaussian elimination. Transform the matrix into one of its REFs
    void toRef(size_t startRow = 0, size_t startCol = 0) noexcept
    {
        if (startRow >= rows() || startCol >= cols())
        {
            return;
        }

        if (_isZero(at(startRow, startCol)))
        {
            // The first element in the current row is zero
            std::optional<size_t> pivotIndex = findPivotElement(startCol, startRow);
            if (!pivotIndex.has_value())
            {
                // We have a column of zeros
                toRef(startRow, startCol + 1);
                return;
            }
            swapRows(startRow, *pivotIndex);
        }
        multiplyRow(startRow, T(1) / at(startRow, startCol));
        makeZerosBelow(startRow, startCol);

        // Run recursively
        toRef(startRow + 1, startCol + 1);
    }


    /// Backward Gaussian elimination. The second step of transforming the matrix into RREF
    void toRrefBackward() noexcept
    {
        for (size_t i = 0; i < rows(); ++i)
        {
            // Can't loop over `row` directly since it's unsigned: row >= 0 is always true
            size_t row = rows() - i - 1;
            std::optional<size_t> leadingIndex = findLeadingElement(row);
            if (!leadingIndex.has_value())
            {
                continue;
            }

            makeZerosAbove(row, *leadingIndex);
        }
    }


    /// Transform to RREF
    void toRref() noexcept
    {
        toRef();
        toRrefBackward();
    }


private:
    std::optional<size_t> findPivotElement(size_t col, size_t startRow = 0) const noexcept
    {
        for (size_t row = startRow; row < rows(); ++row)
        {
            if (!_isZero(at(row, col)))
            {
                return {row};
            }
        }

        return std::nullopt;
    }


    std::optional<size_t> findLeadingElement(size_t row, size_t startCol = 0) const noexcept
    {
        for (size_t col = startCol; col < cols(); ++col)
        {
            if (at(row, col) != 0)
            {
                return {col};
            }
        }

        return std::nullopt;
    }


    void makeZerosAbove(size_t row, size_t col) noexcept
    {
        for (size_t rowAbove = 0; rowAbove < row; ++rowAbove)
        {
            addMultipleOfRow(rowAbove, row, -at(rowAbove, col) / at(row, col));
        }
    }


    void makeZerosBelow(size_t row, size_t col) noexcept
    {
        for (size_t rowBelow = row + 1; rowBelow < rows(); ++rowBelow)
        {
            addMultipleOfRow(rowBelow, row, -at(rowBelow, col) / at(row, col));
        }
    }


    bool _isZero(const T& element) const noexcept
    {
        return ScalarTraits<T>::isZero(element);
    }


    size_t _numRows;
    size_t _numCols;
    std::vector<T> _data;
};

}


#endif  // MATRIX_MATRIX_HPP_INCLUDED

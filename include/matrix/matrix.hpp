#ifndef MATRIX_MATRIX_HPP_INCLUDED
#define MATRIX_MATRIX_HPP_INCLUDED


#include <cmath>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <iostream>
#include <unordered_map>
#include <variant>
#include <vector>


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
        return std::abs(x) < 1e-7;
    }
};


template <>
struct ScalarTraits<double>
{
    static bool isZero(double x) noexcept
    {
        return std::abs(x) < 1e-7;
    }
};


template <>
struct ScalarTraits<long double>
{
    static bool isZero(long double x) noexcept
    {
        return std::abs(x) < 1e-7;
    }
};


std::string makeVar(size_t index)
{
    return "y_{" + std::to_string(index) + "}";
}


template <typename T>
class Matrix
{
public:
    using Literal = std::variant<std::string, T>;
    using Atom = std::pair<Literal, T>;
    using Linear = std::vector<Atom>;
    using Solution = std::vector<Linear>;
    using MaybeSolution = std::optional<Solution>;


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


    /// Concatenate two matrices (horizontally)
    Matrix<T> concatenatedWith(const Matrix<T>& other) const
    {
        if (rows() != other.rows())
        {
            throw std::logic_error("Mismatched matrix sizes");
        }

        Matrix<T> result(rows(), cols() + other.cols());
        for (size_t row = 0; row < rows(); ++row)
        {
            for (size_t col = 0; col < cols(); ++col)
            {
                result.at(row, col) = at(row, col);
            }

            for (size_t col = 0; col < other.cols(); ++col)
            {
                result.at(row, cols() + col) = other.at(row, col);
            }
        }
        return result;
    }


    /// Solve multiple systems of linear equations
    std::vector<MaybeSolution> solveSle(size_t numSystems)
    {
        Matrix<T> sle = *this;
        sle.toRref();
        std::vector<std::optional<Solution>> sols;
        for (size_t sys = 0; sys < numSystems; ++sys)
        {
            sols.push_back(([&]() -> MaybeSolution {
                if (!sle.isConsistent(numSystems, sys))
                {
                    return std::nullopt;
                }

                // row -> col
                std::unordered_map<size_t, size_t> leadingIndices;
                for (size_t row = 0; row < rows(); ++row)
                {
                    std::optional<size_t> leadingIndex = sle.findLeadingElement(row);
                    if (leadingIndex.has_value())
                    {
                        leadingIndices.emplace(*leadingIndex, row);
                    }
                }

                Solution sol;
                for (size_t col = 0; col < cols() - numSystems; ++col)
                {
                    if (leadingIndices.count(col) == 0)
                    {
                        sol.push_back(Linear{{makeVar(col), 1}});
                    }
                    else
                    {
                        sol.push_back(Linear());
                        sol.back().emplace_back(
                            sle.at(leadingIndices.at(col), sle.cols() - numSystems + sys),
                            1
                        );
                        for (size_t freeCol = col + 1; freeCol < cols() - numSystems; ++freeCol)
                        {
                            const T& element = sle.at(leadingIndices.at(col), freeCol);
                            if (_isZero(element))
                            {
                                continue;
                            }
                            sol.back().emplace_back(
                                makeVar(freeCol),
                                -element
                            );
                        }
                    }
                }
                return sol;
            })());
        }
        return sols;
    }


    bool isConsistent(size_t numSystems, size_t sys) const noexcept
    {
        for (size_t row = 0; row < rows(); ++row)
        {
            bool isZeroRow = true;
            for (size_t col = 0; col < cols() - numSystems; ++col)
            {
                if (!_isZero(at(row, col)))
                {
                    isZeroRow = false;
                    break;
                }
            }

            if (isZeroRow && !_isZero(at(row, cols() - numSystems + sys)))
            {
                return false;
            }
        }
        return true;
    }


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
            if (!_isZero(at(row, col)))
            {
                return {col};
            }
        }

        return std::nullopt;
    }

private:


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

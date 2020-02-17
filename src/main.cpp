#include <matrix/matrix.hpp>

#include <iostream>


using matrix::Matrix;


template <typename T>
std::istream& operator>>(std::istream& in, Matrix<T>& matrix)
{
    for (size_t row = 0; row < matrix.rows(); ++row)
    {
        for (size_t col = 0; col < matrix.cols(); ++col)
        {
            in >> matrix.at(row, col);
        }
    }
    return in;
}


template <typename T>
std::ostream& operator<<(std::ostream& out, const Matrix<T>& matrix)
{
    for (size_t row = 0; row < matrix.rows(); ++row)
    {
        for (size_t col = 0; col < matrix.cols(); ++col)
        {
            out << matrix.at(row, col) << '\t';
        }
        out << '\n';
    }

    return out;
}


template <typename T>
void printLatex(std::ostream& out, const T& x)
{
    out << x;
}


template <typename T>
void printLatex(std::ostream& out, const Matrix<T>& matrix)
{
    out << "\\begin{bmatrix}\n";
    for (size_t row = 0; row < matrix.rows(); ++row)
    {
        out << "    ";
        for (size_t col = 0; col < matrix.cols(); ++col)
        {
            printLatex(out, matrix.at(row, col));
            out << " & ";
        }
        out << "\\\\\n";
    }
    out << "\\end{bmatrix}\n";
}


int main()
{
    size_t rows, cols;
    std::cin >> rows >> cols;
    Matrix<double> matrix(rows, cols);
    std::cin >> matrix;
    matrix.toRref();
    std::cout << matrix << std::endl;
    //printLatex(std::cout, matrix);
}

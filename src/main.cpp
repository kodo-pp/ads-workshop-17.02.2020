#include <matrix/matrix.hpp>

#include <iostream>
#include <set>


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


template <typename T>
void printSolution(std::ostream& out, const typename Matrix<T>::MaybeSolution& sol)
{
    if (!sol.has_value())
    {
        out << "No solutions\n";
        return;
    }

    size_t n = sol->size();
    std::set<std::string> vars;

    out << "\\begin{cases}\n";
    for (size_t i = 0; i < n; ++i)
    {
        out << "x_{" << std::to_string(i) << "} = ";
        bool first = true;
        for (const typename Matrix<T>::Atom& atom : (*sol)[i])
        {
            if (!first)
            {
                out << " + ";
            }
            first = false;
            auto [literal, coefficient] = atom;
            if (std::abs(coefficient - 1) > 1e-7) {
                out << coefficient << " \\cdot ";
            }
            std::visit([&out, &vars](const auto& x) {
                using ActualType = std::decay_t<decltype(x)>;
                if constexpr(std::is_same_v<ActualType, std::string>)
                {
                    vars.insert(x);
                }
                out << x;
            }, literal);
        }
        out << "\\\\\n";
    }
    out << "\\end{cases}\\quad\n";
    for (const std::string& var : vars)
    {
        out << ", " << var;
    }

    if (!vars.empty())
    {
        out << " \\in \\mathbb{R}";
    }
}


int main()
{
    size_t numEquations, numVars, numSystems;
    std::cin >> numEquations >> numVars >> numSystems;
    Matrix<double> matrix(numEquations, numVars + numSystems);
    std::cin >> matrix;

    using MaybeSolution = typename decltype(matrix)::MaybeSolution;
    std::vector<MaybeSolution> sols = matrix.solveSle(numSystems);
    for (const MaybeSolution& sol : sols)
    {
        printSolution<double>(std::cout, sol);
        std::cout << "\\\\\n";
    }
}

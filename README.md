# SLE solver (ADS Workshop 17.02.2020)
Uses Gaussian elimination to solve multiple systems of linear equations.

## Input format
```
M N K
[matrix M × (N+K)]
```

M = the number of equations, N = the number of variables, K = the number of systems.

Multiple systems share the coefficients but differ in the free terms. Each row of free terms
constitutes a separate SLE.

## Output format
LaTeX, can be safely pasted into Jupyter Notebook


## License
[Apache2](LICENSE)


## Examples

- [input1.txt](input1.txt)
- [input2.txt](input2.txt)
- [input3.txt](input3.txt)
- [input4.txt](input4.txt)

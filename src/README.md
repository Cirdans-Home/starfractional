# Code to generate matrix equivalents of ⋆ convolution

| **Routine**                            | **Description**   |  
| -------------------------------------  | ----------------- |
| `Bd = genBasisMatrix(d,M)`             | Generates the representation in the the basis of orthonormal Legendre polynomials of the Heaviside $\Theta$ function |
| `Firstrow = Toeplitzpart(d,m)`         | Toeplitz part of the $\Theta$ representation                   |
| `[Firstcol,LastRow] = Hankelpart(d,m)` | Hankel part of the $\Theta$ representation                  |
| `[F, k] = genCoeffMatrix_SP_interval(fV,M,I,trunc) ` | Computes the coefficient matrix of $f(t)\Theta(t-s)$ with a smooth function f and the Heaviside theta function Theta, in a basis of orthonormal Legendre polynomials |

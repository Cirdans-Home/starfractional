# Scalar Examples

The `caputo_legendre_scalar_examples.m` script allows to replicate the examples discussed in Section~6.1 of the Paper. In particular, the examples solve the cases:
1. $y^{(\alpha)}(t) = F y(t)$ with $t \in [0,2]$ and $y(0) = y_0 = 1$,
2. $y^{(\alpha)}(t) = t y(t)$ with $t \in [0,2]$ and $y(0) = y_0 = 1$ for two different values ​​of $\alpha$, i.e. $\alpha = \frac{1}{2}$ and $\alpha = {1}{3}$, where the solution can be written using the closed-form $\star$ formalism as a sum of hypergeometric functions.

The script allows you to select the number of Legendre coefficients of the 
solution to retain to control the truncation error.  
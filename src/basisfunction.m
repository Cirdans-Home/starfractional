function A = basisfunction(m,alpha)
%BASISFUNCTION Computes the m x m basis matrix in the shifted and
%normalized Legendre basis. 
%   INPUT: m is the number of basis functions to be used
%          alpha is the order of the Fractional Derivative
%   OUTPUT: A m x m matrix containing the basis expansion.
%The code uses the stable evaluation of the Riemann-Liouville of the
%Legendre polynomials in term of the Gauss-Hypergeometric function. The
%outer quadrature is computed by means of the Gauss-Legendre-Lobatto
%quadrature on 2m nodes.
%

I = [-1,1];
A = zeros(m,m);
[xi,omega] = lobpts(2*m);
fracpolval = legfracintgauss(alpha,m,xi);

for i=1:m
    Li = legpoly(i-1,I,'norm');
    for j=1:m
        A(i,j) = omega*(Li(xi).*fracpolval(j,:).');
    end
end


end


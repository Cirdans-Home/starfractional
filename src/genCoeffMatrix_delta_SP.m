function [F, k] = genCoeffMatrix_delta_SP(fV,M,I)
%GENCOEFFMATRIX Computes the coefficient matrix of f(t)*Theta(t-s), with a
%smooth function f and the Heaviside theta function Theta, in a basis of
%orthonormal Legendre polynomials
%INPUT:
%   fv = cell with the functions f(t) of interest
%   M = size of basis
%OUTPUT:
%   F = cell of the MxM matrices containing the Fourier coefficients of f(t)*Theta(t-s)
%NOTE:
%   Requires chebfun for the computation of the Legendre coefficients of f

% MODIFIED by Pozza 31/5/2024

t0 = I(1);
tfinal = I(2);

nf = length(fV);
c = cell(nf,1);
k = zeros(nf,1);
F = cell(nf,1);

for j=1:nf
    f = fV{j};
    c{j} = cheb2leg(chebcoeffs(chebfun(@(t) f(t0 + (t+1)*(tfinal-t0)/2)*(tfinal-t0)/2)));
    k(j) = length(c{j});
    F{j} = sparse(M,M);
end
for d = 0:max(k)-1
    Bd = genBasisMatrix(d,M-1);
    for j=1:nf
        if d < k(j)
            F{j} = F{j} + c{j}(d+1)*Bd*sqrt(2)/sqrt(2*d+1);
        end
    end
end
end
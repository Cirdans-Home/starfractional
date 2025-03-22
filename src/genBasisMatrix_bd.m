function Bd = genBasisMatrix_bd(d,M)
%GENBASISMATRIX generates the basismatrix of orthonormal Legendre
%polynomials times theta in a basis spanned by orthonormal Legendre
%polynomials
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   M = size of the basis
%OUPUT:
%   Bd = MxM matrix containing the coefficients in the basis

Z = diag(ones(M-1,1),-1)+diag(-ones(M-1,1),1); Z(1,1) = 1; Z(M+1,M) = 1;
fact = sqrt(2*(0:M-1)+1)'./sqrt(2*(0:M-1)+1);% Constant factor
[Firstcol,LastRow] = Hankelpart(d,M+1);
ToeplRow = Toeplitzpart(d,M+1);
H = hankel(Firstcol,LastRow); H = H(1:M,1:M+1);
T = toeplitz(ToeplRow); T = T(1:M,1:M+1);

Bd = sqrt(2*d+1)/sqrt(8)*fact.*((H.*T)*Z);

end
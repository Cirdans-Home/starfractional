function Lfrac = legfracintgauss(alpha,n,x)
%LEGFRACINTGAUSS Computes the Riemann-Liouville Fractional Integral of 
%order alpha of all the shifted and scaled Legendre polynomials up to 
%degree n and evaluates them on x using Gauss hypergeometric function
%formula.

Lfrac = zeros(n+1,length(x));     % Will contain the polynomials on output
sc = ((1+x).^alpha)/gamma(1+alpha);
for i=0:n
    Lfrac(i+1,:) = (-1)^i.*sc.*sqrt((2*i+1)/2)...
        .*hypergeom([i+1,-i],1+alpha,(1+x)/2);
end

end
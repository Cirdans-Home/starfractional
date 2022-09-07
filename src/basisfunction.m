function A = basisfunction(m,alpha)
%BASISFUNCTION Computes the m x m basis matrix in the shifted and
%normalized Legendre basis
%   Detailed explanation goes here

A = zeros(m,m);

for i=1:m
    Li = legpoly(i-1,[0,1],'norm');
    for j=1:m
        Laj = @(x) polint(alpha,j-1,x);
        A(i,j) = integral(@(x) Li(x).*Laj(x),0,1);
    end
end


end


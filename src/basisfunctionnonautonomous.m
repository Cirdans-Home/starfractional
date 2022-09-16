function A = basisfunctionnonautonomous(m,alpha,f,varargin)
%BASISFUNCTION Computes the m x m basis matrix in the shifted and
%normalized Legendre basis
%   Detailed explanation goes here

if nargin == 3
    I = varargin{1};
else
    I = [0,1];
end
A = zeros(m,m);
[xi,omega] = lobpts(2*m);
xi = ((I(2)-I(1))*xi + I(2)+I(1))/2;
omega = omega*(I(2)-I(1))/2;
lobq = @(F) omega*F(xi);

for i=1:m
    Li = legpoly(i-1,I,'norm');
    for j=1:m
        Laj = @(x) (I(2)-I(1)).^(alpha-1)*sqrt(I(2)-I(1))*polint(alpha,j-1,x/(I(2)-I(1)));
        %A(i,j) = integral(@(x) Li(x).*Laj(x),I(1),I(2));
        A(i,j) = lobq(@(x) f(x).*Li(x).*Laj(x));
    end
end


end


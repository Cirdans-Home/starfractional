function p = polint(alpha,j,x)
%POLINT Computes the alphath integral of the jth Legendre polynomial in x
%   alpha order of the fractional integral
%   j order of the scaled and shifted Legendre polynomial
%   x evaluation point

p = @(i) sqrt(2*j+1)*(-1)^(j-i)*factorial(j+i)...
    /(factorial(j-i)*factorial(i)*gamma(alpha+i+1));
rhoj   = p(j);
for i=j:-1:1
    rhojm1 = rhoj*x + p(i-1);
    rhoj = rhojm1;
end
p = rhoj*x^alpha;

end


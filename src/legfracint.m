function Lfrac = legfracint(alpha,n,x)
%LEGFRACINT Computes the Riemann-Liouville Fractional Integral of order
%alpha of all the shifted and scaled Legendre polynomials up to degree n
%and evaluates them on x

Lfrac = zeros(n+1,length(x));     % Will contains the polynomials on output

% Functions containing the recursion coefficients for the Legendre
% polynomials
b = @(j) sqrt(4-j^(-2));
a = @(j) 2*b(j);
d = @(j) ((j-1)/j)*sqrt((2*j+1)/(2*j-3));

for k = 1:length(x)
    Lfracx = zeros(n+1,n+1);  % Lower triangular matrix with the recursion
    xval = x(k);
    for j=0:n % enumerates rows
        for i=j:n % enumerates columns
            if j == 0
                Lfracx(i+1,j+1) = xval^(alpha+i)/(alpha+i);
            elseif j == 1
                Lfracx(i+1,j+1) = (a(1)*xval - b(1))*Lfracx(i,j)...
                    -a(1)*Lfracx(i+1,j);
            else
                Lfracx(i+1,j+1) = (a(j)*xval - b(j))*Lfracx(i,j)...
                    -a(j)*Lfracx(i+1,j) -d(j)*Lfracx(i-1,j-1);
            end
        end
    end
    Lfrac(:,k) = diag(Lfracx)/gamma(alpha); 
end



end
function Firstrow = Toeplitzpart(d,m)
%TOEPLITZPART Generates Toeplitz matrix appearing in formula of Lemma 0.2
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   M = size of the basis
%OUPUT:
%   Firstrow = vector of length M, which is the first row of the Toeplitz
%   matrix
Firstrow = zeros(1,m);
for alpha = 0:m-1
    if mod(alpha+d,2) == 0 && alpha<=d % d+alpha even
        temp = 1;
        count = 1;
        for j = 1:d+alpha
            if j<=(d-alpha)/2                
                temp = temp*(((d+alpha)/2+j)/(2*j))^2;
                count = count + 1;
            end
            if count<=d
                temp = temp/2^2;
                count = count + 1;
            end
            if j<=alpha
                jj = alpha-j+1; % reversing order for stability
                temp = temp*(d+jj)/(d-alpha+jj);
            end
        end
        Firstrow(alpha+1) = temp*2;
    end
end
end
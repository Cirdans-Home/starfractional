function [Firstcol,LastRow] = Hankelpart(d,m)
%HANKELPART Generates Hankel matrix appearing in formula of Lemma 0.2
%INPUT:
%   d = degree of the orthonormal Legendre polynomial to be expanded
%   M = size of the basis
%OUPUT:
%   Firstcol = vector of length M, which is the first column of the Hankel
%   matrix
%   LastRow = vector of length M, which is the last row of the Hankel matrix
Firstcol = zeros(m,1);
LastRow = zeros(1,m);

for sumCol1 = 0:m-1
    if sumCol1>=d && mod(d+sumCol1,2)==0
        temp = 1/(d+sumCol1+1);
        for j=1:d
            temp = temp*(-d+sumCol1+2*j)/(-d+sumCol1+2*j-1);
        end
        Firstcol(sumCol1+1) = temp;
    end
end

for sumRowm = m-1:2*m-2
    if mod(d+sumRowm,2)==0
        temp = 1/(d+sumRowm+1);
        for j=1:d
            temp = temp*(-d+sumRowm+2*j)/(-d+sumRowm+2*j-1);
        end
        LastRow(sumRowm-m+2) = temp;
    end
end
end

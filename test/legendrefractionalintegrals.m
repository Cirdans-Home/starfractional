%% Tests for Riemann-Liouville integrals of Legendre polynomials.

clear; clc; close all;
addpath("../src/")

x = linspace(0,1,10);
n = 14;
alpha = 0.5;

% This routine is unstable for Polynomials of high-degree
values = zeros(n,length(x));
for i=1:n+1
    for j=1:length(x)
        values(i,j) = polint(alpha,i-1,x(j));
    end
end

% This routine uses the three-term recurrences for the integrals
Lfrac = legfracint(alpha,n,x);

disp(values-Lfrac)
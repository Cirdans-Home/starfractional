%% Check of the Gegenbauer-Integral formula
% We check the formula obtained by using Gegenbauer polynomials agains the
% adaptive quadrature rule by Matlab. The absolute error difference between
% the two is of the order of the tolerance used for the adaptive
% quadrature.
% 
clear; clc; close all;

x = linspace(-1,1,100);
n = 30;
alpha = 0.5;
IalphaL = legfracintgauss(alpha,n,x);

% Brute-force integral
IalphaLt = zeros(n+1,length(x));
for i=0:n
    P = legpoly(i,[-1,1],"norm");
    for j=2:length(x)
        f = @(t) P(t).*(x(j)-t).^(alpha-1);
        IalphaLt(i+1,j) = integral(f,-1,x(j),"AbsTol",1e-12,"RelTol",1e-9)...
            /gamma(alpha);
    end
end

%% Plot
figure(1)
for i=0:n
   yyaxis left
   plot(x,IalphaL(i+1,:),'b-',...
       x,IalphaLt(i+1,:),'r--','LineWidth',2);
   yyaxis right
   semilogy(x,abs(IalphaL(i+1,:)-IalphaLt(i+1,:)),'k-');
   pause(1/25)
end

%% 3D Plot
figure(2)
[C,h] = contour(x,0:30,log10(abs(IalphaL-IalphaLt)),[-13,-14,-15],'k');
clabel(C,h)
xlabel('x')
ylabel('n')
zlabel('Absolute error')
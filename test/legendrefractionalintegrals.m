%% Tests for Riemann-Liouville integrals of Legendre polynomials.

clear; clc; close all;
addpath("../src/")

x = linspace(0,1,100);
n = 30;
alpha = 0.5;
I = [0,1];

%% This routine is unstable for Polynomials of high-degree
tic;
values = zeros(n,length(x));
for i=1:n+1
    for j=1:length(x)
        values(i,j) = polint(alpha,i-1,x(j));
    end
end
thorner = toc;

figure(1)
IalphaL = @(i) (I(2)-I(1)).^(alpha-1)*sqrt(I(2)-I(1))*...
    polint(alpha,i,x/(I(2)-I(1)));
L = legpoly(10,I,'norm');
plot(x,L(x),'r-',x,IalphaL(10),'k:','LineWidth',2);
legend("Legendre","Horner");
title(sprintf('degree = %d',10))
figure(2)
L = legpoly(23,I,'norm');
plot(x,L(x),'r-',x,IalphaL(23),'k:','LineWidth',2);
legend("Legendre","Horner");
title(sprintf('degree = %d',23))

% Horners coefficients
j = 23;
p = @(i) sqrt(2*j+1)*(-1)^(j-i)*factorial(j+i)...
    /(factorial(j-i)*factorial(i)*gamma(alpha+i+1));
pval = arrayfun(p,0:j-1);
figure(3)
ind = 0:j-1;
semilogy(ind(pval > 0),pval(pval > 0),'r+',...
    ind(pval < 0),-pval(pval < 0),'bx')
legend({'Positive','Negative'},'Location','southeast')
axis tight

% Higham Bound
tic;
values = zeros(n,length(x));
bounds = zeros(n,length(x));
for i=1:n+1
    for j=1:length(x)
        [values(i,j),bounds(i,j)] = polintwithbound(alpha,i-1,x(j));
    end
end
thorner = toc;

figure(4)
[ORD,XVAL] = meshgrid(x,0:30);
mesh(ORD,XVAL,log10(bounds));
ylabel('Polynomial Degree','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
zlabel('Running error bound $\log_{10}(\cdot)$','Interpreter','latex')

%% This routine uses the three-term recurrences for the integrals
Lfrac = legfracint(alpha,n,x);

f = @() legfracint(alpha,n,x);
trecurrence = timeit(f,1);

fprintf('Time Horner %e (s)\n',thorner);
fprintf('Time Three-terms recurrence %e (s)\n',trecurrence);
fprintf('||Horner - Three-terms|| = %e\n',norm(values-Lfrac));

%% Check of the obtained functions
I = [0,1];
[xi,~] = lobpts(6*n);
x = ((I(2)-I(1))*xi + I(2)+I(1))/2;
Lfrac = legfracint(alpha,n,x);

for i=0:30
    L = legpoly(i,I,'norm');
    IalphaL = (I(2)-I(1)).^(alpha-1)*sqrt(I(2)-I(1))*...
        polint(alpha,i,x/(I(2)-I(1)));
    figure(2)
    subplot(2,2,1)
    plot(x,L(x),'r-',...
        x,Lfrac(i+1,:),'b--','LineWidth',2);
    legend("Legendre","Three-Terms Recurrence");
    subplot(2,2,2)
    plot(x,L(x),'r-',...
        x,IalphaL,'k:','LineWidth',2);
    legend("Legendre","Horner");
    title(sprintf('degree = %d',i))
    subplot(2,2,[3,4])
    semilogy(x,abs(IalphaL.'-Lfrac(i+1,:)),'k--')
    ylabel("Absolute Error")
    pause(1/30)
end

%% USING GAUSS-HYPERGEOMETRIC FUNCTION
% This is the most stable

x = linspace(-1,1,100);
n = 30;
alpha = 0.5;
IalphaL = legfracintgauss(alpha,n,x);
I = [-1,1];
Lfrac  = @(i) (I(2)-I(1)).^(alpha-1)*sqrt(I(2)-I(1))*polint(alpha,i,0.5*x+0.5);
for i=0:30
    L = legpoly(i,I,"norm");
    figure(3)
    subplot(1,2,1)
    plot(x,L(x),'r-',...
        x,IalphaL(i+1,:),'b--',...
        x,Lfrac(i),'k-','LineWidth',2);
    legend("Legendre","Gauss-Hypergeometric","Horner");
    title(sprintf('degree = %d',i))
    subplot(1,2,2)
    plot(x,abs(IalphaL(i+1,:)-Lfrac(i)),'k--','LineWidth',2);
    pause()
end

%% 
figure(6)
[nval,xval] = meshgrid(0:30,x);
mesh(nval,xval,IalphaL.');
xlabel('$n$','Interpreter','latex')
ylabel('$x$','Interpreter','latex')
zlabel('${}^{RL}I_{[0,x]}^{0.5}L_n(x)$','Interpreter','latex')




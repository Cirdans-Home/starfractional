%% Tests for Riemann-Liouville integrals of Legendre polynomials.

clear; clc; close all;
addpath("../src/")

x = linspace(0,1,10);
n = 30;
alpha = 0.5;

%% This routine is unstable for Polynomials of high-degree
tic;
values = zeros(n,length(x));
for i=1:n+1
    for j=1:length(x)
        values(i,j) = polint(alpha,i-1,x(j));
    end
end
thorner = toc;

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
for i=0:n
    L = legpoly(i,I,'norm');
    IalphaL = (I(2)-I(1)).^(alpha-1)*sqrt(I(2)-I(1))*...
        polint(alpha,i,x/(I(2)-I(1)));
    figure(1)
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
    pause()
end
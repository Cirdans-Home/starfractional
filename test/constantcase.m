%% CONSTANT CASE
% Solution of the constant coefficient FODE
% y^{(\alpha)}(t) = \lambda y(t), \lambda < 0, alpha \in (0,1), t \in [0,1]
% y^{(0)} = 1
% With true solution given by the Mittag-Leffler function
% y(t) = E_{\alpha,1}(lambda t^\alpha) 
clear; clc; close all;
addpath('../src');  % Add to the path the gallery function

%% Definition of the system
lambda = -1;
y0 = 1;
alpha = 0.5;

%% Discretization parameters
M = 5; % Number of basis polynomials

%% Builiding the discretization
I = eye(M);
FW = lambda*basisfunction(M,alpha);
F = lambda*basisfunction(M,1);
H = basisfunction(M,alpha);
H = spdiags(spdiags(H,-10:10),-10:10,M,M);
Lpol = legpoly(0:M-1,[0,1],'norm');
rhs = F*Lpol(0).';

A = I - FW;
x = A\rhs;
x = H*x;

figure(1)
semilogy(abs(x),'*');
xlabel('Coefficient')
ylabel('|x|')


%% Reassembling the solution
trunc = 14;
Lpol = legpoly(0:trunc-1,[0,1],'norm');
c = x(1:trunc);
y = y0 + Lpol*c;

t = linspace(0,1,100);
figure(2)
subplot(1,2,1);
plot(t,y(t),'r',...
    t,ml(lambda*t.^alpha,alpha),'k--','LineWidth',2)
xlabel('t')
ylabel('y(t)')
subplot(1,2,2);
semilogy(t,abs(y(t)-ml(lambda*t.^alpha,alpha)),'r-',...
    t,abs(y(t)-ml(lambda*t.^alpha,alpha))./abs(ml(lambda*t.^alpha,alpha)),...
    'k--','LineWidth',2);
legend('Absolute Error','Relative Error')



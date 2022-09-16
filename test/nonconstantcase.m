%% NON-CONSTANT CASE
% Solution of the constant coefficient FODE
% y^{(\alpha)}(t) = \lambda(t) y(t), \lambda < 0, alpha \in (0,1), t \in [0,1]
% y^{(0)} = 1
clear; clc; close all;
addpath('../src');  % Add to the path the gallery function

%% Definition of the system
lambda = @(t) sin(2*pi*t);
y0 = 1;
alpha = 0.3;
Int = [0,1];

fdefun = @(t,y) lambda(t).*y;
Jfdefun = @(t,y) lambda(t);
t0 = Int(1);
tfinal = Int(2);
h = 0.01;
[TT,YY] = flmm2(alpha,fdefun,Jfdefun,t0,tfinal,y0,h);

%% Discretization parameters
M = 15; % Number of basis polynomials

%% Builiding the discretization
I = eye(M);
FW = basisfunctionnonautonomous(M,alpha,lambda,Int);
F = basisfunctionnonautonomous(M,1,lambda,Int);
H = basisfunction(M,alpha,Int);
%H = spdiags(spdiags(H,-10:10),-10:10,M,M);
Lpol = legpoly(0:M-1,Int,'norm');
rhs = F*Lpol(Int(1)).';

A = I - FW;
x = A\rhs;
x = H*x;

figure(1)
semilogy(abs(x),'*');
xlabel('Coefficient')
ylabel('|x|')


%% Reassembling the solution
trunc = 15;
Lpol = legpoly(0:trunc-1,Int,'norm');
c = x(1:trunc);
y = y0 + Lpol*c;

t = linspace(Int(1),Int(2),100);
figure(2)
subplot(1,2,1);
plot(TT,y(TT),'r',...
    TT,YY,'k--','LineWidth',2)
xlabel('t')
ylabel('y(t)')
subplot(1,2,2);
semilogy(TT,abs(y(TT)-YY),'r-',...
    TT,abs(y(TT)-YY)./abs(YY),...
    'k--','LineWidth',2);
legend('Absolute Error','Relative Error')



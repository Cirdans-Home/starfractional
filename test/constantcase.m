%% CONSTANT CASE
% Solution of the constant coefficient FODE
% y^{(\alpha)}(t) = \lambda y(t), \lambda < 0, alpha \in (-1,1), t \in [0,1]
% y^{(0)} = 1
% With true solution given by the Mittag-Leffler function
% y(t) = E_{\alpha,1}(lambda (t-1)^\alpha) 
clear; clc; close all;
addpath('../src');  % Add to the path the gallery function

%% Input Parser
opts.Interpreter = 'tex';
userinput = inputdlg({'\alpha','Number of basis elements (M)'},...
    'Constant case inputs',[1,50;1,50],{'0.5','40'},opts);

%% Definition of the system
lambda = -1;
y0 = 1;
alpha = str2double(userinput{1});
Int = [-1,1];

%% Discretization parameters
M = str2double(userinput{2}); % Number of basis polynomials

%% Builiding the discretization
I = eye(M);
FW = lambda*basisfunction(M,alpha,Int);
F = lambda*basisfunction(M,1,Int);
H = basisfunction(M,alpha,Int);
Lpol = legpoly(0:M-1,Int,'norm');
rhs = F*Lpol(Int(1)).';

A = I - FW;
x = A\rhs;
x = H*x;

figure(1)
semilogy(abs(x),'*');
xlabel('Coefficient')
ylabel('|x|')

%% How many coefficients do we need to retain?
userinput = inputdlg("Number of coefficients to keep:",...
    "How many coefficients do we keep?",[1,50],{num2str(M)},opts);

%% Reassembling the solution
trunc = min(str2double(userinput{1}),M);
Lpol = legpoly(0:trunc-1,Int,'norm');
c = x(1:trunc);
y = y0 + Lpol*c;

t = linspace(Int(1),Int(2),100);
figure(2)
subplot(1,2,1);
plot(t,y(t),'r',...
    t,ml(lambda*(t-Int(1)).^alpha,alpha)*y0,'k--','LineWidth',2)
xlabel('t')
ylabel('y(t)')
subplot(1,2,2);
semilogy(t,abs(y(t)-ml(lambda*(t-Int(1)).^alpha,alpha)*y0),'r-',...
    t,abs(y(t)-ml(lambda*(t-Int(1)).^alpha,alpha)*y0)./abs(ml(lambda*t.^alpha,alpha)*y0),...
    'k--','LineWidth',2);
legend('Absolute Error','Relative Error')



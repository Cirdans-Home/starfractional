%% Scalar test problems

clear; clc; close all;
%% Examples
fprintf("Examples:\n")
fprintf("1) Mittag-Leffler function\n")
fprintf("2) Non autonomous: explicit solution alpha = 1/2\n")
fprintf("3) Non autonomous: explicit solution alpha = 1/3\n")
fprintf("4) Non autonomous: explicit solution alpha = 1\n")
n_example = input("Select the example: ");
% al -> alpha, [l1, l2] -> time interval, u0 -> initial value, fV{1} -> coeff
% sol -> explicit solution, m -> number of Legendre coeff
switch  n_example
    case 1
        % Example 1: Mittag-Leffler function
        al = 0.7; % free choice for 0<al<=1
        u0 = 1;
        sol = @(t) mlf(al,1,-(t).^al,10);
        fV{1} = @(t) -1 + 0*t;
        l1=0;
        l2=2;
        m = 200;
    case 2
        % Example 2: non-autonomous example with explicit solution
        % alpha = 1/2
        fV{1} = @(t) t;
        l1=0;
        l2=2;
        al = 1/2;
        u0 = 1;
        % True solution for alpha = 1/2
        sol = @(t) u0*(4*t.^(3/2)/(3*sqrt(pi)).*...
            hypergeom([1, 4/3],[7/6, 3/2],t.^3/3)...
            + hypergeom(5/6,2/3,t.^3/3)); 
        m = 100;
    case 3
        % Example 2: non-autonomous example with explicit solution
        % alpha = 1/3
        fV{1} = @(t) t;
        l1=0;
        l2=2;
        al = 1/3;
        u0 = 1;
        % True solution for al = 1/3
        sol = @(t) u0*(hypergeom([7/12,11/12],[1/2,3/4],t.^4/4) ...
              + 63*sqrt(3)*t.^(8/3)*...
              gamma(1/3)/(160*pi).*...
              hypergeom([1,5/4,19/12],[7/6,17/12,5/3],t.^4/4) ...
              + 9*sqrt(3)*t.^(4/3)*gamma(2/3)/(8*pi).*...
              hypergeom([11/12,1,5/4],[5/6,13/12,4/3],t.^4/4)); 
        m = 100;
    case 4
        % Example 2: non-autonomous example with explicit solution
        % alpha = 1 -> ODE
        fV{1} = @(t) t;
        l1=0;
        l2=2;
        al = 1;
        u0 = 1;
        % True solution for alpha = 1
        sol = @(t) u0*exp(t.^2/2); 
        m = 50;
end

%% Star solution
trunc = 1;
[Fdeltaleg,~] = genCoeffMatrix_delta_SP(fV,m,[l1,l2]);
[Flegvec, ~] = genCoeffMatrix_SP_interval(fV,m,[l1,l2],trunc);
H = Flegvec{2}*(l2-l1)/2;
F = Fdeltaleg{1}*2/(l2-l1);

phi = zeros(m,1);
for i=1:m
    phi(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);
end

% The slow decay of this matrix determines the accuracy of the 
% approximation:
Hal = H^al;    
FHal = F*Hal;

% Hence we use Hal just once to introduce less less truncation errors:
I = eye(m,m);
R = (I - Hal*F);
x = R\(sqrt(2)*eye(m,1))*2/(l2-l1);
uleg = x*u0;

%% regularization legendre
figure('Position',[172 472 1560 419])
subplot(1,3,1)
semilogy(abs(uleg),'k+')
title('Computed Legendre coeff of the solution')
k = input('Insert the number of Legendre coefficients you want to keep: ');
uleg = uleg(1:k);


%% Solution in terms of Chebyshev expansion
% From (normalized) Legendre to Chebyshev coefficients
cheb_c = leg2cheb(full(uleg),'norm');
u_approx = chebfun(cheb_c, 'coeffs');


%% sol
tleg = linspace(-1,1,m);  % Discr. interval for Legendre comparison [-1,1]
xleg = (tleg+1)*(l2-l1)/2 + l1;
sol_leg = sol(xleg);
subplot(1,3,2)
plot(xleg, real(sol_leg),'rs',xleg, real(u_approx(tleg)), 'k+') 
title('solution')
legend('analytic','star')
subplot(1,3,3)
semilogy(xleg,abs((u_approx(tleg)-sol_leg)./sol_leg),'k+')
title('relative error')

fprintf("The maximum absolute relative error %e\n",...
    max(abs((u_approx(tleg)-sol_leg)./sol_leg)));


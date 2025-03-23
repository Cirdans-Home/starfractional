%% Fractional Schroedinger equation
% Discretization of the time-dependent Schroedinger equation in 2D using
% the MATLAB PDE Toolbox for the spatial discretization and FODE methods
% for the time integration.

% The time-dependent Schroedinger equation is given by
% i * psi_t = -1/2 * (psi_xx + psi_yy) + V * psi
% where psi is the wave function, V is the potential, and the subscripts
% denote partial derivatives.

% The problem is defined on a two-dimensional square domain
% Ω = [−2,2] × [−2,2].
% The potential V (x) is set to zero within the interior of the
% subdomain Ωp = [−1,1]×[−1,1] and is constant outside, specifically
% V(x) = 0 for x ∈ int(Ωp) and V(x) = 10 for x ∈ cl(Ω/Ωp).
% The initial wavefunction at time t = 0 is chosen to be of Gaussian
% type, given by ψ(0,x) = e^{−|x|^2/2}. For simplicity, all constants
% are normalized to 1.

% 18/02/2025
% The solution is also computed using the star-approach

clear; clc; close all;
addpath('src/')

% We discretize the spatial domain Ω using the MATLAB PDE Toolbox.

%% Plot?
fprintf("Plot the Solution?\n")
fprintf("1 - no plot\n")
fprintf("0 - plot\n")
noplot = input("Plot = ");

%% Create a PDE model
model = createpde();

%% Define the geometry
gd = [3,4,-2,2,2,-2,-2,-2,2,2]';
ns = char('R1');
ns = ns';
sf = 'R1';
g = decsg(gd,sf,ns);

%% Include the geometry in the model
geometryFromEdges(model,g);

%% Plot the geometry to verify
f1 = figure("Position",[675 598 861 388]);
subplot(1,2,1)
pdegplot(model,'EdgeLabels','on');
title('Geometry with Edge Labels');
axis tight

%% Generate the mesh
generateMesh(model,'Hmax',0.3,'GeometricOrder','quadratic');

%% Plot the mesh to verify
figure(f1)
subplot(1,2,2)
pdemesh(model);
title('Mesh');
axis tight

%% Define the coefficients of the PDE
alpha = -1;
while alpha < 0 || alpha > 1 || ~isnumeric(alpha)
    fprintf("Insert fractional order in (0,1)\n")
    alpha = input("alpha = ");
end
specifyCoefficients(model,'m',0,...
    'd',1i^alpha,...
    'c',1/2,...
    'a',@V,...
    'f',0);

%% Define the initial condition
setInitialConditions(model,@u0);

%% Define the boundary conditions
applyBoundaryCondition(model,'dirichlet',...
    'Edge',1:model.Geometry.NumEdges,'u',0);

%% Now we want to perform the time integration using a fractional integrator

% Write the time-dependent Schroedinger equation as a system of first-order
% ODEs. We firs extract the initial condition as a vector

xnode = model.Mesh.Nodes(1,:)';
ynode = model.Mesh.Nodes(2,:)';
u0vec = u0(struct('x',xnode,...
    'y',ynode));

% Define the time span
t0=0;
tf=5;  % 5
tspan = [0,tf];

fprintf("How to impose boundary conditions:\n")
fprintf("1) stiff-spring\n")
fprintf("2) nullspace\n")
bcmethod = input("Boundary choice: ");

% Produce all the necessary FEM matrices
if bcmethod == 1
    FEM = assembleFEMatrices(model,"stiff-spring");
    odefun = @(t,u) FEM.M\(-FEM.Ks*u -FEM.Fs);
    A = -FEM.M\FEM.Ks;
    B = -FEM.M\FEM.Fs;
elseif bcmethod == 2

    FEM = assembleFEMatrices(model,"nullspace");
    odefun = @(t,u) FEM.B*(FEM.M\(- FEM.Kc*(FEM.B'*u) - FEM.Fc)) + FEM.ud;
    A = -FEM.B*(FEM.M\(FEM.Kc*FEM.B'));
    B = -FEM.B*(FEM.M\FEM.Fc) + FEM.ud;
else
    error("Unknown BC method.")
end



%% Solve the FDE by passing it to flmm2
tic
h = 1e-2;
fprintf("Select method: \n")
fprintf("1) Trapezoidal\n")
fprintf("2) Newton-Gregory\n")
fprintf("3) F-BDF2\n")
fprintf("4) F-Predictor/Corrector\n")
method = input("method = ");
if method == 1 || method == 2 || method == 3
    jfun = @(t,u) A;
    [tlist,u] = flmm2(alpha, @(t,u) odefun(t,u),...
        @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
        [],method);
elseif method == 4
    [tlist,u] = fde12(alpha, @(t,u) odefun(t,u),tspan(1),tspan(2),...
        u0vec,h,[],10);
else
    error("Unknown method!")
end
time_method = toc;

%% Visualization
if ~noplot
    h3 = figure("Position",[219 453 1515 346]);
    pause()
    for i = 1:numel(tlist)
        subplot(1,3,1);
        pdeplot(model,'XYData',real(u(:,i)));
        title(['(Real) Solution at t = ',num2str(tlist(i))]);
        subplot(1,3,2);
        pdeplot(model,'XYData',imag(u(:,i)));
        title(['(Imaginary) Solution at t = ',num2str(tlist(i))]);
        subplot(1,3,3);
        pdeplot(model,'XYData',abs(u(:,i)).^2);
        title(['Probability Density at t = ',num2str(tlist(i))]);
        drawnow;
        pause(1/100);
    end

end
%% Solve the ODE by star-approach
% The solution requires a precomputed basis for the Legendre expansion
% this can be achieved by the gen_schur_hside.m function in the src
% subfolder.
tic
m = 1000; % possible values: 100, 500, 1000
max_it = 80;
int = [t0,tf];
[u_leg,X] = star_frac_arnoldi(A,u0vec,alpha,int,max_it,m);
time_star = toc;

figure(4)
clf,
semilogy(abs(u_leg),'bo')
title('Abs value of the computed legendre coeff of the sol')

k = input('Insert the number of Legendre coefficients you want to keep: ');
u_leg = u_leg(1:k,:);


figure(5)
clf,
semilogy(svd(X),'kx')
title('Singular values of the solution of the star matrix eq.')

%% Chebfun representation of the solution
% From (normalized) Legendre to Chebyshev coefficients
cheb_c = leg2cheb(full(u_leg),'norm');
u_star = chebfun(cheb_c, 'coeffs');

fprintf('Error at t=%1.1d: %1.3e \n', tf, norm(u_star(1)-u(:,end).'))

%% Visualization
if ~noplot
    tlist_star = -1:0.01:1;
    Uplot = zeros(size(A,1),length(tlist_star));
    h3 = figure("Position",[219 453 1515 346]);
    pause()
    for i = 1:numel(tlist_star)
        Uplot(:,i) = u_star(tlist_star(i),:).';
    end
    tlist_star = (tlist_star +1)*(tf-t0)/2 + t0;
    for i = 1:numel(tlist_star)
        subplot(1,3,1);
        pdeplot(model,'XYData',real(Uplot(:,i)));
        title(['(Real) Solution at t = ',num2str(tlist_star(i))]);
        subplot(1,3,2);
        pdeplot(model,'XYData',imag(Uplot(:,i)));
        title(['(Imaginary) Solution at t = ',num2str(tlist_star(i))]);
        subplot(1,3,3);
        pdeplot(model,'XYData',abs(Uplot(:,i)).^2);
        title(['Probability Density at t = ',num2str(tlist_star(i))]);
        drawnow;
        % pause(1/1000);
    end

    set(h3,'color','white');
end

%% Auxiliary functions
% The following function defines the potential V(x) as described above.

function V = V(x,~)
%V - Potential function for the time-dependent Schroedinger equation, since
%the potential is piecewise constant, we can define it as a function of
% x.x and x.y,where x.x and x.y are the x and y coordinates of the point
% x and are passed to the function as vectors. We do it in a loop to check
% if the point is inside the region or not.
n1 = 1;
nr = numel(x.x);
V = zeros(n1,nr);
for i = 1:nr
    if x.x(i) >= -1 && x.x(i) <= 1 && x.y(i) >= -1 && x.y(i) <= 1
        V(1,i) = 0;
    else
        V(1,i) = 10.0;
    end
end
end
% The following function defines the initial condition for the wave
% function as described above.

function u0 = u0(x)
%u0 - Initial condition for the time-dependent Schroedinger equation
u0 = exp(-((x.x).^2 + (x.y).^2)/2);
end


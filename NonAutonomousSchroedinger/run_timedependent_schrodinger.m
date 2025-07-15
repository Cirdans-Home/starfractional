clear; clc; close all;

addpath('src')
addpath('chebfun');
fid = fopen('output_table.txt','a+');
fprintf(fid,"Date: %s\n\n",string(datetime));
doflmm = true;

% Parameters
for Hmax = [0.3,0.1]
    alpha = 0.3;
    h = 1e-5;
    t0=0;
    tf=1;
    tspan = [0,tf];

    %% Create the PDE model
    model = createpde();
    %% Define the geometry
    gd = [3,4,-2,2,2,-2,-2,-2,2,2]';
    ns = char('R1');
    ns = ns';
    sf = 'R1';
    g = decsg(gd,sf,ns);
    %% Include the geometry in the model
    geometryFromEdges(model,g);
    generateMesh(model,'Hmax',Hmax);
    %% Define the coefficients
    ep= 0.1;
    w = 5;
    ft = @(t) 0.5*(1+ep*sin(w*pi^2*t));
    %% Specify PDE conditions
    specifyCoefficients(model,'m',0,...
        'd',1i^alpha,...
        'c',1/2,...
        'a',@(x,t) V(x,t,ft),...
        'f',0);
    %% Set the initial condition
    setInitialConditions(model,@u0);
    %% Define the boundary conditions
    applyBoundaryCondition(model,'dirichlet',...
        'Edge',1:model.Geometry.NumEdges,'u',0);
    %% Accurate test with fractional integrator
    xnode = model.Mesh.Nodes(1,:)';
    ynode = model.Mesh.Nodes(2,:)';
    u0vec = u0(struct('x',xnode,...
        'y',ynode));
    % Define the time span
    state.time = 0;
    % Assemble FEM matrices
    str = "nullspace";
    state.time = 0;
    FEM = assembleFEMatrices(model,"nullspace",state);
    K0 = FEM.Kc;
    M0 = FEM.M;
    B0 = FEM.B;
    state.time = 1;
    FEM = assembleFEMatrices(model,"nullspace",state);
    K1 = FEM.Kc;
    M1 = FEM.M;
    odefun = @(t,u) odefunt(t,u,FEM,K0,K1-K0,ft,"nullspace");

    %% Solve with FLMM2
    filename = sprintf("fracsolve_HMAX%1.3f_h%1.2e.mat",Hmax,h);
    if ~exist(filename,"file")
        jfun = @(t,u) jfunt(t,FEM,K0,K1-K0,ft,str);
        [tlist,u] = flmm2(alpha, @(t,u) odefun(t,u),...
            @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
            []);
        savefile = sprintf("src/%s",filename);
        save(savefile,"tlist","u",'-v7.3')
    else
        load(filename)
    end

    %% Solve the ODE by star-approach
    for m = [3000,3500] % [100,250,500,750,1000,1500,2000,2500]
	filename = sprintf("src/schur_hside_%d.mat",m);
	if ~exist(filename,"file")
    		warning("Asked for a non existing basis size: precomputing it");
    		gen_schur_hside(m,'schur_hside');
	end
    	tic
        if abs(ft(1)-ft(0)) > 1e-13
            fD{1} = @(t) (ft(t)-ft(0))/(ft(1)-ft(0));  % scalar function of the time-dependent part
        else
            fD{1} = @(t) ft(0) + 0*t;
            warning('division by zero in f(t) interp')
        end        

        trunc = 0;  % No automatic truncation of the last part of the coeff matrices
        [FlegD, ~] = genCoeffMatrix_delta_SP(fD,m,[t0,tf]); % It generates the f(t)delta(t-s) coeff matrix
        F2D = sparse(FlegD{1})*2/(tf-t0); % coeff matrix of fF{1}(t)delta(t-s)
    
        % Load from file
        name = ['schur_hside_',num2str(m),'.mat'];
        load(name)
        Usch = schur_hside.U;
        Tsch = schur_hside.T*(tf-t0)/2;
        phisch = schur_hside.Uphi;
        Hleg = schur_hside.Hleg*(tf-t0)/2;


        if verLessThan('matlab','24.2')
            Tal = expm(alpha*logm(Tsch));
            T1mal = expm((1-alpha)*logm(Tsch));
        else
            Tal = Tsch^alpha;
            T1mal = Tsch^(1-alpha);
        end
        Hal = Usch*Tal*Usch';
        H1mal = Usch*T1mal*Usch';
        % AA{1} = full(-(B0*(M0\(K0*B0'))).'); AA{2} = full(-(B0*(M0\((K1-K0)*B0'))).');
        FF{1} = Hal; FF{2} = F2D*Hal;  % coefficient matrix in the resolvent

        % Legendre poly vector at t=t0
        phi = zeros(m,1);
        for i=1:m
            phi(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);
        end
        H1aphi = H1mal*phi;

        d= length(u0vec);

        % Solution of the system by iterative method
        maxit = 20;     % max n. it. (outer iterations)
        tol = 1e-7;     % tol estimated residual stop crit. (outer iterations)
        it_arn = ceil(2*log(d)); % it. Arnoldi (inner iterations)
        trunc = 1e-9;  % svd truncation tol (outer iterations)

        [u_leg_left, u_leg_right, svec, res, flag, left, right] = it_LR_block_arnoldi(B0, M0, K0, K1, u0vec, H1aphi, maxit, tol, it_arn, trunc, FF, Tal, Usch, Hal, t0, tf);

        time_star = toc;
        u_leg = u_leg_left*u_leg_right.';
        X = right*left.';
        %% To compute the error with check the results
        err = NaN(m,1);
        for k = 1:m
            cheb_c = leg2cheb(full(u_leg(1:k,:)),'norm');
            u_star = chebfun(cheb_c, 'coeffs');
            err(k) = norm(u_star(1)-u(:,end).');
        end
        [err,k] = min(err);
        fprintf("%d & %d & %1.2e & %1.2e & %1.2f & %d \\\\\n",...
            m,k,err,err/norm(u(:,end).'),time_star,rank(X));
        fprintf(fid,"%d & %d & %1.2e & %1.2e & %1.2f & %d \\\\\n",...
            m,k,err,err/norm(u(:,end).'),time_star,rank(X));
    end

    %% Repeated solution by FLMM2
    if doflmm
    for h = [1e-1,1e-2,1e-3,1e-4]
        tic;
        jfun = @(t,u) jfunt(t,FEM,K0,K1-K0,ft,str);
        [tlist,uh] = flmm2(alpha, @(t,u) odefun(t,u),...
            @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
            []);
        time_flmm2 = toc;
        abserr = norm(uh(:,end)-u(:,end));
        relerr = abserr/norm(u(:,end));
        fprintf("%1.2e & %1.2f & %1.2e & %1.2e \\\\\n",h,time_flmm2,abserr,relerr);
        fprintf(fid,"%1.2e & %1.2f & %1.2e & %1.2e \\\\\n",h,time_flmm2,abserr,relerr);
    end
    else
	fprintf("Skipping flmm\n");
	fprintf(fid,"\\\\\n");
    end
end

fclose(fid);

%% Auxiliary functions
function V = V(x,t,f)
%V - Potential function for the time-dependent Schroedinger equation,
% I HAVE GENERALIZED THIS TO THE 2D CASE...
%  Walker-Preston model of a diatomic molecule in a strong laser field
%  - see https://www.ehu.eus/ccwmuura/research/SplitMagdef.pdf
%
n1 = 1;
nr = numel(x.x);
V = zeros(n1,nr);
for i = 1:nr
    if x.x(i) >= -1 && x.x(i) <= 1 && x.y(i) >= -1 && x.y(i) <= 1
        V(1,i) = 0;
    else
        V(1,i) = 10.0;
    end
    V(1,i) = V(1,i) + f(t.time); % A*cos(w*t.time);
end
end

function u0 = u0(x)
%u0 - Initial condition for the time-dependent Schroedinger equation
u0 = exp(-((x.x).^2 + (x.y).^2)/2);
end

function u1 = odefunt(t,u,FEM,K0,Kft,ft,str)
state.time = t;
switch str
    case "nullspace"
        Kc = (Kft)*(ft(t)-ft(0))/(ft(1)-ft(0)) + K0;
        u1 = FEM.B*(FEM.M\(- Kc*(FEM.B'*u) - FEM.Fc)) + FEM.ud;
    otherwise
	error("Non existing case");
end
end

function u1 = jfunt(t,FEM,K0,Kft,ft,str)
switch str
    case "nullspace"
        Kc = (Kft)*(ft(t)-ft(0))/(ft(1)-ft(0)) + K0;
        u1 = FEM.B*(FEM.M\(- Kc*FEM.B'));
    otherwise
	error("Non existing case");
end
end

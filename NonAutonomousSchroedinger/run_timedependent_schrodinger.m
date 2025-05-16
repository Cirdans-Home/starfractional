clear; clc; close all;

addpath('src')
addpath('chebfun');
fid = fopen('output_table.txt','a+');
fprintf(fid,"Date: %s\n\n",string(datetime));

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
    for m = [100,500,1000]
        tic
        fD{1} = @(t) (ft(t)-ft(0))/(ft(1)-ft(0));  % scalar function of the time-dependent part
        fT{1} = @(t) 1 + 0*t; % It will produce the Heaviside function

        trunc = 0;  % No automatic truncation of the last part of the coeff matrices
        [FlegT, ~] = genCoeffMatrix_SP_interval(fT,m,[t0,tf],trunc); % It generates the f(t)Theta(t-s) coeff matrix
        H = sparse(FlegT{2}*(tf-t0)/2); % Heaviside function coeff matrix
        [FlegD, ~] = genCoeffMatrix_delta_SP(fD,m,[t0,tf]); % It generates the f(t)delta(t-s) coeff matrix
        F2D = sparse(FlegD{1})*2/(tf-t0); % coeff matrix of fF{1}(t)delta(t-s)

        if verLessThan('matlab','24.2')
            Hal = expm(alpha*logm(H));       % Heaviside to the alpha
            H1mal = expm((1-alpha)*logm(H)); % Heaviside to the 1-alpha
        else
            Hal = H^alpha;  % Heaviside to the alpha
            H1mal = H^(1-alpha);  % Heaviside to the 1-alpha
        end
        % AA{1} = full(-(B0*(M0\(K0*B0'))).'); AA{2} = full(-(B0*(M0\((K1-K0)*B0'))).');
        FF{1} = Hal; FF{2} = F2D*Hal;  % coefficient matrix in the resolvent

        % Legendre poly vector at t=t0
        phi = zeros(m,1);
        for i=1:m
            phi(i,1) = (-1)^(i-1)*sqrt((2*(i-1)+1)/2);
        end
        H1aphi = H1mal*phi;

        d= length(u0vec);

        % Solution of the system by gmres
        % aprod = @(V) reshape(V - FF{1}*V*AA{1} - FF{2}*V*AA{2},m*d,1);
        tol = 1e-3;
        maxit = 300;
        aprod = @(V) reshape((V.' + B0*(M0\( K0*B0'*(V.'*FF{1}.') + (K1-K0)*B0'*(V.'*FF{2}.')))).',m*d,1);
        [x,flag,relres,nit,res] = gmres(@(v)aprod(reshape(v,m,d)), kron(u0vec,H1aphi),[],tol,maxit);
        X = reshape(x,m,d);  % solution of the corresponding matrix equation
        u_leg = Hal*X*2/(tf-t0); % leg coeff of the solution
        time_star = toc;
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
    for h = [1e-1,1e-2,1e-3,1e-4]
        tic;
        [tlist,uh] = flmm2(alpha, @(t,u) odefun(t,u),...
            @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
            []);
        time_flmm2 = toc;
        abserr = norm(uh(:,end)-u(:,end));
        relerr = abserr/norm(u(:,end));
        fprintf("%1.2e & %1.2f & %1.2e & %1.2e \\\\\n",h,time_flmm2,abserr,relerr);
        fprintf(fid,"%1.2e & %1.2f & %1.2e & %1.2e \\\\\n",h,time_flmm2,abserr,relerr);
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
    case "stiff-spring"
        % u1= FEM.M\(-FEM.Ks*u -FEM.Fs);
    case "nullspace"
        Kc = (Kft)*(ft(t)-ft(0))/(ft(1)-ft(0)) + K0;
        u1 = FEM.B*(FEM.M\(- Kc*(FEM.B'*u) - FEM.Fc)) + FEM.ud;
end
end

function u1 = jfunt(t,FEM,K0,Kft,ft,str)
% state.time = t;
% FEM = assembleFEMatrices(model,str,state);
switch str
    case "stiff-spring"
        % u1= FEM.M\(-FEM.Ks);
    case "nullspace"
        Kc = (Kft)*(ft(t)-ft(0))/(ft(1)-ft(0)) + K0;
        u1 = FEM.B*(FEM.M\(- Kc*FEM.B'));
end
end
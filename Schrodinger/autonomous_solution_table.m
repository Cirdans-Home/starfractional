%% Run the experiment for the table for the autonomous Schroedinger eq

clear; clc; close all;
addpath('src')
addpath('chebfun');

alpha = 0.5;
method = 3; % F-BDF2
tspan = [0,1];

fid = fopen('output_table.txt','a+');
fprintf(fid,"Date: %s\n\n",string(datetime));

for Hmax = [0.3,0.1,0.05]
    fprintf("Hmax = %1.2f\n",Hmax)
    fprintf(fid,"Hmax = %1.2f\n",Hmax);
    %% Spatial discretization
    model = createpde();
    gd = [3,4,-2,2,2,-2,-2,-2,2,2]';
    ns = char('R1');
    ns = ns';
    sf = 'R1';
    g = decsg(gd,sf,ns);
    geometryFromEdges(model,g);
    generateMesh(model,'Hmax',Hmax,'GeometricOrder','quadratic');
    specifyCoefficients(model,'m',0,...
        'd',1i^alpha,...
        'c',1/2,...
        'a',@V,...
        'f',0);
    setInitialConditions(model,@u0);
    applyBoundaryCondition(model,'dirichlet',...
        'Edge',1:model.Geometry.NumEdges,'u',0);
    FEM = assembleFEMatrices(model,"stiff-spring");
    odefun = @(t,u) FEM.M\(-FEM.Ks*u -FEM.Fs);
    A = -FEM.M\FEM.Ks;
    B = -FEM.M\FEM.Fs;
    fprintf("DoFs = %d\n",length(A))
    
    xnode = model.Mesh.Nodes(1,:)';
    ynode = model.Mesh.Nodes(2,:)';
    u0vec = u0(struct('x',xnode,...
        'y',ynode));
    
    %% Benchmark Solution
    % This is going to be slow!
    
    h = 1e-5;
    filename = sprintf("fracsolve_HMAX%1.3f_h%1.2e.mat",Hmax,h);
    if ~exist(filename,"file")
        jfun = @(t,u) A;
        [tlist,u] = flmm2(alpha, @(t,u) odefun(t,u),...
            @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
            [],method);
        savefile = sprintf("src/%s",filename);
        save(savefile,"tlist","u",'-v7.3')
    else
        load(filename)
    end
    
    
    for m = [100,500,1000,1500,2000]
        %% We then run the solution with star-lanczos
        % we pass a cutoff of -1 to get all the coefficient and check for the best
        % cutoff outside.
        k = -1;
        max_it = m;
        int = [tspan(1),tspan(2)];
        [u_leg,~,X,time_star] = starsolve_autonomous(m,max_it,int,alpha,A,u0vec,k);
        
        kvals = 1:m;
        abserror = NaN(length(kvals),1);
        relerror = abserror;
        for i = 1:length(kvals)
            cheb_c = leg2cheb(full(u_leg(1:kvals(i),:)),'norm');
            u_star = chebfun(cheb_c, 'coeffs');
            abserror(i) = norm(u_star(tspan(2))-u(:,end).');
            relerror(i) = abserror(i)/norm(u(:,end));
        end
        % Find minimum
        [~,minind] = min(relerror);
        k = kvals(minind);
        abserror = abserror(minind);
        relerror = relerror(minind);
        
        fprintf('star & %d & %d & %1.2e & %1.2e & %d & %1.2f \n',...
            m,k,abserror,relerror,rank(X),time_star)
        fprintf(fid,'star & %d & %d & %1.2e & %1.2e & %d & %1.2f \n',...
            m,k,abserror,relerror,rank(X),time_star);
    end
    % Check with F-BDF2 at different values of h
    for h = [1e-1,1e-2,1e-3,1e-4]
        jfun = @(t,u) A;
        tic;
        [tlist,uh] = flmm2(alpha, @(t,u) odefun(t,u),...
            @(t,u) jfun(t,u),tspan(1),tspan(2),u0vec,h,...
            [],method);
        timeflmm = toc;
        abserror = norm(uh(:,end)-u(:,end));
        relerror = abserror/norm(u(:,end));
        fprintf('flmm2 & %1.2e & %f & %1.2e & %1.2e \n',h,timeflmm,abserror,relerror)
        fprintf(fid,'flmm2 & %1.2e & %f & %1.2e & %1.2e \n',h,timeflmm,abserror,relerror);
    end
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

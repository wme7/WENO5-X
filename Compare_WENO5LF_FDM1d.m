%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Solving 1-D wave equation using various implementations of 5th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%                 du/dt + df/dx = S, for x \in [a,b]
%                  where f = f(u): linear/nonlinear
%                     and S = s(u): source term
%
%             coded by Manuel Diaz, manuel.ade'at'gmail.com 
%           Biomedical Simulation Laboratory, NHRI, 2017.05.20
%                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [WENO5-JS] Jiang, Guang-Shan, and Chi-Wang Shu. "Efficient implementation
%            of weighted ENO schemes." JCP 126.1 (1996): 202-228.  
% [WENO5-M ] Henrick, Andrew K., Tariq D. Aslam, and Joseph M. Powers.
%            "Mapped weighted essentially non-oscillatory schemes:
%            achieving optimal order near critical points." JCP 207.2
%            (2005): 542-567. 
% [WENO5-Z ] Borges, Rafael, et al. "An improved weighted essentially
%            non-oscillatory scheme for hyperbolic conservation laws." 
%            JCP 227.6 (2008): 3191-3211. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: Conservative finite difference implementations of the method of
% lines (MOL) using WENO5 associated with SSP-RK33 time integration method. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Special Thx to Dr. Lusher D. for suggesting comparing these formulations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;

%% Parameters
     nx = 0080;	% number of cells
    CFL = 0.20;	% Courant Number
   tEnd = 0.40; % End time

fluxfun='buckley'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1.0; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
    case 'buckley' % Buckley-Leverett, CFL_max: 0.20 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end

% Build discrete domain
a=-1; b=1; dx=(b-a)/nx; x=a+dx/2:dx:b; 

% Build IC
ICcase=2;  % {1}Testing, {2}Costum ICs
switch ICcase
    case 1 % Testing IC
        u0=TestingIC(x);  % Jiang and Shu IC
    case 2 % Guassian IC
        u0=CommonIC(x,9)-1; % cases 1-10 <- check them out!
    otherwise
        error('IC file not listed');
end

% Plot range
dl=0.1; plotrange=[a,b,min(u0)-dl,max(u0)+dl];

%% Solver Loop

% load initial conditions
t=0; it=0; u=u0; 

tic
while t < tEnd
	% Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end
    
	% Update time and iteration counter
    t=t+dt; it=it+1;
    
    % RK Initial step
    uo = u;
    
    % 1st stage
    dF = WENO5JSresAdv1d(u,flux,dflux,S,dx); 
    u = uo-dt*dF;
    
    % 2nd Stage
    dF = WENO5JSresAdv1d(u,flux,dflux,S,dx);
    u = 0.75*uo+0.25*(u-dt*dF);

    % 3rd stage
    dF = WENO5JSresAdv1d(u,flux,dflux,S,dx);
    u = (uo+2*(u-dt*dF))/3;
end
cputime=toc; fprintf('WENO-JS cputime: %g\n',cputime);

% save result and clear space
u_JS = u; clear t it u;

% load initial conditions
t=0; it=0; u=u0; 

tic
while t < tEnd
	% Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end
    
	% Update time and iteration counter
    t=t+dt; it=it+1;
    
    % RK Initial step
    uo = u;
    
    % 1st stage
    dF = WENO5MresAdv1d(u,flux,dflux,S,dx);
    u = uo-dt*dF;
    
    % 2nd Stage
    dF = WENO5MresAdv1d(u,flux,dflux,S,dx);
    u = 0.75*uo+0.25*(u-dt*dF);

    % 3rd stage
    dF = WENO5MresAdv1d(u,flux,dflux,S,dx);
    u = (uo+2*(u-dt*dF))/3;
end
cputime=toc; fprintf('WENO-M  cputime: %g\n',cputime);

% save result and clear space
u_M = u; clear t it u;

% load initial conditions
t=0; it=0; u=u0; 

tic
while t < tEnd
	% Update/correct time step
    dt=CFL*dx/max(abs(u)); if t+dt>tEnd, dt=tEnd-t; end
    
	% Update time and iteration counter
    t=t+dt; it=it+1;
    
    % RK Initial step
    uo = u;
    
    % 1st stage
    dF = WENO5ZresAdv1d(u,flux,dflux,S,dx);
    u = uo-dt*dF;
    
    % 2nd Stage
    dF = WENO5ZresAdv1d(u,flux,dflux,S,dx);
    u = 0.75*uo+0.25*(u-dt*dF);

    % 3rd stage
    dF = WENO5ZresAdv1d(u,flux,dflux,S,dx);
    u = (uo+2*(u-dt*dF))/3;
end
cputime=toc; fprintf('WENO-Z  cputime: %g\n',cputime);

% save result and clear space
u_Z = u; clear t it u;

%% Final Plot
plot(x,u0,'-',x,u_JS,'o',x,u_M,'s',x,u_Z,'d'); axis(plotrange);
legend('Initial Condition','WENO5-JS','WENO5-M','WENO5-Z'); legend boxoff;
title('WENO5 - cell averages plot','interpreter','latex','FontSize',18);
xlabel('$\it{x}$','interpreter','latex','FontSize',14);
ylabel({'$\it{u(x)}$'},'interpreter','latex','FontSize',14);

%% Conclusion
%
% The comparisons show that the WENO-JS has been superseded by WENO5-Z.
% A more extensive comparison for system of equations has been presented in
% Brehm, Christoph, et al. "A comparison of higher-order finite-difference
% shock capturing schemes." Computers & Fluids 122 (2015): 184-208. 
%
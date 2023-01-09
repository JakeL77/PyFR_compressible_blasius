clear;clc;close all

% Jake Lewis 2022
% Imperial College London

% Notation
% nu is the wall normal similarity coordinate, n the wall normal physical
% coordinate
% y is the similarity solution vector (see Oz and Kara 2021 for more
% details on notation)

%% USER SET FLOW CONDITIONS

c_2 = 110.4; % K, for Sutherland's
Pr = 0.71; % Prandtl number
gamma = 1.4; % adiabatic ratio

T_e = 54.8; % K, freestream temp
M_e = 7.73; % freestream Mach
u_e = 7.73; % freestream velocity
mu_e = 6.24194e-5; % freestream viscosity
rho_e = 1.4; % freestream density
T_w = 354.8; % K, wall temperature; unused if adiabatic

%x_0 = 0.57678; % distance along flat plate
x_0 = 0.2883;

wallCondition = "isothermal";
%wallCondition = "adiabatic";

%% USER SET NEWTON METHOD PARAMETERS

y30InitialGuess = 0.5; % initial guess for y30 boundary condition; 0.1 often works 
% for supersonic, 0.5 for hyp

% note y40InitialGuess unused if isothermal, y50InitialGuess unused if
% adiabatic

y40InitialGuess = 10.0; % initial guess for y40 boundary condition; 3.0 often works 
% for supersonic, 10.0 for hyp; as this is actually the adiabatic
% temperature ratio at the wall, we would of course expect this to increase
% with Mach; maybe put in a semi-analytical method to estimate this to
% accelerate convergence?

y50InitialGuess = 3.0; % initial guess for y50 boundary condition

derivativeIncrement = 1e-10; % delta for forward differencing 
% derivative approximation in bl_generator

newtonTol = 1e-9; % tolerance for Newton iterations

nuEnd = 20; % upper integration limit in similarity coordinate

%% GENERATING BOUNDARY LAYER PROFILE

% nu, uBar, rhoBar gives the profile in similarity solution coordinates, n,
% u, rho gives the solution in physical coordinates
[nu,uBar,rhoBar,n,u,rho] = bl_generator(y30InitialGuess,y40InitialGuess,y50InitialGuess,...
    derivativeIncrement,newtonTol,nuEnd,c_2,T_e,Pr,gamma,M_e,mu_e,rho_e,u_e,T_w,x_0,wallCondition);

%% CALCULATING V

% calculates vertical velocity profile using continuity; generates a
% profile (nInc,uInc, rhoInc) a distance deltaX further along the plate 
% than the original profile, then numerically integrates this using 
% continuity integrator to give v, which is the vertical velocity profile 
% in physical coordinates

deltaX = 1e-4;
[~,~,~,nInc,uInc,rhoInc] = bl_generator(y30InitialGuess,y40InitialGuess,y50InitialGuess,...
    derivativeIncrement,newtonTol,nuEnd,c_2,T_e,Pr,gamma,M_e,mu_e,rho_e,u_e,T_w,x_0+deltaX,wallCondition);

v = continuity_integrator(n,u,rho,nInc,uInc,rhoInc,deltaX);

%% GETTING BL PROPERTIES

% calculates some useful boundary layer properties; incompressible
% displacement thickness (deltaStarIC), incompressible momentum thickness
% (thetaIC), incompressible shape factor (HIC), and compressible
% displacement thickness (deltaStarC)
% May add compressible momentum thickness and shape factor
[deltaStarIC,thetaIC,HIC,deltaStarC,thetaC] = bl_properties(n,uBar,rhoBar);

% calculates Reynolds numbers based on these parameters
ReThetaIC = thetaIC*rho_e*u_e/mu_e;
ReDeltaStarIC = deltaStarIC*rho_e*u_e/mu_e;
ReDeltaStarC = deltaStarC*rho_e*u_e/mu_e;
ReThetaC = thetaC*rho_e*u_e/mu_e;

% finding delta 99
delta99 = interp1(u/u(end),n,0.99);

%% FITTING

% generates fits using num2fit, which utilises the MATLAB curve fitting
% toolbox; velFit, rhoFit, vFit contain anonymous functions, while
% velFitString etc contain strings suitable for input to PyFR
if wallCondition == "adiabatic"
    [velFit, velFitString] = fitArray(n,u,u_e,'sin8',[],0.0);
    [rhoFit, rhoFitString] = fitArray(n,rho,rho_e,'generalLogistic',[1230,66,1.43,7.6],0.0);
    [vFit,vFitString] = fitArray(n,v,v(end),'generalLogisticUpper',[750,38,10,13,0],0.0);
elseif wallCondition == "isothermal"
    [velFit, velFitString] = fitArray(n,u,u_e,'splitPolySqrt',[0.1,-4.2,130,4.3],0.025);
    [rhoFit, rhoFitString] = fitArray(n,rho,rho_e,'splitPolyGeneralLogistic',[1230,66,23,1.43,0],0.025);
    [vFit,vFitString] = fitArray(n,v,v(end),'splitPolyGeneralLogistic',[750,38,10,13,0],0.025);
end

%% CALCULATING FIRST NODE HEIGHT

% calculates a first node height based on y+, using the gradient at the
% wall (NOT the gradient of the fit)
% in PyFR the height of the first cell will not be this first node point,
% but will depend on the order of the solution polynomial and the point set
% used
muBarWall = (rhoBar(1)^1.5)*((1.0+c_2/T_e)/(rhoBar(1)+c_2/T_e));
muWall = muBarWall*mu_e;
dudyWall = (u(2)-u(1))/(n(2)-n(1));
tauWall = muWall*dudyWall;
uTauWall = sqrt(tauWall/rho(1));
targetYplus = 0.5;
nodeHeightWall = muWall*targetYplus/(uTauWall*rho(1));

%% OUTPUT
fprintf('Generated boundary layer properties:\n')
fprintf('Re(theta, incompressible) = %.3f\n',ReThetaIC)
fprintf('Re(delta, incompressible) = %.3f\n',ReDeltaStarIC)
fprintf('Re(delta, compressible) = %.3f\n',ReDeltaStarC)
fprintf('Re(theta, compressible) = %.3f\n',ReThetaC)
fprintf('Shape factor (incompressible) = %.3f\n',HIC)
fprintf('99%% thickness = %.3f\n',delta99)
fprintf('\n')
fprintf('Generated fit strings for PyFR:\n')
fprintf('u:\n%s\n',velFitString)
fprintf('v:\n%s\n',vFitString)
fprintf('rho:\n%s\n',rhoFitString)

%% PLOTTING

figure()
hold on
title('Streamwise velocity profile','FontSize',35)
grid minor
nPlot = 0:1e-5:1.5*delta99; % physical coord range to plot fit over
ylim([0 1.5])
xlim([0 1.1])
plot(velFit(nPlot)/u_e,nPlot/delta99,'-k')
scatter(u/u_e,n/delta99,'x')
legend('Fit','Numerical Solution')
set(gca,'FontSize',20)
ylabel('$y/\delta_{99}$','Interpreter','latex','FontSize',35)
xlabel('$U/U_e$','Interpreter','latex','FontSize',35)
hold off

figure()
hold on
title('Static density profile','FontSize',35)
grid minor
nPlot = 0:1e-5:1.5*delta99;
ylim([0 1.5])
xlim([0.5*rho(1) 1.1])
plot(rhoFit(nPlot)/rho_e,nPlot/delta99,'-k')
scatter(rho/rho_e,n/delta99,'x')
legend('Fit','Numerical Solution')
set(gca,'FontSize',20)
xlabel('$\rho/\rho_e$','Interpreter','latex','FontSize',35)
ylabel('$y/\delta_{99}$','Interpreter','latex','FontSize',35)
hold off

figure()
hold on
title('Wall normal velocity profile','FontSize',35)
grid minor
nPlot = 0:1e-5:0.08;
ylim([0 1.5])
xlim([0 1.1])
plot(vFit(nPlot)/v(end),nPlot/delta99,'-k')
scatter(v/v(end),n/delta99,'x')
legend('Fit','Numerical Solution')
set(gca,'FontSize',20)
xlabel('$V/V_e$','Interpreter','latex','FontSize',35)
ylabel('$y/\delta_{99}$','Interpreter','latex','FontSize',35)
hold off

figure()
hold on
title('Combined plot (of fits)','FontSize',35)
grid minor
nPlot = 0:1e-5:delta99*1.5;
ylim([0 1.5])
xlim([0 1.1])
plot(velFit(nPlot)/u_e,nPlot/delta99,'-r')
plot(rhoFit(nPlot)/rho_e,nPlot/delta99,'-b')
plot(vFit(nPlot)/v(end),nPlot/delta99,'-k')
set(gca,'FontSize',20)
xlabel('$U/U_e,\rho/\rho_e,V/V_e$','Interpreter','latex','FontSize',35)
ylabel('$y/\delta_{99}$','Interpreter','latex','FontSize',35)
legend('$U/U_e$','$\rho/\rho_e$','$V/V_e$','Interpreter','latex','FontSize',35)
hold off
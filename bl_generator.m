function [nu,uBar,rhoBar,n,u,rho] = bl_generator(y30InitialGuess,y40InitialGuess,...
    derivativeIncrement,newtonTol,nuEnd,c_2,T_e,Pr,gamma,M_e,mu_e,rho_e,u_e,x_0)

% Using Newton's method to obtain correct y30 and y40 boundary conditions

% getting initial errors
[~,y] = bl_integrator(y30InitialGuess,y40InitialGuess,nuEnd,c_2,T_e,Pr,gamma,M_e);
y2old = y(end,2);
y4old = y(end,4);
newtonError2 = 1-y2old;
newtonError4 = 1-y4old;

y30old = y30InitialGuess;
y40old = y40InitialGuess;

% Newton loop
while abs(newtonError2) > newtonTol || abs(newtonError4) > newtonTol
    % increments boundary conditons slightly to get derivatives 
    [~,y] = bl_integrator(y30old+derivativeIncrement,y40old,nuEnd,c_2,T_e,Pr,gamma,M_e);
    y2new1 = y(end,2);
    y4new1 = y(end,4);
    
    [~,y] = bl_integrator(y30old,y40old+derivativeIncrement,nuEnd,c_2,T_e,Pr,gamma,M_e);
    y2new2 = y(end,2);
    y4new2 = y(end,4);
    
    % solves for new boundary conditions increments
    derivativeMatrix = [(y2new1-y2old)/derivativeIncrement (y2new2-y2old)/derivativeIncrement;
        (y4new1-y4old)/derivativeIncrement (y4new2-y4old)/derivativeIncrement];
    dy0 = derivativeMatrix \ [1.0-y2old;1.0-y4old];
    
    % updates boundary conditions
    y30old = y30old+dy0(1);
    y40old = y40old+dy0(2);
    
    % finds new errors
    [nu,y] = bl_integrator(y30old,y40old,nuEnd,c_2,T_e,Pr,gamma,M_e);
    y2old = y(end,2);
    y4old = y(end,4);
    newtonError2 = 1-y2old;
    newtonError4 = 1-y4old;
end

%% FINAL PARAMETERS

% sets found boundary conditions
y30Final = y30old;
y40Final = y40old;

% sets final similarity coordinate vector and corresponding
% solution vector
nuFinal = nu;
yFinal = y;

% needed for coordinate transform
rhoArrayIntermediate = rho_e./yFinal(:,4);

% housekeeping
clear nu y;

%% TRANSFORMING COORDS
% calculating streamwise similarity coordinate
s = mu_e*rho_e*u_e*x_0;
% using this to transform similarity wall normal coordinate into physical
% wall normal coordinate n; nu contains the wall normal similarity values
% used in the integration
[nu,n] = coord_transform_integrator(s,u_e,nuFinal,rhoArrayIntermediate);

% interpolating solution vector onto physical coordinate vector (as yFinal
% corresponds to nuFinal, not nu)
% this probably needs a bit of a cleanup to make clearer
yTransformed = interp1(nuFinal,yFinal,nu);

% assigning physically relevant variable names to solution vector
% components
uBar = yTransformed(:,2);
rhoBar = yTransformed(:,4);
u = yTransformed(:,2)*u_e;
rho = rho_e./yTransformed(:,4);

end
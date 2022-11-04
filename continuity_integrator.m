function v = continuity_integrator(n,u,rho,nInc,uInc,rhoInc,deltaX)
% uses two boundary layer profiles separated by slight streamwise
% difference to obtain vertical velocity profile applying continuity

% has to interpolate the incremented profile so it is mapped to the same
% physical coordinates of the original profile
uInc = interp1(nInc,uInc,n);
rhoInc = interp1(nInc,rhoInc,n);

% forward differencing to obtain gradient
drhoudx = (uInc.*rhoInc-u.*rho)/deltaX;

% numerically integrating
rhov = zeros(length(n),1);
for i = 2:length(n)
    rhov(i) = rhov(i-1)-(n(i)-n(i-1))*drhoudx(i-1); 
end
v = rhov./rho;
end
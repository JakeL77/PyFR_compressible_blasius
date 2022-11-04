function [deltaStarIC,thetaIC,HIC,deltaStarC] = bl_properties(n,uBar,rhoBar)
% numerically integrating for incompressible bl quantities
deltaStarIntegrandIC = 1-uBar;
deltaStarIC = trapz(n,deltaStarIntegrandIC);

thetaIntegrandIC = (uBar).*(1-(uBar));
thetaIC = trapz(n,thetaIntegrandIC);
HIC = deltaStarIC/thetaIC;

% numerically integrating for compressible bl quantities
deltaStarIntegrandC = (1-uBar./rhoBar);
deltaStarC = trapz(n,deltaStarIntegrandC);
end
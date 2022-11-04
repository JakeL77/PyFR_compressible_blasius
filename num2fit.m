function [velFit, velFitString, rhoFit, rhoFitString, vFit, vFitString] = num2fit(n,velocityXArray,rhoArray,vArray,u_e,rho_e,v_e)

% uses MATLAB curve fitting toolbox; 8 term fourier series works well at 
% low Mach, but at higher Mach less so; a different fit may be desirable
fourierVelFit = fit(n,velocityXArray,'fourier8');

% obtains coefficient values
velFitCoeffs = coeffvalues(fourierVelFit);
a0Vel = velFitCoeffs(1);
a1Vel = velFitCoeffs(2);
b1Vel = velFitCoeffs(3);
a2Vel = velFitCoeffs(4);
b2Vel = velFitCoeffs(5);
a3Vel = velFitCoeffs(6);
b3Vel = velFitCoeffs(7);
a4Vel = velFitCoeffs(8);
b4Vel = velFitCoeffs(9);
a5Vel = velFitCoeffs(10);
b5Vel = velFitCoeffs(11);
a6Vel = velFitCoeffs(12);
b6Vel = velFitCoeffs(13);
a7Vel = velFitCoeffs(14);
b7Vel = velFitCoeffs(15);
a8Vel = velFitCoeffs(16);
b8Vel = velFitCoeffs(17);
wVel = velFitCoeffs(18);

% builds up fit function
fourierVelFitFun = @(x) a0Vel + a1Vel*cos(x*wVel) + b1Vel*sin(x*wVel) + ...
    a2Vel*cos(2*x*wVel) + b2Vel*sin(2*x*wVel) + a3Vel*cos(3*x*wVel) + b3Vel*sin(3*x*wVel) + ...
    a4Vel*cos(4*x*wVel) + b4Vel*sin(4*x*wVel) + a5Vel*cos(5*x*wVel) + b5Vel*sin(5*x*wVel) + ...
    a6Vel*cos(6*x*wVel) + b6Vel*sin(6*x*wVel) + a7Vel*cos(7*x*wVel) + b7Vel*sin(7*x*wVel) + ...
    a8Vel*cos(8*x*wVel) + b8Vel*sin(8*x*wVel);
velFit = @(x) fourierVelFitFun(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
    u_e*0.5.*(tanh(1e100*(x-n(end)))+1);

aVelCoeffs = velFitCoeffs(2:2:16);
bVelCoeffs = velFitCoeffs(3:2:17);

% constructs string suitable for PyFR
velFitString = ['(',num2str(a0Vel),'+'];
for i = 1:7
    velFitString = [velFitString,num2str(aVelCoeffs(i)),'*cos(',num2str(wVel),'*',num2str(i),'*y)+',...
        num2str(bVelCoeffs(i)),'*sin(',num2str(wVel),'*',num2str(i),'*y)+'];
end
velFitString = [velFitString,num2str(aVelCoeffs(8)),'*cos(',num2str(wVel),'*',num2str(8),'*y)+',...
        num2str(bVelCoeffs(8)),'*sin(',num2str(wVel),'*',num2str(8),'*y)'];
velFitString = [velFitString,')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
    num2str(u_e),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];

% performs same process for density
fourierRhoFit = fit(n,rhoArray,'fourier8');
rhoFitCoeffs = coeffvalues(fourierRhoFit);
a0rho = rhoFitCoeffs(1);
a1rho = rhoFitCoeffs(2);
b1rho = rhoFitCoeffs(3);
a2rho = rhoFitCoeffs(4);
b2rho = rhoFitCoeffs(5);
a3rho = rhoFitCoeffs(6);
b3rho = rhoFitCoeffs(7);
a4rho = rhoFitCoeffs(8);
b4rho = rhoFitCoeffs(9);
a5rho = rhoFitCoeffs(10);
b5rho = rhoFitCoeffs(11);
a6rho = rhoFitCoeffs(12);
b6rho = rhoFitCoeffs(13);
a7rho = rhoFitCoeffs(14);
b7rho = rhoFitCoeffs(15);
a8rho = rhoFitCoeffs(16);
b8rho = rhoFitCoeffs(17);
wrho = rhoFitCoeffs(18);
fourierRhoFitFun = @(x) a0rho + a1rho*cos(x*wrho) + b1rho*sin(x*wrho) + ...
    a2rho*cos(2*x*wrho) + b2rho*sin(2*x*wrho) + a3rho*cos(3*x*wrho) + b3rho*sin(3*x*wrho) + ...
    a4rho*cos(4*x*wrho) + b4rho*sin(4*x*wrho) + a5rho*cos(5*x*wrho) + b5rho*sin(5*x*wrho) + ...
    a6rho*cos(6*x*wrho) + b6rho*sin(6*x*wrho) + a7rho*cos(7*x*wrho) + b7rho*sin(7*x*wrho) + ...
    a8rho*cos(8*x*wrho) + b8rho*sin(8*x*wrho);
rhoFit = @(x) fourierRhoFitFun(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
    rho_e*0.5.*(tanh(1e100*(x-n(end)))+1);

rhoFitString = '';

aRhoCoeffs = rhoFitCoeffs(2:2:16);
bRhoCoeffs = rhoFitCoeffs(3:2:17);
rhoFitString = ['(',num2str(a0rho),'+'];
for i = 1:7
    rhoFitString = [rhoFitString,num2str(aRhoCoeffs(i)),'*cos(',num2str(wrho),'*',num2str(i),'*y)+',...
        num2str(bRhoCoeffs(i)),'*sin(',num2str(wrho),'*',num2str(i),'*y)+'];
end
rhoFitString = [rhoFitString,num2str(aRhoCoeffs(8)),'*cos(',num2str(wrho),'*',num2str(8),'*y)+',...
        num2str(bRhoCoeffs(8)),'*sin(',num2str(wrho),'*',num2str(8),'*y)'];
rhoFitString = [rhoFitString,')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
    num2str(rho_e),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];

% and for wall normal velocity
fourierVFit = fit(n,vArray,'fourier8');
vFitCoeffs = coeffvalues(fourierVFit);
a0v = vFitCoeffs(1);
a1v = vFitCoeffs(2);
b1v = vFitCoeffs(3);
a2v = vFitCoeffs(4);
b2v = vFitCoeffs(5);
a3v = vFitCoeffs(6);
b3v = vFitCoeffs(7);
a4v = vFitCoeffs(8);
b4v = vFitCoeffs(9);
a5v = vFitCoeffs(10);
b5v = vFitCoeffs(11);
a6v = vFitCoeffs(12);
b6v = vFitCoeffs(13);
a7v = vFitCoeffs(14);
b7v = vFitCoeffs(15);
a8v = vFitCoeffs(16);
b8v = vFitCoeffs(17);
wv = vFitCoeffs(18);
fourierVFitFun = @(x) a0v + a1v*cos(x*wv) + b1v*sin(x*wv) + ...
    a2v*cos(2*x*wv) + b2v*sin(2*x*wv) + a3v*cos(3*x*wv) + b3v*sin(3*x*wv) + ...
    a4v*cos(4*x*wv) + b4v*sin(4*x*wv) + a5v*cos(5*x*wv) + b5v*sin(5*x*wv) + ...
    a6v*cos(6*x*wv) + b6v*sin(6*x*wv) + a7v*cos(7*x*wv) + b7v*sin(7*x*wv) + ...
    a8v*cos(8*x*wv) + b8v*sin(8*x*wv);
vFit = @(x) fourierVFitFun(x)*0.5.*(tanh(-1e100*(x-n(end)))+1)+...
    v_e*0.5.*(tanh(1e100*(x-n(end)))+1);

vFitString = '';

avCoeffs = vFitCoeffs(2:2:16);
bvCoeffs = vFitCoeffs(3:2:17);
vFitString = ['(',num2str(a0v),'+'];
for i = 1:7
    vFitString = [vFitString,num2str(avCoeffs(i)),'*cos(',num2str(wv),'*',num2str(i),'*y)+',...
        num2str(bvCoeffs(i)),'*sin(',num2str(wv),'*',num2str(i),'*y)+'];
end
vFitString = [vFitString,num2str(avCoeffs(8)),'*cos(',num2str(wv),'*',num2str(8),'*y)+',...
        num2str(bvCoeffs(8)),'*sin(',num2str(wv),'*',num2str(8),'*y)'];
vFitString = [vFitString,')*0.5*(tanh(-10000000*(y-',num2str(n(end)),'))+1)+',...
    num2str(v_e),'*0.5*(tanh(10000000*(y-',num2str(n(end)),'))+1)'];

end
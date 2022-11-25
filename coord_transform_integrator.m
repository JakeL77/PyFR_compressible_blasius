function [nu,n] = coord_transform_integrator(s,u_e,nuFinal,rhoArray,nuEnd)
% uses ode45 to numerically integrate to obtain physical coordinates from
% similarity coordinates
    options = odeset('RelTol',1e-10,'AbsTol',1e-10);
    [nu,n] = ode45(@(nu,n) coord_transform_derivative(nu,n,s,u_e,nuFinal,rhoArray),[0 nuEnd],0,options);
end
function [nu,n] = coord_transform_integrator(s,u_e,nuFinal,rhoArray)
% uses ode45 to numerically integrate to obtain physical coordinates from
% similarity coordinates
    options = odeset('RelTol',eps,'AbsTol',eps);
    [nu,n] = ode45(@(nu,n) coord_transform_derivative(nu,n,s,u_e,nuFinal,rhoArray),[0 10],0,options);
end
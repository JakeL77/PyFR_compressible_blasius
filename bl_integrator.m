function [nu,y] = bl_integrator(y30,y40,nuEnd,c_2,T_e,Pr,gamma,M_e)

    % uses MATLAB ode45 to numerically integrate a similarity solution
    % boundary layer profile for some give y30, y40 boundary conditions
    % note nuEnd is the end of the similarity solution domain (edge of BL), 
    % and should be large enough that freestream conditions are valid there

    nuInterval = [0 nuEnd];
    y0 = [0,0,y30,y40,0];
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [nu,y] = ode45(@(nu,y) similarity_derivative(nu,y,c_2,T_e,Pr,gamma,M_e),nuInterval,y0,...
       options);
end
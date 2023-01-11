function [nu,y] = bl_integrator(y30,y40,y50,nuEnd,c_2,T_e,T_w,Pr,gamma,M_e,wall_condition,viscosityLaw)

    % uses MATLAB ode45 to numerically integrate a similarity solution
    % boundary layer profile for some give y30, y40 boundary conditions
    % note nuEnd is the end of the similarity solution domain (edge of BL), 
    % and should be large enough that freestream conditions are valid there

    nuInterval = [0 nuEnd];
    if wall_condition == "adiabatic"
        y0 = [0,0,y30,y40,0];
    elseif wall_condition == "isothermal"
        y0 = [0,0,y30,T_w/T_e,y50];
    end
    options = odeset('RelTol',1e-8,'AbsTol',1e-8);
    [nu,y] = ode45(@(nu,y) similarity_derivative(nu,y,c_2,T_e,Pr,gamma,M_e,viscosityLaw),nuInterval,y0,...
       options);
end
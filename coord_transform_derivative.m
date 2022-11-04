function dndnu = coord_transform_derivative(nu,n,s,u_e,nuPoints,rhoPoints)
% integrated by coord_transform_integrator
    rho = interp1(nuPoints,rhoPoints,nu);
    dndnu = sqrt(2*s)/(u_e*rho);
end
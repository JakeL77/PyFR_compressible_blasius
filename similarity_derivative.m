function dydnu = similarity_derivative(nu,y,c_2,T_e,Pr,gamma,M_e)
    % integrated when bl_integrator called
    dydnu = zeros(5,1);
    dydnu(1) = y(2);
    dydnu(2) = y(3);
    dydnu(3) = -y(3)*( y(5)/(2*y(4)) - y(5)/(y(4)+c_2/T_e) )-...
        y(1)*y(3)*( (y(4)+c_2/T_e)/(sqrt(y(4))*(1+c_2/T_e)) );
    dydnu(4) = y(5);
    dydnu(5) = -y(5)^2*(1/(2*y(4))-1/(y(4)+c_2/T_e))-...
        Pr*y(1)*y(5)/sqrt(y(4))*(y(4)+c_2/T_e)/(1+c_2/T_e)-...
        (gamma-1)*Pr*M_e^2*y(3)^2;
end
function [C_p] = power_factor(P_rotor, rho, V_1, A_rotor)
% step 11 of BEM. calculates power factor only. Inputs (P_rotor, rho, V_1,
% A_rotor) outputs C_p

C_p = P_rotor /(0.5*rho*V_1^3*A_rotor);


end


% Part 5b of BEM iteration inputs (C_L, C_D, phi)

function C_r = rotational_coefficiant(C_L, C_D, phi)
% angle should be in rads
% computes the rotational component of the coefficiants (ref BEM slide)
C_r = C_L * sin(phi) - C_D * cos(phi);
end
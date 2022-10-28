%% Find azimuthal average
function [U_rel] = azimuthal_relative_velocity(U_0, yaw, omega, r, a_prime, azimuthal_angle)
    U_rel = sqrt((U_0*cos(yaw) - U_0*(1+K * r_R * sin(azimuthal_angle)))^2 + (omega * r * (1+a_prime) - U_0 *sin(yaw) * cos(azimuthal_angle))^2);
end

%% Find azimuthal average
% function [U_rel] = azimuthal_relative_velocity(U_0, yaw, omega, r_R, R, a_prime, azimuthal_angle, K)
%     U_rel = sqrt((U_0*cos(yaw) - U_0*(1+K * r_R * sin(azimuthal_angle)))^2 + (omega * (R * r_R) * (1+a_prime) - U_0 *sin(yaw) * cos(azimuthal_angle))^2);
% end

function [U_rel] = azimuthal_relative_velocity(U_0, yaw, u)
    U_rel = (U_0*cos(yaw)- u);
end
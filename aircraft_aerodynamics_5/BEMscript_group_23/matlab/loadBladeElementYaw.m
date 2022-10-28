%% Find Loads for each blade element
function [F, P] = loadBladeElementYaw(R, u, yaw, U_0, skew_angle)
    rho = 1.225; %kg/m^3
    F = rho * pi * R^2 * 2 * u * sqrt((U_0 * cos(yaw)-u)^2 + (U_0 * sin(yaw))^2);
    P = rho * pi * R^2 * 2 * u * sqrt((U_0 * cos(yaw)-u)^2 + (U_0 * sin(yaw))^2) * (U_0 * cos(yaw) - u);
    K = coleman_model(skew_angle);
end

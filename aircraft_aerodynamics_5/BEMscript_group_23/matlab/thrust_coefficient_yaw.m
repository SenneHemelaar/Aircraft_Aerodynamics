%% find skew angle from induction factor and yaw
function [CT] = thrust_coefficient_yaw(a, a_initial, yaw)
    CT = 4 * a * sqrt(1-a_initial * (2 * cos(yaw)-a_initial));
%% find skew angle from induction factor and yaw
function [CT] = thrust_coefficient(a)
%     CT = 4 * a * (cos(yaw) + sin(yaw)* tan(skew_angle/2)-a * sec(skew_angle/2)^2);
    CT = 4 * a *(1-a);
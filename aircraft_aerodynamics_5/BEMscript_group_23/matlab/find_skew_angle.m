%% find skew angle from induction factor and yaw
function [skew_angle] = find_skew_angle(a, yaw)
    skew_angle = (0.6*a + 1) * yaw;
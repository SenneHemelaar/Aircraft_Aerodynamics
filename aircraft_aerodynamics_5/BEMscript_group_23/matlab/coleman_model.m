%% Find Loads for each blade element
function [K] = coleman_model(skew_angle)
   K =  2 * tan(skew_angle/2);
end

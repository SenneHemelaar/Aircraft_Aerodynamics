function [J, M] = optimization_function(Forces, q, k)
%OPTIMIZATION_FUNCTION
S = [Forces.strip1(:,8);Forces.strip2(:,8)] .* [Forces.strip1(:,4); Forces.strip2(:,4)] * q;
M = sum(S .* [Forces.strip1(:,2);Forces.strip2(:,2)]);
J  = k * Forces.CDind/ 0.0151 + (1-k) *  M/ 2.693084059735463e+06;
end


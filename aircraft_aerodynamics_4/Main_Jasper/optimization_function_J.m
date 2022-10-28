function [J, M] = optimization_function_J(Forces, q, k, x, winglet)
L = 14;
%OPTIMIZATION_FUNCTION
S1 = Forces.strip1(:,8).* Forces.strip1(:,4)*q;
M1 = sum(S1 .* Forces.strip1(:,2));

M2 = 0;
M2x = 0;
if winglet
    %S2 = Forces.strip2(:,8).* Forces.strip2(:,4)*q;
    %M2 = sum(S2 .* Forces.strip2(:,2));
    S2 = Forces.strip2(:,7).* Forces.strip2(:,4)*q;
    M2z = sum(S2 .* Forces.strip2(:,2) * cos(pi/2-x(2)));
    M2x = sum((Forces.strip2(:,2)-L)/tan(x(2)) .*S2 *sin(pi/2-x(2)));
    M2 = sum([M2z M2x]);
end 

if isnan(M2x)
    M = M1;
else
    M = sum([M1 M2]);
end

J  = k * (Forces.CDind/ 0.0151) + (1-k) *  (M/ 2.853024411148614e+06);
end


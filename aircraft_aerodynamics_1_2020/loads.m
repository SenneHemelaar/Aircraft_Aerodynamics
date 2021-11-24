function [CL, dp] = loads(gamma, rho, U_inf, chord, N)
%LOADS

DL = zeros(1,N);
dp = zeros(1,N);
for i=1:N
   DL(i) = rho*U_inf*gamma(i);
   dp(i) = DL(i)/chord/N;
end
L=sum(DL);

CL=L/(0.5*rho*U_inf^2*chord);
    
end


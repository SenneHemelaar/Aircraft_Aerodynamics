function [CL, Dp] = calculating_loads(gamma,rho,U_0,chord,N)

DL = zeros(1,N);
Dp = zeros(1,N);

for i=1:N
   DL(i) = rho*U_0*gamma(i);
   Dp(i) = DL(i)/chord/N;
end

L=sum(DL);

CL=L/(0.5*rho*U_0^2*chord);

end


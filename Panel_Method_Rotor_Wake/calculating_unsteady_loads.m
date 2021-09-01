function [CL, Dp] = calculating_unsteady_loads(gamma, v_g, rho,U_0,chord,N, alpha, dt)

DL = zeros(1,N);
Dp = zeros(1,N);

% length_i = chord /N;
length_i = sqrt((v_g.x(2:end) - v_g.x(1:end-1)).^2 + (v_g.z(2:end) - v_g.z(1:end-1)).^2);

for i=1:N
    Dp(i) = rho * (U_0 * gamma(i)/length_i(i) + (gamma(i,end) - gamma(i,end-1))/dt);
end
L = sum(Dp*length_i');

% for i=1:N
%    DL(i) = rho*U_0*gamma(i);
%    Dp(i) = DL(i)/chord/N;
% end


CL=L/(0.5*rho*U_0^2*chord);

end


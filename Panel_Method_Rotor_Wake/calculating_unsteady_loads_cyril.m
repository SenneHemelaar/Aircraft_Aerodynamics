function [CL, Dp] = calculating_unsteady_loads_cyril(gamma,t_n, dt, v_g, rho, U_0, chord, N)


DL = zeros(1,N);
Dp = zeros(1,N);
D_panel = zeros(1,N);

for i=1:N % loop over panels
    
    D_panel(i) = sqrt((v_g.x(i+1)-v_g.x(i))^2 + (v_g.z(i+1)-v_g.z(i))^2) ; 
    
    if t_n == 1
        Dp(i) = rho*(U_0*gamma(i,t_n)/D_panel(i) + (sum(gamma(1:i,t_n)))/dt);
    else
        Dp(i) = rho*(U_0*gamma(i,t_n)/D_panel(i) + (sum(gamma(1:i,t_n))-sum(gamma(1:i,t_n-1)))/dt);
    end 
   
    DL(i) = Dp(i)*(v_g.x(i+1)-v_g.x(i));
   
end %end loop over panels

L=sum(DL);

CL=L/(0.5*rho*U_0^2*chord);

end

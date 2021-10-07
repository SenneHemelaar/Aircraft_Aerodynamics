function [L,d_P] = loads(gamma,U_inf,rho,af_geo)
a = [af_geo.a_U,af_geo.a_L];

d_L = rho*U_inf.*gamma;
L = sum(d_L);

d_P = rho*U_inf.*gamma./a;
    
end

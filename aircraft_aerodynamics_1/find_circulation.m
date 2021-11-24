function [Gamma,a] = find_circulation(af_geo,n_panels,U_inf)

Gamma = zeros(1,2*n_panels);

% Control Points
CP_xL = af_geo.CP_xL; CP_zL = af_geo.CP_zL;
CP_xU = af_geo.CP_xU; CP_zU = af_geo.CP_zU;
CP_x = [CP_xL CP_xU]; CP_z = [CP_zL CP_zU];
% CP_x = flip(CP_x); CP_z = flip(CP_z);
% Vortex Points
VP_xL = af_geo.VP_xL; VP_zL = af_geo.VP_zL;
VP_xU = af_geo.VP_xU; VP_zU = af_geo.VP_zU;
VP_x = [VP_xL VP_xU]; VP_z = [VP_zL VP_zU];
% VP_x = flip(VP_x); VP_z = flip(VP_z);
% Normal Vectors
N_xL = af_geo.N_xL; N_zL = af_geo.N_zL;
N_xU = af_geo.N_xU; N_zU = af_geo.N_zU;
N_x = [N_xL N_xU]; N_z = [N_zL N_zU];
% N_x = flip(N_x); N_z = flip(N_z);
% Build arrays
u=zeros(2*n_panels,2*n_panels);
w=zeros(2*n_panels,2*n_panels);
a=zeros(2*n_panels,2*n_panels);
RHS=zeros(1,2*n_panels);

for i=1:2*n_panels        % Loop over control points  
    for j=1:2*n_panels    % Loop over vortex points 
        [u(i,j),w(i,j)]= VOR2D(1,CP_x(i),CP_z(i),VP_x(j),VP_z(j));
        a(i,j)=u(i,j)*N_x(i)+w(i,j)*N_z(i);
    end
    if i<=n_panels
        RHS(i)=-(U_inf*N_x(i)+U_inf*N_z(i));
    else
        RHS(i)=-(U_inf*N_x(i)-U_inf*N_z(i));
    end
end

Gamma=-linsolve(a,RHS'); 

end


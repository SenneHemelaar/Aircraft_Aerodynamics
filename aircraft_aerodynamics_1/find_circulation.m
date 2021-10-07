function [Gamma] = find_circulation(af_geo)

% Control Points
CP_xL = af_geo.CP_xL; CP_zL = af_geo.CP_zL;
CP_xU = af_geo.CP_xU; CP_zU = af_geo.CP_zU;
CP_x = [CP_xL CP_xU]; CP_z = [CP_zL CP_zU];

% Vortex Points
VP_xL = af_geo.VP_xL; VP_zL = af_geo.VP_zL;
VP_xU = af_geo.VP_xU; VP_zU = af_geo.VP_zU;
VP_x = [VP_xL VP_xU]; VP_z = [VP_zL VP_zU];

% Normal vortices
u=zeros(1,N);
w=zeros(1,N);
a=zeros(1,N);
RHS=zeros(1,N);

for i=1:n_panels        % Loop over control points  
    for j=1:n_panels    % Loop over vortex points 
        [u(i,j),w(i,j)]= VOR2D(1,CP_x(i),CP_z(i),VP_x(j),VP_z(j));
        a(i,j)=u(i,j)*n_x(i)+w(i,j)*n_z(i);
    end
    RHS(i)=-(U_inf*n_x(i)+W_inf*n_z(i));
end

end


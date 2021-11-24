function [gamma] = find_circulation(af_geo, N, u_inf, w_inf)
%FIND_CIRCULATION
u=zeros(1,N);
w=zeros(1,N);
a=zeros(1,N);
RHS=zeros(1,N);

for i=1:N % loop over vortex points  
    for j=1:N % loop over vortex points
        [u(i,j),w(i,j)]= VOR2D(af_geo.CP_x(i), af_geo.CP_z(i),...
                               af_geo.VP_x(j), af_geo.VP_z(j), 1);
        a(i,j)=u(i,j)*af_geo.N_x(i)+w(i,j)*af_geo.N_z(i);
    end
    RHS(i)=-(u_inf*af_geo.N_x(i)+w_inf*af_geo.N_z(i));
end

% solve the linear system
RHS = RHS';    % transpose to make suitable for linear solve
gamma=linsolve(a,RHS); 

end


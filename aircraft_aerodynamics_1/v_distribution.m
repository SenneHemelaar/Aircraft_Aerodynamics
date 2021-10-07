function [vx,vz,vabs] = v_distribution(U_inf,grid_density,plot_limit,gamma,af_geo)

gamma_x = [af_geo.VP_xL,flip(af_geo.VP_xU)];
gamma_z = [af_geo.VP_zL,flip(af_geo.VP_zU)];

%constructing grid on plot area
vx = zeros(grid_density,grid_density); %x velocity at all grid points
vz = zeros(grid_density,grid_density); %z velocity at all grid points
v_abs = zeros(grid_density,grid_density); %absolute velocity at all grid points
x_grid = linspace(plot_limit(1),plot_limit(2),grid_density);
z_grid = linspace(plot_limit(3),plot_limit(4),grid_density);



%Calculating the velocity induced by the 
for i = 1:grid_density
    for j = 1:grid_density
        for k = 1:length(gamma)
            [vx(i,j),vz(i,j)] = VOR2D(gamma(k),x_grid(i),z_grid(j),gamma_x(k),gamma_z(k));
            vx(i,j) = vx(i,j) + U_inf;
            vz(i,j) = vz(i,j);
            vabs(i,j) = sqrt(vx(i,j)^2 + vz(i,j)^2);
        end
    end   
end

end
function [vx,vz] = v_distribution(U_inf,grid_density,plot_limit,gamma,af_geo)

gamma_x = [af_geo.VP_xU,af_geo.VP_xL];
gamma_z = [af_geo.VP_zU,af_geo.VP_zL];

%constructing grid on plot area
vx = ones(grid_density,grid_density); %x velocity at all grid points
vz = ones(grid_density,grid_density); %z velocity at all grid points
x_grid = linspace(plot_limit(1),plot_limit(2),grid_density);
z_grid = linspace(plot_limit(3),plot_limit(4),grid_density);



%Calculating the velocity induced by the 
for i = 1:grid_density
    for j = 1:grid_density
        for k = 1:n_panels
            [vx(i,j),vz(i,j)] = VOR2D(gamma(k),x_grid(i),z_grid(j),gamma_x(k),gamma_z(k));
            vx(i,j) = vx(i,j) + U_inf;
            vz(i,j) = vz(i,j);
        end
    end   
end

end
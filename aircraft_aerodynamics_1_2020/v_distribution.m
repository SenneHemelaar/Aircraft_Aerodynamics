function [fields] = v_distribution(grid_density, plot_limit, gamma, af_geo, u_0, w_0, rho)
%V_DISTRIBUTION Calculates the v_field field

%constructing grid on plot area
vx     = zeros(grid_density,grid_density) + u_0;
vz     = zeros(grid_density,grid_density) + w_0;
vabs   = zeros(grid_density,grid_density);
p      = zeros(grid_density,grid_density);
x_grid = linspace(plot_limit(1),plot_limit(2),grid_density);
z_grid = linspace(plot_limit(3),plot_limit(4),grid_density);

% Calculating induced 
for i = 1:grid_density
    for j = 1:grid_density
        for k = 1:length(gamma)
            [u,w] = VOR2D(x_grid(i), z_grid(j), af_geo.VP_x(k),...
                          af_geo.VP_z(k), gamma(k));
            vx(i,j) = vx(i,j) + u;
            vz(i,j) = vz(i,j) + w;
        end
        vabs(i,j) = sqrt(vx(i,j).^2 + vz(i,j).^2);
        p(i,j) = (0.5*rho*(u_0^2+w_0^2) -0.5*rho*vabs(i,j)^2)/(0.5*rho*(u_0^2+w_0^2)) ;

    end   
end

fields.vx   = vx;
fields.vz   = vz;
fields.vabs = vabs;
fields.p    = p;
end


function plot_velocity(grid_density, plot_limit, fields, af_geo)
%PLOT_VELOCITY Summary of this function goes here
%   Detailed explanation goes here
x_grid = linspace(plot_limit(1),plot_limit(2),grid_density);
z_grid = linspace(plot_limit(3),plot_limit(4),grid_density);

figure(1)
hold on; box on
cMap=jet(256);
[cb,h] = contourf(x_grid,z_grid,fields.vabs');
set(h, 'edgecolor','none');
colormap(cMap);
cb=colorbar;
cb.Label.String = 'velocity [m/s]';
xlabel('x');
ylabel('z');
plot(af_geo.x, af_geo.z, 'k-')

end


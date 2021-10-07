function plot_velocity(vx,vz,vabs,grid_density,plot_limit,af_geo)

x_grid = linspace(plot_limit(1),plot_limit(2),grid_density);
z_grid = linspace(plot_limit(3),plot_limit(4),grid_density);

figure(1)
hold on
contourf(x_grid,z_grid,vabs')


plot(af_geo.xU,af_geo.zU,'bo-')
hold on
plot(af_geo.xL,af_geo.zL,'bo-')
axis equal
plot(af_geo.CP_xU,af_geo.CP_zU,'go')
plot(af_geo.VP_xU,af_geo.VP_zU,'ro')
plot(af_geo.CP_xL,af_geo.CP_zL,'go')
plot(af_geo.VP_xL,af_geo.VP_zL,'ro')
plot(af_geo.xC,af_geo.zC,'k--')

end
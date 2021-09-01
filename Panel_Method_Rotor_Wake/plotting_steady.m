%% Plotting Results for Steady Panel Method 


figure(1)
plot(v_g.x_c,abs(CL(10)-0.2086),'-x')
hold on 
plot(v_g.x_c,abs(CL(15)-0.7786),'-o')
xlabel('x/c');
ylabel('C_L');
hold on;



A = readtable('naca2205_aseq_-6_90_6.txt');
cl_aoa = figure(2)
plot(alphas_deg, CL,'-x');
hold on;
plot(A{:,1},A{:,2});
xlabel('AoA');
ylabel('C_L');
xlim([-10 20]);
legend('Panel Method: Thin Airfoil with 2% camber, max camber at 0.2*x/c', 'XFoil Data: NACA 2205','Location','Best')
grid on; 
% saveas(cl_aoa, 'figures/steady_CL_vs_aoa.png');

figure(3)
cMap=jet(256); %set the colomap using the "jet" scale
[cb,h]=contourf(x_field,z_field,abs_vel_field);
set(h, 'edgecolor','none');
colormap(cMap);
cb=colorbar;
cb.Label.String = 'velocity [m/s]';
xlabel('x');
ylabel('z');
hold on;
plot(v_g.x, v_g.z,'-b');
box on;

figure(4)
cMap=jet(256); %set the colomap using the "jet" scale
[cb,h]=contourf(x_field,z_field,pres_field);
set(h, 'edgecolor','none');
colormap(cMap);
cb=colorbar;
cb.Label.String = '\Delta C_P';
xlabel('x');
ylabel('z');
hold on;
plot(v_g.x, v_g.z, '-b');
box on;


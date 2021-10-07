%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              MAIN: PANEL METHOD             %%%
%%% Authors: Jasper van Dongen & Senne Hemelaar %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all

%% Configurable parameters
n_panels = 40;
NACA = '0016';      % Choose desired 4 digit NACA airfoil NACA'INPUT'
U_inf = 5;
plot_limit = [-2 2 -2.5 2.5]; %x min x max y min y max
grid_density = 100;             %number of grid points on nxn matrix
rho = 1.23;          % [kg/m^3] density air


%% Generate Geometry 
af_geo = generate_geometry(n_panels,NACA);  % af_geo stands for airfoil geometry

%Test-plot to check geometry
figure(1)
plot(af_geo.xU,af_geo.zU,'bo-')
hold on
plot(af_geo.xL,af_geo.zL,'bo-')
axis equal
plot(af_geo.CP_xU,af_geo.CP_zU,'go')
plot(af_geo.VP_xU,af_geo.VP_zU,'ro')
plot(af_geo.CP_xL,af_geo.CP_zL,'go')
plot(af_geo.VP_xL,af_geo.VP_zL,'ro')
plot(af_geo.xC,af_geo.zC,'k--')

%% Find cirulation gamma
gamma = find_circulation(af_geo,n_panels,U_inf);

%% Calculating forces
[vx,vz,vabs] = v_distribution(U_inf,grid_density,plot_limit,gamma,af_geo);
L = loads(U_inf,gamma,rho,af_geo);

%% Plotting the results 
plot_velocity(vx,vz,vabs,grid_density,plot_limit,af_geo);


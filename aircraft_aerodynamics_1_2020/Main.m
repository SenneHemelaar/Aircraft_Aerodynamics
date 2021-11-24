%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%       MAIN: PANEL METHOD THIN AIRFOIL       %%%
%%% Authors: Jasper van Dongen & Senne Hemelaar %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all

%% Configurable parameters
U_inf = 10;                         % free stream velocity
alpha = 0;                          % angle of attack
u_0 = U_inf*cos(alpha);             % u component free stream velocity      
w_0 = U_inf*sin(alpha);             % w component free stream velocity
plot_limit = [-1 2 -1 1];           % x min x max y min y max
grid_density = 25;                 % number of grid points on nxn matrix
rho = 1.23;                         % density

%% Geometry Generator
N = 20;                             % number of panels
chord = 1;                          % mean aerodynamic chord                                      
camber = 0.02;                      % max camber
x_max_camber = 0.2;                 % location max camber
af_geo = generate_geometry(chord, N, camber, x_max_camber);

%% Find cirulation gamma
gamma = find_circulation(af_geo, N, u_0, w_0);

%% Calculating velocity field and loads
fields = v_distribution(grid_density, plot_limit, gamma, af_geo, u_0, w_0, rho);
[L, dp] = loads(gamma, rho, U_inf, chord, N);

%% Plotting the results
velocity  = true;
pres_dist = true;

if velocity
    plot_velocity(grid_density, plot_limit, fields, af_geo);
end
if pres_dist
    plot_pres_dist(chord, dp)
end


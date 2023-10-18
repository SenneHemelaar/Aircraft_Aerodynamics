%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  MAIN: PANEL METHOD THIN AIRFOIL                    %%%
%%%             Authors: Jasper van Dongen & Senne Hemelaar             %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all

% Configurable parameters
U_inf = 10;                         % free stream velocity
alpha_max = deg2rad(0);            % maximum angle of attack [rad]
alpha_min = deg2rad(0);             % minimum angle of attack [rad]
n_alpha = 30;                       % number of angles of attack to evaluate
alphas = linspace(alpha_min,...
         alpha_max, n_alpha);       % list of angles of attack
plot_limit = [-1 2 -1 1];           % x min x max y min y max for fields
grid_density = 25;                  % number of grid points on nxn matrix
rho = 1.23;                         % density

% Geometry Generator
N = 30;                             % number of panels
chord = 1;                          % mean aerodynamic chord                                      
camber = [0, 0.02, 0.05];           % max camber list
x_max_camber = 0.2;                 % location max camber
for i = 1:length(camber)
af_geo{i} = generate_geometry(chord, N, camber(i), x_max_camber);
end

for i = 1:length(alphas)            % loop over angles of attack
    for j = 1:length(camber)        % loop over max camber values
        alpha  = alphas(i);         % angle of attack
        u_0 = U_inf*cos(alpha);     % u component free stream velocity      
        w_0 = U_inf*sin(alpha);     % w component free stream velocity

        % Find cirulation gamma
        gamma{i,j} = find_circulation(af_geo{j}, N, u_0, w_0);

        % Calculating velocity field and loads
        fields{i,j} = v_distribution(grid_density, plot_limit,...
                      gamma{i,j}, af_geo{j}, u_0, w_0, rho);
        [CL(i,j), dp{i,j}, dCp{i,j}] = loads(gamma{i,j}, rho, U_inf,...
                                       chord, N);
    end
end

% Loading validation/ verification data
CL_ver = table2array(readtable('xf-n0009sm-il-50000'));
CL_val = load('0012.abbottdata.cl.dat');

CP_val_00 = textread('Pressure_Data_Alpha_00_Re_6mil.txt');
CP_val_10 = textread('Pressure_Data_Alpha_10_Re_6mil.txt');
CP_val_15 = textread('Pressure_Data_Alpha_15_Re_6mil.txt');

dCP_val_15 = (CP_val_15(24:end,2)) - flip(CP_val_15(1:22,2));
dCP_val_10 = (CP_val_10(24:end,2)) - flip(CP_val_10(1:22,2));

% From here all plotting procedures follow

figure(5)
plot(CP_val_15(:,1),CP_val_15(:,2))

%% Plotting the results

% % Lift polar comparison with open literature
% figure(1)
% box on; hold on; grid on;
% plot(rad2deg(alphas),CL(:,1),'k-')
% plot(CL_ver(:,1),CL_ver(:,2),'k--')
% plot(CL_val(:,1),CL_val(:,2),'k:','LineWidth',2)
% xlabel('$\alpha$','Interpreter','latex')
% ylabel('$C_L$','Interpreter','latex')
% title('')
% legend('Zero-Thickness Panel Method','NACA-0009 XFOIL','NACA-0012 Experimental Data','Interpreter','latex','Location','NorthWest')
% xlim([-20, 20])
% set(gcf,'position',[300,300,950,450])

% % Pressure distribution comparison with open literature
% x1 = linspace(0+0.5*(chord/20), chord-0.5*(chord/20), length(dCp{1}));
% x2 = linspace(0+0.5*(chord/20), chord-0.5*(chord/20), length(dCP_val_15));
% figure(2)
% hold on; grid on; box on;
% plot(x1, dCp{1},'k')
% plot(x2, -dCP_val_10,'k:','LineWidth',2)
% xlabel('$x/c$','Interpreter','latex')
% ylabel('$\Delta C_P$','Interpreter','latex')
% set(gcf,'position',[100,300,700,500])
% legend('Zero-Thickness Panel Method','NACA-0012 Experimental Data','Interpreter','latex')
% figure(3)
% hold on; grid on; box on;
% plot(x1, dCp{end},'k')
% plot(x2, -dCP_val_15,'k:','LineWidth',2)
% xlabel('$x/c$','Interpreter','latex')
% ylabel('$\Delta C_P$','Interpreter','latex')
% legend('Zero-Thickness Panel Method','NACA-0012 Experimental Data','Interpreter','latex')
% set(gcf,'position',[950,300,700,500])

% % Lift polar Varying camber
% figure(5)
% box on; hold on; grid on;
% plot(rad2deg(alphas),CL(:,3),'k:','LineWidth',2)
% plot(rad2deg(alphas),CL(:,2),'k--')
% plot(rad2deg(alphas),CL(:,1),'k-')
% xlabel('$\alpha$','Interpreter','latex')
% ylabel('$C_L$','Interpreter','latex')
% set(gcf,'position',[300,300,950,450])
% legend('Max Camber = 0.00','Max Camber = 0.02','Max Camber = 0.05','Interpreter','latex','Location','NorthWest')

% % Pressure distribution Varying camber
% x1 = linspace(0+0.5*(chord/20), chord-0.5*(chord/20), length(dCp{1}));
% figure(6)
% hold on; grid on; box on;
% % plot(x1, dCp{1,1},'k')
% plot(x1, dCp{1,2},'k--')
% plot(x1, dCp{1,3},'k:','LineWidth',2)
% xlabel('$x/c$','Interpreter','latex')
% ylabel('$\Delta C_P$','Interpreter','latex')
% set(gcf,'position',[100,300,700,500])
% legend('Max Camber = 0.00','Max Camber = 0.02','Max Camber = 0.05','Interpreter','latex')
% figure(7)
% hold on; grid on; box on;
% % plot(x1, dCp{end,1},'k')
% plot(x1, dCp{end,2},'k--')
% plot(x1, dCp{end,3},'k:','LineWidth',2)
% xlabel('$x/c$','Interpreter','latex')
% ylabel('$\Delta C_P$','Interpreter','latex')
% set(gcf,'position',[950,300,700,500])
% legend('Max Camber = 0.00','Max Camber = 0.02','Max Camber = 0.05','Interpreter','latex')


velocity   = true;
if velocity
    plot_velocity(grid_density, plot_limit, fields{1,3}, af_geo{3});
end

%}

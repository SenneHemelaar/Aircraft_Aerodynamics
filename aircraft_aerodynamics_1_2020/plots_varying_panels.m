close all; clear; clc;

alpha_max = deg2rad(-15);            % maximum angle of attack [rad]
alpha_min = deg2rad(15);           % minimum angle of attack [rad]
n_alpha = 30;                       % number of angles of attack to evaluate
alphas = linspace(alpha_min,...
         alpha_max, n_alpha);       % list of angles of attack
CL_10  = load('CL_10_panels.mat');
CL_30  = load('CL_30_panels.mat');
dCp_10 = load('dCp_10_panels.mat');
dCp_30 = load('dCp_30_panels.mat');

% figure(1)
% hold on; grid on; box on;
% plot(rad2deg(alphas),CL_10.CL(:,3),'k-')
% plot(rad2deg(alphas),CL_30.CL(:,3),'k:','LineWidth',1.5)
% xlim([-5 5])
% xlabel('$\alpha$','Interpreter','latex')
% ylabel('$C_L$','Interpreter','latex')
% legend('10 Panels','30 Panels','Interpreter','latex','Location','NorthWest')
% set(gcf,'position',[300,300,950,450])

% x1 = linspace(0+0.5*(1/20), 1-0.5*(1/20), 10);
% x2 = linspace(0+0.5*(1/20), 1-0.5*(1/20), 30);
% figure(2)
% hold on; grid on; box on;
% plot(x1, dCp_10.dCp{1,3},'k-')
% plot(x2, dCp_30.dCp{1,3},'k:','LineWidth',1.5)
% xlabel('$x$','Interpreter','latex')
% ylabel('$\Delta C_p$','Interpreter','latex')
% legend('10 Panels','30 Panels','Interpreter','latex','Location','NorthEast')
% set(gcf,'position',[100,300,700,500])
% figure(3)
% hold on; grid on; box on;
% plot(x1, dCp_10.dCp{end,3},'k-')
% plot(x2, dCp_30.dCp{end,3},'k:','LineWidth',1.5)
% xlabel('$x$','Interpreter','latex')
% ylabel('$\Delta C_p$','Interpreter','latex')
% legend('10 Panels','30 Panels','Interpreter','latex','Location','NorthEast')
% set(gcf,'position',[950,300,700,500])


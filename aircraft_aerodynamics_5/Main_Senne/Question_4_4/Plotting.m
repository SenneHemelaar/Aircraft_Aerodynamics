clc; close all; clear;
C = linspecer(5);

% select only one of the settings to be true
settings.varying_N = true;
settings.varying_D = true;

if settings.varying_D
% Varying D, Constant CT
results_D_CT_015 = load('results_D_CT_0.15.mat'); results_D_CT_015 = results_D_CT_015.results;
results_D_CT_01 = load('results_D_CT_0.1.mat'); results_D_CT_01 = results_D_CT_01.results;
results_D_CT_005 = load('results_D_CT_0.05.mat'); results_D_CT_005 = results_D_CT_005.results;

f1 = figure(1);
box on; hold on; grid on;
f1.Position = [200 200 700 450];
plot(results_D_CT_015.D, results_D_CT_015.eff,'Color',C(1,:),'linewidth',2)
plot(results_D_CT_01.D, results_D_CT_01.eff-0.005,'Color',C(2,:),'linewidth',2)
plot(results_D_CT_005.D, results_D_CT_005.eff,'Color',C(3,:),'linewidth',2)
xlabel('$D$ [m]','interpreter','latex','fontsize',13)
ylabel('$\eta$','interpreter','latex','fontsize',13)
ylim([0.3 0.9])
legend('$C_T=0.15$','$C_T=0.1$','$C_T=0.05$','interpreter','latex','fontsize',12,'location','southwest')

% Varying D, Constant CP
results_D_CP_02 = load('results_D_CP_0.2.mat'); results_D_CP_02 = results_D_CP_02.results;
results_D_CP_015 = load('results_D_CP_0.15.mat'); results_D_CP_015 = results_D_CP_015.results;
results_D_CP_01 = load('results_D_CP_0.1.mat'); results_D_CP_01 = results_D_CP_01.results;

f2 = figure(2);
box on; hold on; grid on;
f2.Position = [1000 200 700 450];
plot(results_D_CP_02.D, results_D_CP_02.eff,'Color',C(1,:),'linewidth',1.8)
plot(results_D_CP_015.D, results_D_CP_015.eff,'Color',C(2,:),'linewidth',1.8)
plot(results_D_CP_01.D, results_D_CP_01.eff,'Color',C(3,:),'linewidth',1.8)
xlabel('$D$ [m]','interpreter','latex','fontsize',13)
ylabel('$\eta$','interpreter','latex','fontsize',13)
ylim([0.3 0.9])
legend('$C_P=0.2$','$C_P=0.15$','$C_P=0.1$','interpreter','latex','fontsize',12,'location','southwest')
end

if settings.varying_N
% Varying N, Constant CT
results_N_CT_015 = load('results_N_CT_0.15.mat'); results_N_CT_015 = results_N_CT_015.results;
results_N_CT_01 = load('results_N_CT_0.1.mat'); results_N_CT_01 = results_N_CT_01.results;
results_N_CT_005 = load('results_N_CT_0.05.mat'); results_N_CT_005 = results_N_CT_005.results;

f1 = figure(3);
box on; hold on; grid on;
f1.Position = [200 200 700 450];
plot(results_N_CT_015.N(2:end), results_N_CT_015.eff(2:end),'s-','Color',C(1,:),'linewidth',2)
plot(results_N_CT_01.N(2:end), results_N_CT_01.eff(2:end),'s-','Color',C(2,:),'linewidth',2)
plot(results_N_CT_005.N, results_N_CT_005.eff,'s-','Color',C(3,:),'linewidth',2)
xlim([2 10]); ylim([0.3 0.9])
xlabel('$N$','interpreter','latex','fontsize',13)
ylabel('$\eta$','interpreter','latex','fontsize',13)
legend('$C_T=0.15$','$C_T=0.1$','$C_T=0.05$','interpreter','latex','fontsize',12,'location','southwest')


% Varying N, Constant CP
results_N_CP_02 = load('results_N_CP_0.2.mat'); results_N_CP_02 = results_N_CP_02.results;
results_N_CP_015 = load('results_N_CP_0.15.mat'); results_N_CP_015 = results_N_CP_015.results;
results_N_CP_01 = load('results_N_CP_0.1.mat'); results_N_CP_01 = results_N_CP_01.results;

f2 = figure(4);
box on; hold on; grid on;
f2.Position = [1000 200 700 450];
plot(results_N_CP_02.N(2:end), results_N_CP_02.eff(2:end),'s-','Color',C(1,:),'linewidth',1.8)
plot(results_N_CP_015.N(2:end), results_N_CP_015.eff(2:end),'s-','Color',C(2,:),'linewidth',1.8)
plot(results_N_CP_01.N(2:end), results_N_CP_01.eff(2:end),'s-','Color',C(3,:),'linewidth',1.8)
xlim([2 10]); ylim([0.3 0.9])
xlabel('$N$','interpreter','latex','fontsize',13)
ylabel('$\eta$','interpreter','latex','fontsize',13)
legend('$C_P=0.2$','$C_P=0.15$','$C_P=0.1$','interpreter','latex','fontsize',12,'location','southwest')
end


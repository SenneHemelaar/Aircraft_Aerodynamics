clc; close all; clear;
C = [108/256, 194/256, 74/256; ...
     237/256, 104/256, 66/256; ...
     0, 166/256, 214/256;      ...
     224/256, 60/256, 49/256];
Thrust_settings = [500,600,1200];
Power_settings = [18000,20000,40000];

D_T_str = ["results_D_T_",'results_D_T_','results_D_T_'];
D_P_str = ["results_D_P_",'results_D_P_','results_D_P_'];
N_T_str = ["results_N_T_",'results_N_T_','results_N_T_'];
N_P_str = ["results_N_P_",'results_N_P_','results_N_P_'];

for i = 1:length(Thrust_settings)
    D_T_str(i) = D_T_str(i) + num2str(Thrust_settings(i)) + ".mat";
    D_P_str(i) = D_P_str(i) + num2str(Power_settings(i)) + ".mat";
    N_T_str(i) = N_T_str(i) + num2str(Thrust_settings(i)) + ".mat";
    N_P_str(i) = N_P_str(i) + num2str(Power_settings(i)) + ".mat";
end




% select only one of the settings to be true
settings.varying_N = true;
settings.varying_D = true;



% Varying D, Constant T
results_D_T_500 = load(D_T_str(1)); results_D_T_500 = results_D_T_500.results;
results_D_T_600 = load(D_T_str(2)); results_D_T_600 = results_D_T_600.results;
results_D_T_700 = load(D_T_str(3)); results_D_T_700 = results_D_T_700.results;

results_D_T_500.D(results_D_T_500.D == 0) = NaN;
results_D_T_600.D(results_D_T_600.D == 0) = NaN;
results_D_T_700.D(results_D_T_700.D == 0) = NaN;

f1 = figure(1);
box on; hold on; grid minor;
f1.Position = [200 200 700 450];
plot(results_D_T_500.D, results_D_T_500.eff,'o-','Color',C(1,:),'linewidth',2)
plot(results_D_T_600.D, results_D_T_600.eff,'o-','Color',C(2,:),'linewidth',2)
plot(results_D_T_700.D, results_D_T_700.eff,'o-','Color',C(3,:),'linewidth',2)
xlabel('$D$ [m]','interpreter','latex','fontsize',13)
ylabel('$\eta$ [-]','interpreter','latex','fontsize',13)
% ylim([0.3 0.9])
legend('$T=500 [N]$','$T=600 [N]$','$T=700 [N]$','interpreter','latex','fontsize',12,'location','northwest')

% Varying D, Constant CP
results_D_P_20000 = load(D_P_str(1)); results_D_P_20000 = results_D_P_20000.results;
results_D_P_22000 = load(D_P_str(2)); results_D_P_22000 = results_D_P_22000.results;
results_D_P_24000 = load(D_P_str(3)); results_D_P_24000 = results_D_P_24000.results;

results_D_P_20000.D(results_D_P_20000.D == 0) = NaN;
results_D_P_22000.D(results_D_P_22000.D == 0) = NaN;
results_D_P_24000.D(results_D_P_24000.D == 0) = NaN;

f2 = figure(2);
box on; hold on; grid minor;
f2.Position = [1000 200 700 450];
plot(results_D_P_20000.D, results_D_P_20000.eff,'o-','Color',C(1,:),'linewidth',1.8)
plot(results_D_P_22000.D, results_D_P_22000.eff,'o-','Color',C(2,:),'linewidth',1.8)
plot(results_D_P_24000.D, results_D_P_24000.eff,'o-','Color',C(3,:),'linewidth',1.8)
xlabel('$D$ [m]','interpreter','latex','fontsize',13)
ylabel('$\eta$ [-]','interpreter','latex','fontsize',13)

legend('$P=18000 [W]$','$P=20000 [W]$','$P=22000 [W]$','interpreter','latex','fontsize',12,'location','south')


% Varying N, Constant CT
N_T_1 = load(N_T_str(1)); N_T_1  = N_T_1.results;
N_T_2 = load(N_T_str(2)); N_T_2 = N_T_2.results;
N_T_3 = load(N_T_str(3)); N_T_3 = N_T_3.results;

N_T_1.N(N_T_1.N == 0) = NaN;
N_T_2.N(N_T_2.N == 0) = NaN;
N_T_3.N(N_T_3.N == 0) = NaN;

f3 = figure(3);
box on; hold on; grid minor;
f3.Position = [200 200 700 450];
plot(N_T_1 .N, N_T_1 .eff,'o-','Color',C(1,:),'linewidth',2)
plot(N_T_2.N, N_T_2.eff,'o-','Color',C(2,:),'linewidth',2)
plot(N_T_3.N, N_T_3.eff,'o-','Color',C(3,:),'linewidth',2)
% xlim([2 10]); ylim([0.3 0.9])
xlabel('$N$ [-]','interpreter','latex','fontsize',13)
ylabel('$\eta$ [-]','interpreter','latex','fontsize',13)
legend('$T=500 [N]$','$T=600 [N]$','$T=700 [N]$','interpreter','latex','fontsize',12,'location','southwest')


% Varying N, Constant CP
N_P_1 = load(N_P_str(1)); N_P_1  = N_P_1.results;
N_P_2 = load(N_P_str(2)); N_P_2 = N_P_2.results;
N_P_3 = load(N_P_str(3)); N_P_3 = N_P_3.results;

N_P_1.N(N_P_1.N == 0) = NaN;
N_P_2.N(N_P_2.N == 0) = NaN;
N_P_3.N(N_P_3.N == 0) = NaN;

f4 = figure(4);
box on; hold on; grid minor;
f4.Position = [1000 200 700 450];
plot(N_P_1.N, N_P_1.eff,'o-','Color',C(1,:),'linewidth',1.8)
plot(N_P_2.N, N_P_2.eff,'o-','Color',C(2,:),'linewidth',1.8)
plot(N_P_3.N, N_P_3.eff,'o-','Color',C(3,:),'linewidth',1.8)
% xlim([2 10]); ylim([0.3 0.9])
xlabel('$N$ [-]','interpreter','latex','fontsize',13)
ylabel('$\eta$ [-]','interpreter','latex','fontsize',13)
legend('$P=18000 [W]$','$P=20000 [W]$','$P=22000 [W]$','interpreter','latex','fontsize',12,'location','southwest')


saveas(f1,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\D_vs_eta_constCT_varCP.png")
saveas(f2,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\D_vs_eta_varCP_constCT.png") 
saveas(f3,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\N_vs_eta_constCT_varCP.png")
saveas(f4,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\N_vs_eta_varCP_constCT.png") 
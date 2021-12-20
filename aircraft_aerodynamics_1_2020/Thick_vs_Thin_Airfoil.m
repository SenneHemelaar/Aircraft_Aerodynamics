clc; clear; close all;

Data = table2array(readtable('Assignment_Data'));

dCp_Data_05 = (Data(81:end,5)) - flip(Data(1:80,5));
dCp_05 = load('dCp_Thin_Symmetrical_Airfoil_AoA_5.mat');

x2 = linspace(0+0.5*(1/30), 1-0.5*(1/30), 30);
figure(1)
hold on; grid on; box on;
plot(x2, dCp_05.dCp{1,1},'k-')
plot(flip(Data(1:80,2)),dCp_Data_05,'k:','LineWidth',1.5)
xlabel('$x/c$','Interpreter','latex')
ylabel('$\Delta C_p$','Interpreter','latex')
legend('Thin Airfoil','NACA-0015','Interpreter','latex','Location','NorthEast')
set(gcf,'position',[300,300,950,450])

CL = (1/(length(dCp_Data_05))) * sum(dCp_Data_05);
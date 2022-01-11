clear all; close all

Forces = load('forces_varying_phiw'); Forces = Forces.Forces;
phiw_list = linspace(0,90,91);

Forces_list = zeros(91,1);
for i = 1:length(phiw_list)
    Forces_list(i) = Forces{i}.CDind;
end
    


% plotting
figure(1)
hold on; grid on; box on
plot(phiw_list,Forces_list,'k')
xlabel('$\phi_w$','Interpreter','latex')
ylabel('$C_D$','Interpreter','latex')
set(gcf,'position',[300,300,950,450])
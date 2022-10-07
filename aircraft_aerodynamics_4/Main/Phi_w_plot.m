clear; close all

Forces = load('forces_varying_phiw'); Forces = Forces.Forces;
phiw_list = linspace(0,90,91);

Forces_list = zeros(91,1);
for i = 1:length(phiw_list)
    Forces_list(i) = Forces{i}.CDind;
end

figure(1)
hold on; grid on; box on
plot(phiw_list,Forces_list,'k','LineWidth',1.5)
xlabel('$\phi_w $ [deg]','Interpreter','latex','FontSize',16)
ylabel('${C_D}_i$','Interpreter','latex','FontSize',15)
set(gcf,'position',[300,300,950,450])
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
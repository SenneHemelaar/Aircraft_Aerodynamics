clear; close all

M = load('moments_varying_phiw');
M = (M.M)*(10^-6);
phiw_list = linspace(0,90,91);



figure(1)
hold on; grid on; box on
plot(phiw_list,M,'k','LineWidth',1.5)
xlabel('$\phi_w $ [deg]','Interpreter','latex','FontSize',16)
ylabel('$M_r$  [MNm]','Interpreter','latex','FontSize',14)
set(gcf,'position',[300,300,950,450])
% ylh = get(gca,'ylabel');
% gyl = get(ylh);                                                         % Object Information
% ylp = get(ylh, 'Position');
% set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
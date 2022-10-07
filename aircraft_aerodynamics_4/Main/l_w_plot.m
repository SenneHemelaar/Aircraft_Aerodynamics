clear; close all

Forces = load('forces_varying_lw'); Forces = Forces.Forces;
lw_list = linspace(0.03,0.09,100);

Forces_list = zeros(91,1);
for i = 1:length(lw_list)
    Forces_list(i) = Forces{i}.CDind;
end

Cd_sp = [0.0151 0.0149 0.0144 0.0148 0.0145 0.0145 0.0144 0.0143 0.0142 0.0141 0.0140 0.0137 0.0139 0.0137 0.0137 0.0136 0.0137 0.0134 0.0133 0.0132 0.0132 0.0130 0.0129 0.0129 0.0128 0.0127];
span_list = linspace(0.03, 0.09, length(Cd_sp));

figure(1)
hold on; grid on; box on
plot(lw_list,Forces_list,'k','LineWidth',1.5)
plot(span_list,Cd_sp,'k--','LineWidth',1.5)
xlabel('$\%\,b$','Interpreter','latex','FontSize',16)
ylabel('${C_D}_i$','Interpreter','latex','FontSize',15)
set(gcf,'position',[300,300,950,450])
legend('Winglet Height','Span Expansion','interpreter','latex','FontSize',12)
ylh = get(gca,'ylabel');
gyl = get(ylh);                                                         % Object Information
ylp = get(ylh, 'Position');
set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
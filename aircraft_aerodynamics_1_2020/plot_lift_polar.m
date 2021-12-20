function plot_lift_polar(alphas, CL, CL_ver, CL_val)
%PLOT_LIFT_POLAR

figure(3)
box on; hold on; grid on;
plot(rad2deg(alphas),CL,'k-')
plot(CL_ver(:,1),CL_ver(:,2),'k--')
plot(CL_val(:,1),CL_val(:,2),'k:','LineWidth',2)
xlabel('$\alpha$','Interpreter','latex')
ylabel('$C_L$','Interpreter','latex')
title('')
legend('Zero-Thickness Panel Method','NACA-0009 XFOIL','NACA-0012 Experimental','Interpreter','latex','Location','NorthWest')
xlim([-20, 20])
set(gcf,'position',[300,300,700,500])

end


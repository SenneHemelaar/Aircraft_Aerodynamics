%% Plotting Unsteady
C = linspecer(6);

% figure(5)
% plot(v_g.x_c,Dp)
% xlabel('x/c');
% ylabel('\Delta C_p');

figure(8)
plot(rad2deg(alpha_list(N-time_period:N-1)), period_lift)
xlabel('\alpha')
ylabel('C_L')
grid on;
hold on;

%error figure, should be fixed
% figure
% semilogx(ceil(1:N/time_period), error_list)
% xlabel('N')
% ylabel('error')
% grid on;
 
line_types = ["-.", "--", "-" ];
line_type = line_types(dt_step);

figure(3)
loglog(ceil(1:N/time_period), error_list, line_type, 'color', [0,0,0] + 0.25*(dt_step-1))
xlabel('N')
ylabel('error')
ylim([1e-4, 1])
grid on;
hold on;


%% logic to find the right indices for the angles of attack
% converged_alphas = alpha_list(end-time_period:end);
% [max_alpha, idx_max_alpha] = max(converged_alphas);
% [min_alpha, idx_min_alpha] = min(converged_alphas);
% 
% % max_alpha_idx_wr_end = -time_period+idx_max_alpha;
% % min_alpha_idx_wr_end = -time_period+idx_min_alpha;
% 
% len = length(converged_alphas);
% low_bound_idx = idx_min_alpha  - ceil(len/4);
% high_bound_idx = idx_max_alpha  - ceil(len/4);
% 
% %when index is negative start from the end and substract it
% if low_bound_idx < 0
%     low_bound_idx = len + low_bound_idx;
% elseif high_bound_idx < 0
%     high_bound_idx = len + high_bound_idx;
% end
% 
% low_bound_alpha = converged_alphas(low_bound_idx);
% high_bound_alpha = converged_alphas(high_bound_idx);
% 
% indices = [idx_max_alpha, idx_min_alpha, high_bound_idx, low_bound_idx];
% aoa_s = [max_alpha, min_alpha, high_bound_alpha, low_bound_alpha];
% 
% %% heat plot
% for i=1:4
%     index = indices(i);
%     aoa = aoa_s(i);
%     figure()
%     
%     vel = squeeze(abs_vel_field(:,:,index));
% %     pres = squeeze(pres_field(:,:,index));
% %     ax1 = subplot(2,2,i);
%     cMap=jet(256); %set the colomap using the "jet" scale
%     [cb,h]=contourf(x_field,z_field,vel, 100);
%     set(h, 'edgecolor','none');
%     colormap(cMap);
%     caxis([8 13])
%     cb=colorbar;
%     cb.Label.String = 'velocity [m/s]';
%     xlabel('x');
%     ylabel('z');
%     hold on;
%     plot(v_g.x, v_g.z,'-k');
%     box on;
% 
%     ax2 = subplot(2,2,2+j);
%     cMap=jet(256); %set the colomap using the "jet" scale
%     [cb,h]=contourf(x_field,z_field,pres);
%     set(h, 'edgecolor','none');
%     colormap(cMap);
%     caxis([-0.3 0.2])
%     cb=colorbar;
%     cb.Label.String = '\Delta C_P';
%     xlabel('x');
%     ylabel('z');
%     title(append('\alpha = ', string(round(rad2deg(aoa)))))
%     hold on;
%     plot(v_g.x, v_g.z, '-b');
%     box on;
% end

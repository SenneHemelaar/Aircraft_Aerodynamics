function plot_pres_dist(chord, dp, dCp)
%PLOT_PRES_DIST
x = linspace(0,chord,length(dp));

figure(2)
hold on; grid on; 
plot(x, dp)
end


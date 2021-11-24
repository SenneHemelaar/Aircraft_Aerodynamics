function plot_pres_dist(chord, dP)
%PLOT_PRES_DIST
x = linspace(0,chord,length(dP));

figure(2)
plot(x, dP)

end


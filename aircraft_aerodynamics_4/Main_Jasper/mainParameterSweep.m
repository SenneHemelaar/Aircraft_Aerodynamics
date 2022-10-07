%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        AVL parameter sweep                          %%%
%%%                     Author: Jasper van Dongen                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

addpath(genpath('runavl'))

blue = [114 147 203]./255;
red = [211 94 96]./255;
steps = 10;

rho = 0.6601;
a = 316.5;
mach = 0.7;

q = 0.5*(a*mach)^2*rho;
k = 0;
% Phi sweep
%%%===========================DesignVector==============================%%%
%     l_w , phi        ,C_w_r   ,lambda_w, Lambda       , e_w_r, e_w_t 
x = [1.68, deg2rad(0), 1.6    , 0.4    , deg2rad(30.0), 0     , 0];


%%%=====================Select parameter to sweep=======================%%%

param_num = 2;
param_lb = deg2rad(0);
param_up = deg2rad(90);


param = linspace(param_lb,param_up,steps);

phi.CDind = zeros(1, steps);
phi.CDff = zeros(1, steps);
phi.M = zeros(1, steps);
phi.x = zeros(1, steps);

for i = 1:steps
    x(param_num) = param(i);
    [~,forces] = runobj(x);
    phi.M(i) = optimization_function(forces, q, k)/1e6;
    phi.CDind(i) = forces.CDind;
    phi.CDff(i) = forces.CDff;
    phi.x(i) = param(i);
end

% l_w @ 0 phi sweep
%%%===========================DesignVector==============================%%%
%     l_w , phi        ,C_w_r   ,lambda_w, Lambda       , e_w_r, e_w_t 
x = [1.68, deg2rad(0), 1.6    , 0.4    , deg2rad(30.0), 0     , 0];


%%%=====================Select parameter to sweep=======================%%%

param_num = 1;
param_lb = 0.1;
param_up = 2;


param = linspace(param_lb,param_up,steps);

l_w.CDind = zeros(1, steps);
l_w.CDff = zeros(1, steps);
l_w.x = zeros(1, steps);

for i = 1:steps
    x(param_num) = param(i);
    [~,forces] = runobj(x);
    l_w.CDind(i) = forces.CDind;
    l_w.CDff(i) = forces.CDff;
    l_w.x(i) = param(i);
end

% l_w @ 90 phi sweep
%%%===========================DesignVector==============================%%%
%     l_w , phi        ,C_w_r   ,lambda_w, Lambda       , e_w_r, e_w_t 
x = [1.68, deg2rad(90), 1.6    , 0.4    , deg2rad(30.0), 0     , 0];


%%%=====================Select parameter to sweep=======================%%%

param_num = 1;
param_lb = 0.1;
param_up = 2;


param = linspace(param_lb,param_up,steps);

tip_ext.CDind = zeros(1, steps);
tip_ext.CDff = zeros(1, steps);
tip_ext.x = zeros(1, steps);

for i = 1:steps
    x(param_num) = param(i);
    [~,forces] = runobj(x);
    tip_ext.CDind(i) = forces.CDind;
    tip_ext.CDff(i) = forces.CDff;
    tip_ext.x(i) = param(i);
end


%%%================================plot=================================%%%

f1 = figure(1);
hold on 
grid on 

plot(rad2deg(phi.x),phi.CDind, 'Color', blue, 'LineWidth', 2)

ylim([min(phi.CDind)*0.95 max(phi.CDind)*1.05])
xlim([min(rad2deg(phi.x)) max(rad2deg(phi.x))])

ylabel('Induced drag C_{D_i} [-]')
xlabel('Cant angle \phi [deg]') 

f2 = figure(2);
hold on 
grid on 

plot(l_w.x,l_w.CDind,'Color', blue, 'LineWidth', 2)
plot(tip_ext.x,tip_ext.CDind, 'Color', red, 'LineWidth', 2)

ylim([min(min(l_w.CDind,tip_ext.CDind))*0.95 ,max(max(phi.CDind,tip_ext.CDind))*1.05])
xlim([min(l_w.x) max(l_w.x)])

ylabel('Induced drag C_{D_i} [-]')
xlabel('Tip extension l_w [m]')

legend('Winglet @ \phi = 0 deg','Winglet @ \phi = 90 deg')

f3 = figure(3);
grid on
hold on 
plot(rad2deg(phi.x),phi.M, 'Color', blue, 'LineWidth', 2)

ylim([min(phi.M)*0.95 max(phi.M)*1.05])
xlim([min(rad2deg(phi.x)) max(rad2deg(phi.x))])

ylabel('Root bending moment [MNm]')
xlabel('Cant angle \phi [deg]') 

%%%=======================saving figures================================%%%


saveas(f1,'cant_angles_Cdi.png')
saveas(f2,'winglet_height_tip_extension_Cdi.png')
saveas(f3,'cant_angles_bending_moment.png')
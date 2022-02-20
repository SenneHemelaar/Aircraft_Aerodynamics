%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    AVL optimisation in Matlab                       %%%
%%%                     Author: Jasper van Dongen                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

addpath(genpath('runavl'))
%%%=============================Diagnostics=============================%%%

%%%=======================Lower and upperbounds=========================%%%

lb = [0.56, deg2rad(10), 0.8, 0.4, deg2rad(0), deg2rad(-10) , deg2rad(-6)];
ub = [2.6, deg2rad(90), 2.0, 1.0, deg2rad(45), deg2rad(6) , deg2rad(+10)];

%%%===========================DesignVector==============================%%%

x0 = [1.58, 50, 1.4, 0.7,22.5, 0, 0];

%%%==============================Options================================%%%

options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.Diagnostics     = 'on';
options.PlotFcns        = {@optimplotx, ...                 % Plot functions
                           @optimplotfval, ...
                           @optimplotconstrviolation, ... 
                           @optimplotfirstorderopt, ...
                           @optimplotstepsize, ...
                           @optimplotfunccount};
                       
%%%=========================RunOptimisation=============================%%%
tic
x = fmincon(@runobj,x0,[],[],[],[],lb,ub,[],options);
toc
%%%=====================ExtractCovergenceHistory========================%%%







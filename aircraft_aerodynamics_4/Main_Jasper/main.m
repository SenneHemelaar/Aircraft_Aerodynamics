%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    AVL optimisation in Matlab                       %%%
%%%                     Author: Jasper van Dongen                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

addpath(genpath('runavl'))

%%%=======================Lower and upperbounds=========================%%%
lb = [0.56, deg2rad(10), 0.8, 0.4, deg2rad(0), deg2rad(-10) , deg2rad(-6)];
ub = [2.6, deg2rad(90), 2.0, 1.0, deg2rad(45), deg2rad(6) , deg2rad(+10)];

%%%===========================DesignVector==============================%%%
x0 = [1.58, deg2rad(50), 1.4, 0.7,deg2rad(22.5), 0, 0];

%%%==============================Options================================%%%
options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.Diagnostics     = 'on';
% options.PlotFcns        = {@optimplotx, ...                 % Plot functions
%                            @optimplotfval, ...
% %                            @optimplotconstrviolation, ... 
%                            @optimplotfirstorderopt, ...
%                            @optimplotstepsize, ...
%                            @optimplotfunccount};
options.DiffMinChange   = 0.01;
                       
%%%=========================RunOptimisation=============================%%%
tic
x = fmincon(@runobj,x0,[],[],[],[],lb,ub,[],options);
toc
%% =====================ExtractCovergenceHistory========================%%%

k_0 = load('x_vector_k_0.mat');   k_0 = k_0.x;
k_05= load('x_vector_k_0.5.mat'); k_05= k_05.x;
k_1 = load('x_vector_k_1.mat');   k_1 = k_1.x;





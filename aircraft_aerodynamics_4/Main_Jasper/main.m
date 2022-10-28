%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                    AVL optimisation in Matlab                       %%%
%%%                     Author: Jasper van Dongen                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

addpath(genpath('runavl'))

%%%==========================INPUTS GUIDE===============================%%%
% lw      : Winglet Length                          [% span]
% Lam     : Winglet Leading Edge Sweep              [degrees]
% phiw    : Winglet Cant Angle                      [degrees]
% cwr     : Winglet Root Chord                      [% main wing chord]
% lam     : Winglet Taper Ratio                     [cwt/cwr]
% epsR    : Winglet Root Twist Angle                [degrees]
% epsT    : Winglet Tip Twist Angle                 [degrees]
% cR      : Main Wing Root Chord                    [m]
% cT      : Main Wing Tip Chord                     [m]
% L       : Main Wing Half-Span                     [m]
% Ncm     : No. of Chordwise Sections Main Wing     [-]
% Ncw     : No. of Chordwise Sections Winglet       [-]
% Nsm     : No. of Spanwise Sections Main Wing      [-]
% Nsw     : No. of Spanwise Sections Winglet        [-]

%%%=============================INPUTS==================================%%%
lw_lb   = 0.02;  lw_ub   = 0.1;    lw_0   = 0.06;
phiw_lb = 10;    phiw_ub = 90;     phiw_0 = 45;
cwr_lb  = 0.4;   cwr_ub  = 1.0;    cwr_0  = 0.8;
lam_lb  = 0.4;   lam_ub  = 1.0;    lam_0  = 0.4;
Lam_lb  = 0;     Lam_ub  = 45;     Lam_0  = 30;
epsR_lb = -6;    epsR_ub = 10;     epsR_0 = 0;
epsT_lb = -10;   epsT_ub = 6;      epsT_0 = 0;
cR      = 5.5;
cT      = 2.0;
L       = 14;

%%%====================PRELIMINARY CALCULATIONS=========================%%%
lw_lb   = lw_lb*2*L;         lw_ub   = lw_ub*2*L;         lw_0   = lw_0*2*L;
Lam_lb  = deg2rad(Lam_lb);   Lam_ub  = deg2rad(Lam_ub);   Lam_0  = deg2rad(Lam_0);
phiw_lb = deg2rad(phiw_lb);  phiw_ub = deg2rad(phiw_ub);  phiw_0 = deg2rad(phiw_0);
cwr_lb  = cwr_lb * cT;       cwr_ub  = cwr_ub * cT;       cwr_0  = cwr_0 * cT;
epsR_lb = deg2rad(epsR_lb);  epsR_ub = deg2rad(epsR_ub);  epsR_0 = deg2rad(epsR_0);
epsT_lb = deg2rad(epsT_lb);  epsT_ub = deg2rad(epsT_ub);  epsT_0 = deg2rad(epsT_0);

%%%=======================lower and upper bounds========================%%%
lb = [lw_lb, phiw_lb, cwr_lb, lam_lb, Lam_lb, epsR_lb, epsT_lb];
ub = [lw_ub, phiw_ub, cwr_ub, lam_ub, Lam_ub, epsR_ub, epsT_ub];

%%%===========================DesignVector==============================%%%
x0 = [lw_0,  phiw_0,  cwr_0,  lam_0,  Lam_0,  epsR_0,  epsT_0];


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
options.DiffMinChange   = 0.01;
                       
%%%=========================RunOptimisation=============================%%%
tic
x = fmincon(@runobj,x0,[],[],[],[],lb,ub,[],options);
toc
%%%=====================ExtractCovergenceHistory========================%%%

% k_0 = load('x_vector_k_0.mat');   k_0 = k_0.x;
% k_05= load('x_vector_k_0.5.mat'); k_05= k_05.x;
% k_1 = load('x_vector_k_1.mat');   k_1 = k_1.x;

%% 

Geometry_plot_J(x)


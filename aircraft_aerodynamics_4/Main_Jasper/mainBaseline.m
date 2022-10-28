%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        AVL parameter sweep                          %%%
%%%                     Author: Jasper van Dongen                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

addpath(genpath('runavl'))

blue = [114 147 203]./255;
red = [211 94 96]./255;

rho = 0.6601;
a = 316.5;
mach = 0.7;

q = 0.5*(a*mach)^2*rho;
k = 0.5;

%%%===========================Baseline wing=============================%%%
%     l_w , phi        ,C_w_r   ,lambda_w, Lambda       , e_w_r, e_w_t 
x = [1.68, deg2rad(0), 1.6    , 0.4    , deg2rad(30.0), -45     , -6];

load('k_1_20102022.mat')

[J0,forces,M0] = runobj(x);
CDind0 = forces.CDind;


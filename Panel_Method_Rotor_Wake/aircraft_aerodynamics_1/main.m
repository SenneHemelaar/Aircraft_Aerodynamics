%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%              MAIN: PANEL METHOD             %%%
%%% Authors: Jasper van Dongen & Senne Hemelaar %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear all; close all

%% Configurable parameters
n_panels = 50;
NACA = '0016';      % Choose desired 4 digit NACA airfoil NACA'INPUT'

%% Generate Geometry 
af_geo = generate_geometry(n_panels,NACA);  % af_geo stands for airfoil geometry

% Test-plot to check geometry
plot(af_geo.xU,af_geo.zU,'bo-')
hold on
plot(af_geo.xL,af_geo.zL,'bo-')
axis equal
plot(af_geo.CP_xU,af_geo.CP_zU,'go')
plot(af_geo.VP_xU,af_geo.VP_zU,'ro')
plot(af_geo.CP_xL,af_geo.CP_zL,'go')
plot(af_geo.VP_xL,af_geo.VP_zL,'ro')
plot(af_geo.xC,af_geo.zC,'k--')

%% Vortex panel method 


%% Plotting the results 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          AVL in Matlab                              %%%
%%%                     Author: Senne Hemelaar                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%%%====================SELECT AVL GEOMETRY FILE=========================%%%
f = 1;
avl_file = strcat('GEOMETRY_',string(f));

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
writeAVL(avl_file, 1.68, 1.6, 0.4, 0,...
         0, deg2rad(30), deg2rad(90), 105, 3.75,...
         28,  10,  22,   8,  22)
%--------avl_file, lw, cwr, lamw, epsr,...
%--------epst, Lamw, phiw, Sref, cref,...
%--------bref, Ncm, Nsm, Ncw, Nsw

%%%============================RUN AVL==================================%%%
runAVL(avl_file, 0.6601, 0.7) 
%------avlfile, rho, Mach

%%%===================EXTRACT FORCES FROM RESULTS=======================%%%
Forces = forces(avl_file);
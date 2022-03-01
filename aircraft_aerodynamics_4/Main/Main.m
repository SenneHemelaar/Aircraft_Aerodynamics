%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          AVL in Matlab                              %%%
%%%                     Author: Senne Hemelaar                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%%%====================SELECT AVL GEOMETRY FILE=========================%%%
f = 1;
avl_file = strcat('GEOMETRY_',string(f));

%%%==================LOADING OPTIMIZATION RESULTS=======================%%%
k_0 = load('x_vector_k_0.mat');   k_0 = k_0.x;
k_05= load('x_vector_k_0.5.mat'); k_05= k_05.x;
k_1 = load('x_vector_k_1.mat');   k_1 = k_1.x;

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
% phiw_list = linspace(0,90,91);

% for i = 1:length(phiw_list)
phiw = 90; %phiw_list(i);
writeAVL(avl_file, k_1(1), k_1(3) , k_1(4), k_1(6),...
         k_1(7), k_1(5), k_1(2), 105, 3.75,...
         28,  10,  22,   8,  22)
%--------avl_file, lw, cwr, lamw, epsr,...
%--------epst, Lamw, phiw, Sref, cref,...
%--------bref, Ncm, Nsm, Ncw, Nsw

%%%============================RUN AVL==================================%%%
runAVL(avl_file, 0.6601, 0.7) 
%------avlfile, rho, Mach

%%%===================EXTRACT FORCES FROM RESULTS=======================%%%
Forces = forces(avl_file);
% end

%%%=====================COMPUTE OBJECT FUNCTION=========================%%%
q = 16193;
k = 1;
J = optimization_function(Forces, q, k);




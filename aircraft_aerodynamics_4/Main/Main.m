%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                          AVL in Matlab                              %%%
%%%                     Author: Senne Hemelaar                          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

%%%====================SELECT AVL GEOMETRY FILE=========================%%%
f = 1;
avl_file = strcat('GEOMETRY_',string(f));

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
lw        = 0.6;
Lam       = 30;
phiw      = 0;
cwr       = 0.8;
lam       = 0.4;
epsR      = 0;
epsT      = 0;
cR        = 5.5;
cT        = 2.0;
L         = 14;
Ncm       = 10;
Ncw       = 8;
Nsm       = 22;
Nsw       = 22;

%%%====================PRELIMINARY CALCULATIONS=========================%%%
lw   = lw*2*L;         Lam  = deg2rad(Lam);
phiw = deg2rad(phiw);  cwr  = cwr * cT;
epsR = deg2rad(epsR);  epsT = deg2rad(epsT);

%%%==================LOADING OPTIMIZATION RESULTS=======================%%%
k_0 = load('x_vector_k_0.mat');   k_0  = k_0.x;
k_05= load('x_vector_k_0.5.mat'); k_05 = k_05.x;
k_1 = load('x_vector_k_1.mat');   k_1  = k_1.x;

lw   = k_1(1);
phiw = k_1(2);
cwr  = k_1(3);
lam  = k_1(4);
Lam  = k_1(5);
epsR = k_1(6);
epsT = k_1(7);

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
writeAVL(avl_file, lw, cwr, lam, epsR,...
         epsT, Lam, phiw, 105, 3.75,...
         28,  Ncm,  Nsm,   Ncw,  Nsw)
%--------avl_file, lw, cwr, lamw, epsr,...
%--------epst, Lamw, phiw, Sref, cref,...
%--------bref, Ncm, Nsm, Ncw, Nsw

%%%============================RUN AVL==================================%%%
runAVL(avl_file, 0.6601, 0.7) 
%------avlfile, rho, Mach

%%%===================EXTRACT FORCES FROM RESULTS=======================%%%
Forces = forces(avl_file);

%%%=====================COMPUTE OBJECT FUNCTION=========================%%%
q = 16193;
k = 1;
[J, M] = optimization_function(Forces, q, k);

lw   = lw/(2*L);
phiw = rad2deg(phiw);
cwr  = cwr/cT;
Lam  = rad2deg(Lam);
epsR = rad2deg(epsR);
epsT = rad2deg(epsT);




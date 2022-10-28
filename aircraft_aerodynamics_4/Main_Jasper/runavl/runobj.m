function [J,Forces,M] = runobj(x)

%%%====================SELECT AVL GEOMETRY FILE=========================%%%
f = 1;
avl_file = strcat('GEOMETRY_',string(f));
winglet = true;
%%%==========================INPUTS GUIDE===============================%%%
% lw      : Winglet Length                          [m]
% Lam     : Winglet Leading Edge Sweep              [rad]
% phiw    : Winglet Cant Angle                      [rad]
% cwr     : Winglet Root Chord                      [m]
% lam     : Winglet Taper Ratio                     [cwt/cwr]
% epsR    : Winglet Root Twist Angle                [rad]
% epsT    : Winglet Tip Twist Angle                 [rad]
% cR      : Main Wing Root Chord                    [m]
% cT      : Main Wing Tip Chord                     [m]
% b       : Main Wing Span                          [m]
% Ncm     : No. of Chordwise Sections Main Wing     [-]
% Ncw     : No. of Chordwise Sections Winglet       [-]
% Nsm     : No. of Spanwise Sections Main Wing      [-]
% Nsw     : No. of Spanwise Sections Winglet        [-]


%%%=========================DEFINE CONSTANTS============================%%%
cR  = 5.5;
cT  = 2.0;
b   = 28;
Ncm = 12;    
Ncw = 10;   
Nsm = 30;    
Nsw = 20;   
S_ref = b*cR - 0.5*b*(cR-cT); 
c_ref = (cR+cT)/2;
rho = 0.6601;
a = 316.5;
mach = 0.7;

q = 0.5*(a*mach)^2*rho;
k = 1.0;

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
writeAVL(avl_file, x(1), x(3) , x(4), x(6),...
         x(7), x(5), x(2), S_ref, c_ref,...
         b,  Ncm,  Nsm,   Ncw,  Nsw, winglet)
%--------avl_file, lw, cwr, lam, epsR,...
%--------epsT, Lam, phiw, Sref, cref,...
%--------bref, Ncm, Nsm, Ncw, Nsw


%%%============================RUN AVL==================================%%%
runAVL(avl_file, rho, mach) 
%------avlfile, rho, Mach

%%%===================EXTRACT FORCES FROM RESULTS=======================%%%
Forces = forces(avl_file, winglet);

%%%=====================COMPUTE OBJECT FUNCTION=========================%%%

[J,M] = optimization_function_J(Forces, q, k, x, winglet);
end

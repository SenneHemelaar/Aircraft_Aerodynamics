function J = runobj(x)

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

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
writeAVL(avl_file, x(1), x(3) , x(4), x(6),...
         x(7), x(5), x(2), 105, 3.75,...
         28,  10,  22,   8,  22)
%--------avl_file, lw, cwr, lam, epsR,...
%--------epsT, Lam, phiw, Sref, cref,...
%--------bref, Ncm, Nsm, Ncw, Nsw

%%%============================RUN AVL==================================%%%
runAVL(avl_file, 0.6601, 0.7) 
%------avlfile, rho, Mach

%%%===================EXTRACT FORCES FROM RESULTS=======================%%%
Forces = forces(avl_file);

%%%=====================COMPUTE OBJECT FUNCTION=========================%%%
q = 16193;
k = 0.5;
% J0 = 3.826932635843151e+06;
J = optimization_function(Forces, q, k);
% J = J/J0;
end

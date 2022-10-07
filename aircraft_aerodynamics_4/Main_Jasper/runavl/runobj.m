function [J,Forces] = runobj(x)

%%%====================SELECT AVL GEOMETRY FILE=========================%%%
f = 1;
avl_file = strcat('GEOMETRY_',string(f));

%%%====================CREATE AVL GEOMETRY FILE=========================%%%
writeAVL(avl_file, x(1), x(3) , x(4), x(6),...
         x(7), x(5), x(2), 105, 3.75,...
         28,  10,  22,   8,  22)
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
k = 0.5;
% J0 = 3.826932635843151e+06;
J = optimization_function(Forces, q, k);
% J = J/J0;
end

function [forces] = complete_AVL_run(x)
    %% set variables
    lw = x(1); cwr= x(2); lamw=x(3); epsr=x(4); 
    epst=x(5); Lamw=x(6); phiw=deg2rad(x(7));
    Sref=105; cref=3.75; bref=28;
    Ncm=40; Nsm=44; 
    Ncw=12; Nsw=22;
    
    
    avl_file = 'wing_w_winglets';
    
    %% run avl functions
    writeConfW(avl_file, lw, cwr, lamw, epsr, epst, Lamw, phiw, Sref, cref, bref, Ncm, Nsm, Ncw, Nsw);
    runAVL(avl_file, 0.6601, 0.7);
    forces = parseForces(avl_file);
    
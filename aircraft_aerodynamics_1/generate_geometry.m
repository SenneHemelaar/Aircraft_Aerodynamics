function af_geo =  generate_geometry(n_panels,NACA)

iaf.designation=NACA;           % Choose desired NACA airfoil
iaf.n=n_panels;                 % Amount of panels
iaf.HalfCosineSpacing=0;        % Cosine-spacing (yes/no)
iaf.wantFile=0;                 % Save file (yes/no)
iaf.datFilePath='./';           % Save folder folder
iaf.is_finiteTE=0;

af_geo = naca4gen(iaf);         % Fetch NACA airfoil generator function

for i = 1:iaf.n
    dist_x(i) = af_geo.xU(i) - af_geo.xU(i+1);
    dist_z(i) = af_geo.zU(i) - af_geo.zU(i+1);
    dist_rU(i) = sqrt(dist_x(i)^2+dist_z(i)^2);
    angle(i)  = atan(dist_z(i)/dist_x(i));
    
    
    % Control points upper side
    CP_xU(i)  = - (0.25*dist_rU(i)*cos(angle(i))) + af_geo.xU(i);
    CP_zU(i)  = - (0.25*dist_rU(i)*sin(angle(i))) + af_geo.zU(i);

    % Vortex elements upper side
    VP_xU(i)  = - (0.75*dist_rU(i)*cos(angle(i))) + af_geo.xU(i);
    VP_zU(i)  = - (0.75*dist_rU(i)*sin(angle(i))) + af_geo.zU(i);
    
    % Normal vectors upper side
    N_xU(i)   =   sin(angle(i));
    N_zU(i)   =   cos(angle(i));
    
    dist_x(i) =   af_geo.xL(i) - af_geo.xL(i+1);
    dist_z(i) =   af_geo.zL(i) - af_geo.zL(i+1);
    dist_rL(i) =   sqrt(dist_x(i)^2+dist_z(i)^2);
    angle(i)  =   atan(dist_z(i)/dist_x(i));

    % Control points lower side
    CP_xL(i)  =   (0.75*dist_rL(i)*cos(angle(i))) + af_geo.xL(i);
    CP_zL(i)  =   (0.75*dist_rL(i)*sin(angle(i))) + af_geo.zL(i);

    % Vortex points lower side
    VP_xL(i)  =   (0.25*dist_rL(i)*cos(angle(i))) + af_geo.xL(i);
    VP_zL(i)  =   (0.25*dist_rL(i)*sin(angle(i))) + af_geo.zL(i);
    
    % normal vectors lower side
    N_xL(i)   =   sin(angle(i));
    N_zL(i)   =   cos(angle(i));
end

% Add to struct
af_geo.CP_xL = CP_xL; af_geo.CP_zL = CP_zL;
af_geo.CP_xU = CP_xU; af_geo.CP_zU = CP_zU;
af_geo.VP_xL = VP_xL; af_geo.VP_zL = VP_zL;
af_geo.VP_xU = VP_xU; af_geo.VP_zU = VP_zU;
af_geo.N_xU  =  N_xU; af_geo.N_zU  =  N_zU;
af_geo.N_xL  =  N_xL; af_geo.N_zL  =  N_zL;
af_geo.a_U   = dist_rU; af_geo.a_L = dist_rL;


end
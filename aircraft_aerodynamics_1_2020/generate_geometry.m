function [af_geo] = generate_geometry(chord, N, camber, x_max_camber)
%GENERATE_GEOMETRY  From this function we obtain the airfoil geometry
    x      = zeros(1, N+1);
    z      = zeros(1, N+1);
    VP_x     = zeros(1, N);
    VP_z     = zeros(1, N);
    CP_x     = zeros(1, N);
    CP_z     = zeros(1, N);
    N_x     = zeros(1, N);
    N_z     = zeros(1, N);
    af_geo = struct('x',x, 'z',z, 'VP_x',VP_x, 'VP_z',VP_z, 'CP_x',CP_x, 'CP_z',CP_z, 'N_x',N_x, 'N_z',N_z);
        
    % when no camber is provided treat as flat plate
    if ~exist('camber', 'var')
        camber = false;
    end
    af_geo.x = linspace(0, 1, N+1);
    dp = chord/N;
    if camber
        if ~exist('x_max_camber', 'var') || x_max_camber == 0
            x_max_camber = 0.2;
        end
        % NACA four-digit formula (From Matlab library)
        for i=1:N+1
            if af_geo.x(i) <= x_max_camber
                af_geo.z(i)= camber/(x_max_camber^2) * ...
                    (2*x_max_camber * af_geo.x(i)-af_geo.x(i)^2);
            	if i >= 2
                    af_geo.VP_z(i-1) = af_geo.z(i)-0.75*norm(af_geo.z(i)-af_geo.z(i-1));
                    af_geo.CP_z(i-1) = af_geo.z(i)-0.25*norm(af_geo.z(i)-af_geo.z(i-1));
                end   
            else
                af_geo.z(i)= camber*(1-x_max_camber)^-2 * ...
                    ((1- 2 * x_max_camber)+ 2*x_max_camber * af_geo.x(i)-af_geo.x(i)^2);
                if i >= 2
                    af_geo.VP_z(i-1) = af_geo.z(i)+0.75*norm(af_geo.z(i)-af_geo.z(i-1));
                    af_geo.CP_z(i-1) = af_geo.z(i)+0.25*norm(af_geo.z(i)-af_geo.z(i-1));
                end 
            end 
        end 
    end
    
    % Dimensionalize
    af_geo.x = af_geo.x * chord;
    af_geo.z = af_geo.z * chord;
    af_geo.VP_x = af_geo.x(2:end)-0.75*dp;
    af_geo.CP_x = af_geo.x(2:end)-0.25*dp;
    af_geo.N_x = -(af_geo.CP_z-af_geo.VP_z)./sqrt((af_geo.CP_x-af_geo.VP_x).^2+(af_geo.CP_z-af_geo.VP_z).^2);
    af_geo.N_z = (af_geo.CP_x-af_geo.VP_x)./sqrt((af_geo.CP_x-af_geo.VP_x).^2+(af_geo.CP_z-af_geo.VP_z).^2);
end


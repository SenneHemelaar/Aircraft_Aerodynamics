function v_g = build_vorticity_geometry(chord, n_panels, camber, location_max_camber)
%BUILD_VORTICITY_GEOMETRY the function to build a 2D geometry of the
% aerofoil to use in the vortex panel method.
    x_array = zeros(1, n_panels+1);
    z_array = zeros(1, n_panels+1);
    xv_array = zeros(1, n_panels);
    zv_array = zeros(1, n_panels);
    xc_array = zeros(1, n_panels);
    zc_array = zeros(1, n_panels);
    xn_array = zeros(1, n_panels);
    zn_array = zeros(1, n_panels);
    v_g = struct('x', x_array, 'z',z_array,'x_v', xv_array, 'z_v', zv_array, 'x_c', xc_array, 'z_c', zc_array,'x_n', xn_array, 'z_n', zn_array);
        
    % when no camber is provided treat as flat plate
    if ~exist('camber', 'var')
        camber = false;
    end
    
    v_g.x = linspace(0, 1, n_panels+1);
    dp = chord/n_panels;
    % logic to set up the z direction based on four-digit NACA series
    if camber
        
        %Make sure the location of maximum camber exists
        if ~exist('location_max_camber', 'var') || location_max_camber == 0
            location_max_camber = 0.2;
        end
        
        %% NACA four-digit formula
        for i=1:n_panels+1
            if v_g.x(i) <= location_max_camber
                v_g.z(i)= camber/(location_max_camber^2) * ...
                    (2*location_max_camber * v_g.x(i)-v_g.x(i)^2);
                
            	if i >= 2
                    v_g.z_v(i-1) = v_g.z(i)-0.75*norm(v_g.z(i)-v_g.z(i-1));
                    v_g.z_c(i-1) = v_g.z(i)-0.25*norm(v_g.z(i)-v_g.z(i-1));
                end   
                
            else
                v_g.z(i)= camber*(1-location_max_camber)^-2 * ...
                    ((1- 2 * location_max_camber)+ 2*location_max_camber * v_g.x(i)-v_g.x(i)^2);
                
                if i >= 2
                    v_g.z_v(i-1) = v_g.z(i)+0.75*norm(v_g.z(i)-v_g.z(i-1));
                    v_g.z_c(i-1) = v_g.z(i)+0.25*norm(v_g.z(i)-v_g.z(i-1));
                end 
            end 

        end 
    end
    
    %% dimensionalize output
    
    v_g.x = v_g.x * chord;
    v_g.z = v_g.z * chord;
    
    
    v_g.x_v = v_g.x(2:end)-0.75*dp;
    v_g.x_c = v_g.x(2:end)-0.25*dp;

    
    v_g.x_n = -(v_g.z_c-v_g.z_v)./sqrt((v_g.x_c-v_g.x_v).^2+(v_g.z_c-v_g.z_v).^2);
    v_g.z_n = (v_g.x_c-v_g.x_v)./sqrt((v_g.x_c-v_g.x_v).^2+(v_g.z_c-v_g.z_v).^2);
   

    %% debug to check geometry of airfoil
%   plot(vortex_geometry.x, vortex_geometry.z);
%   axis equal;
    
end


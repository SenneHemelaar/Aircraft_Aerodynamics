%% Find Loads for each blade element
function [fnorm , ftan, alpha, inflowangle, gamma] =... 
    solveAzimuthalElement(U_0, R, r_R, chord, twist, polar_alpha, polar_cl, polar_cd, a, a_prime, omega, yaw)
    n_azimuthal_sections = 50;
    fnorm_array = zeros(1, n_azimuthal_sections);
    ftan_array = zeros(1, n_azimuthal_sections);
    alpha_array = zeros(1, n_azimuthal_sections);
    inflowangle_array = zeros(1, n_azimuthal_sections);
    gamma_array = zeros(1, n_azimuthal_sections);


     K = coleman_model(find_skew_angle(a, yaw));
     u_i_0 = a * U_0;
    
    for i=1:n_azimuthal_sections
        azimuthal_angle = 2 * pi/n_azimuthal_sections * i; % average azimuthal angle for section
        u = u_i_0 * (1 + K * r_R * sin(azimuthal_angle));
        U_rel = azimuthal_relative_velocity(U_0, yaw, u); %average velocity
        vnorm = U_rel * (1-a); %velocity normal to 
        vtan = U_0* sin(yaw) * cos(azimuthal_angle) + (1+a_prime)*omega*r_R*R;
        vmag2 = vnorm^2 + vtan^2;
        
        inflowangle_array(i) = atan2(vnorm,vtan); %phi
        alpha_array(i) = twist/180*pi + inflowangle_array(i); %alpha

        cl = interp1(polar_alpha, polar_cl, alpha_array(i)); %extract lift coefficient from test data
        cd = interp1(polar_alpha, polar_cd, alpha_array(i)); %extract lift coefficient from test data
        lift = 0.5*vmag2*cl*chord; %compute lift force
        drag = 0.5*vmag2*cd*chord; %compute drag force
        fnorm_array(i) = lift*cos(inflowangle_array(i))+drag*sin(inflowangle_array(i)); %force in flow normal direction
        ftan_array(i) = lift*sin(inflowangle_array(i))-drag*cos(inflowangle_array(i)); %force in flow tangential direction
        gamma_array(i) = 0.5*sqrt(vmag2)*cl*chord;  %circulation
    end
    
    inflowangle = sum(inflowangle_array)/n_azimuthal_sections* 180 / pi;
    alpha = sum(alpha_array)/n_azimuthal_sections * 180 / pi;
    fnorm = sum(fnorm_array)/n_azimuthal_sections;
    gamma = sum(gamma_array)/n_azimuthal_sections;
    ftan = sum(ftan_array)/n_azimuthal_sections;
end
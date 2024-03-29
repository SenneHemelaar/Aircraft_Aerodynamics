%% Find Loads for each blade element
function [fnorm , ftan, gamma, alpha, inflowangle] = loadBladeElement(vnorm, vtan, r_R, chord, twist, polar_alpha, polar_cl, polar_cd)
    vmag2 = vnorm^2 + vtan^2;
    inflowangle = atan2(vnorm,vtan); %phi
    alpha = twist/180*pi + inflowangle; %alpha
    cl = interp1(polar_alpha, polar_cl, alpha);
    cd = interp1(polar_alpha, polar_cd, alpha);
    lift = 0.5*vmag2*cl*chord;
    drag = 0.5*vmag2*cd*chord;
    
    %returned variables
    fnorm = lift*cos(inflowangle)+drag*sin(inflowangle); %force in flow normal direction
    ftan = lift*sin(inflowangle)-drag*cos(inflowangle); %force in flow tangential direction
    gamma = 0.5*sqrt(vmag2)*cl*chord;  %circulation
    alpha = alpha * 180 / pi; %convert back to degrees 
    inflowangle = inflowangle* 180 / pi; %convert back to degrees 
end


function [J] = complete_obj_func(x)
    %% retrieve global variables
    global q
    global k
    global x_progress
    
    %% track development of x
    b = abs(x_progress); % it is a trick which forces min to do not consider your negative numbers as minimum 
    [~,ind] = min(b);
    
    x_progress(ind(1), :) = x;    

    %% compute objective function TODO INCLUDE STRIP2
    forces = complete_AVL_run(x);
    shear = [forces.strip1(:,7);forces.strip2(:,7)] .* [forces.strip1(:,4); forces.strip2(:,4)] * q;
    moment = sum(shear .* [forces.strip1(:,2);forces.strip2(:,2)]);
    J  = k * forces.CDind/ 0.0151 + (1-k) *  moment/ 1.2927e+06
end

% 
% moment/ 3.0499e+06
% 
% ans =
% 
%     1.0474
% 
% forces.CDind/0.0162319
% 
% ans =
% 
%     0.8604

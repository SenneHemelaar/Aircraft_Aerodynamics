%% The vortex panel method

gamma = zeros(n_panels, n_alpha);
Dp = zeros(n_panels, n_alpha);
CL =  zeros(1, n_alpha);


%% Loop over angles of attack

for i=1:n_alpha
    alpha = alphas(i);
    
    gamma(:, i) =  calculating_gamma(v_g, n_panels,u_0(i),w_0(i),v_g.x_n,v_g.z_n);
    
    
    % Calculating velocity and pressure fields

    [abs_vel_field, u_vel_field, w_vel_field, pres_field]  = calculating_fields(gamma(:,15), v_g.x_v, v_g.z_v, x_field, z_field, rho, u_0(15), w_0(15));

    

    % Calculating and storing load values over airfoil 
    [CL_i, Dp_i] = calculating_loads(gamma(:,i),rho,U_0,chord,n_panels);
    
    CL(:,i) = CL_i;
    Dp(:,i) = Dp_i;
end


%% The unsteady vortex panel method



%% Getting the geometry of the aerofoil
v_g = build_vorticity_geometry(chord, n_panels, camber, location_max_camber);


gamma_b = [];     % epmty list for bound vorticep at each panel
Dp = zeros(n_panels, max_timestep);          % empty list for pressure difference
CL =  zeros(1, max_timestep);                % empty list for lift coefficients

error_list = [];

%% begin time loop

time_period = ceil(2 * pi /omega/dt) + 1; % time steps required for one loop

% empty list for strengths of shedded vortices 
%[FIRST IN LIST IS FIRST SHEDDED VORTEX] 
gammas_shedded = [];
x_v_w = zeros(1, max_timestep) + v_g.x(end);
z_v_w = zeros(1, max_timestep) + v_g.z(end);
previous_period_lift = zeros(ceil(time_period),1);



abs_vel_field = zeros(length(x_field), length(z_field),time_period);
pres_field = zeros(length(x_field), length(z_field),time_period); 


for n = 1:max_timestep
    period_lift = zeros(time_period,1);
    for t=1:time_period
        N = (n-1)*time_period + t;
        
        % finding coordinates of shedded vortices and adding them to their list
        %[FIRST IN LIST IS FURTHEST SHEDDED VORTEX] 

        alpha_list = [alpha_list, alpha + deg2rad(5*sin(omega*N*dt))];
        u_0 = U_0 * cos(alpha_list(N));
        w_0 = U_0 * sin(alpha_list(N));

        x_v_w(1:N) = x_v_w(1:N) + u_0*dt;
        z_v_w(1:N) = z_v_w(1:N) + w_0*dt;

        %[FIRST IN LIST IS MOST RECENT SHEDDED VORTEX] 
        x_w = flip(x_v_w(1:N));
        z_w = flip(z_v_w(1:N));

        [gamma_b(:,N), gamma_w] =  calculating_unsteady_gammas(v_g, n_panels, u_0, w_0, N, gammas_shedded, x_w, z_w);
        gammas_shedded = [gammas_shedded, gamma_w];

        % Calculating and storing load values over airfoil 
    %     [CL_i, Dp_i] = calculating_loads(gamma_b(:,t_n),rho,U_0,chord,n_panels);

        % Calculating and storing load values over airfoil 
        if N > 1
            [CL_i, Dp_i] = calculating_unsteady_loads(gamma_b(:,N-1:N), v_g, rho, U_0, chord, n_panels, alpha_list(N), dt);
%             [CL_i, Dp_i] = calculating_unsteady_loads_cyril(gamma_b, N, dt, v_g, rho, U_0, chord, n_panels);
            CL(:,N) = CL_i;
            Dp(:,N) = Dp_i;
            period_lift(t) = CL_i;
        end
    end

    %track error
    error = abs(mean(period_lift) - mean(previous_period_lift))/mean(period_lift);
    error_list = [error_list, error];
    
    %check for convergence
    if error<max_error
        disp('converged');
        
        % Calculating velocity and pressure fields
        for t_n=1:time_period
            t = t_n-1;
            [abs_vel_field(:,:,end-t), u_vel_field, w_vel_field ,pres_field(:,:,end-t)]  = calculating_unsteady_fields(gamma_b(:,end-t), gammas_shedded(1:end-t), v_g.x_v, v_g.z_v, x_field, z_field, rho, u_0, w_0, x_v_w(1:end-t), z_v_w(1:end-t));
        end
        break;
    end
    disp(error)
    previous_period_lift = period_lift;

end % end loop over time
% Main script for assignment 3 meant to contain settings, transfer to the next script and plotting relevant figures. 
clear;clc;close;

%% Settings
% Select (un)steady panel code: 
% true for steady, false for unsteady
steady = true;

% Geometry Settings & Build
n_panels = 10;                                          % number of panels
chord = 1;                                              % mean aerodynamic chord (MAC)                                      
camber = 0.02;                                          % max camber, to model zero camber set to zero
location_max_camber = 0.2;                              % location max camber, no need to set to zero when no camber is applied. 
v_g = build_vorticity_geometry(chord, n_panels, camber, location_max_camber);

% General Flow Parameters
U_0 = 10;                                               % inflow velocity [m/s ]
rho = 1.225;                                            % air density

% Steady Flow

if steady   
    % Steady Flow paramaters
    alpha_max = deg2rad(20);                            % maximum angle of attack [rad]
    alpha_min = deg2rad(-10);                           % minimum angle of attack [rad]
    n_alpha = 30;                                       % number of angles of attack to evaluate
    alphas = linspace(alpha_min, alpha_max, n_alpha);  	% list of angles of attack
    alphas_deg = rad2deg(alphas);                     	% list of aoa's in deg. NOTHING IS DONE WITH THESE VALUES, ONLY DISPLAYING FOR OVERVIEW CALCULATED ANGLES
    u_0 = U_0*cos(alphas);                              
    w_0 = U_0*sin(alphas);
    
    % Field boundaries and and grid density (for plotting)
    x_field = linspace(-1,2, 12);
    z_field = linspace(-1.5, 1.5, 12);
    
    % Running Steady Panel Method
    vortex_panel_method;

    % Plotting Resutls for Steady Panel Method
    plotting_steady;

else
    
% Unsteady Flow

    t_array = [0.1, 0.05, 0.03];
    for dt_step=1:length(t_array)
        % Unsteady Flow Parameters
        dt = 0.01;                                           % timestep
        dt = t_array(dt_step);                           % reduced frequency
        k = 0.1;                                          % reduced frequency
        omega = k*2*U_0/chord;                              % anglular velocity (period)
        alpha = 0;                                      	% angle of attack
        u_0 = U_0*cos(alpha);                               % freestream velocity component in x-direction
        w_0 = U_0*sin(alpha) ;                              % freestream velocity component in z-direction
        max_timestep = 100000;                              	% number of timesteps
        alpha_list = [];                                    % angle of attack at each timestep
        max_error = 1e-3;                                    % max error

        % Field boundaries and and grid density (for plotting)
        x_field = linspace(-1.5,2, 30);
        z_field = linspace(-2, 2, 30);

        % Running Unsteady Panel Method
        vortex_panel_method;

        % Plotting Results for Unsteady Panel Method
        plotting_steady;
    end
end

figure(5)
plot(v_g.x,v_g.z)
hold on
plot(v_g.x_n,v_g.z_n,'o')



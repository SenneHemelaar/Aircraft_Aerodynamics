%% BEM Group assignment

% Group 23's approach at the BEM script as explained in the lectures
clear all;
close all;

%% Settings
variable_TSR = false;
variable_yaw = false;
iteratingDesignParameters = false;
cosine_radial_distribution = false;

%% Assignment variables
pitch = 2; %degrees
da = 10^-9; % step size for da
da_prime = 10^-9; % step size for da_prime
dtheta = 180 / pi; % step for iteration of the twist angle
yaw_angles = [0, 15, 30]; %angles of yaw, non-zero should be tested with TSR = 8
TSR_array = [6, 8, 10]; %Tip Speed Ratio
% annulus_array = [20, 50, 150];
R = 50; %m radius of blade
U_0 = 10;  %m/s, upstream wind velocity
n_blades = 6; %number of blades 
rho = 1.225; %kg/m^3, density of air
P_0 = 101325;


%for i=50:100:250

%% Setting up the grid
annuli = 50; % number of annuli across cross section
if (cosine_radial_distribution) 
    r_R = cosspace(0.2, 1, annuli);
else
    r_R = linspace(0.2, 1, annuli);
end
r = r_R * R;
rootradius_R = 0.1; %measured in r/R
tipradius_R = 1; %measured in r/R

%% Iteration parameters
error_limit = 10^-5; %error limit
alfa = 7; % Angle of attack opti
Cp_limit = 0.45; %Power factor limit for convergence 
max_iterations = 600; %maximum allowed number of iterations

chord_dist = 3 * (1-r_R) + 1; % chord of blade in [m]
twist_dist = -14 * (1-r_R); %twist over the radius of the blade
mean_chord = mean(chord_dist); %average chord of blade

C_p_array = [];

%% Import values

% Importing values from airfoil
airfoilDat = readtable('polar_DU95W180.xlsx');

alfaMat = airfoilDat.Alfa; %vector of angle of attack values 
alfaMat = alfaMat * pi/180; %to radians
clMat = airfoilDat.Cl; %vector of lift coefficient values
cdMat = airfoilDat.Cd; %vector of drag coefficient values. 
cmMat = airfoilDat.Cm; %vector of moment coefficient values. 

%% Execution based on input values
%Array of results in the following order: [a , a_prime, r_R_annulus, f_normal , f_tangential, gamma]

if variable_TSR
    TSR_loop_start = 1;
    TSR_loop_size = length(TSR_array);
else
    TSR_loop_start = 2;
    TSR_loop_size = 2;
    TSR = 8;
end

if variable_yaw
    yaw_loop_start = 1;
    yaw_loop_end = length(yaw_angles);
else
    yaw_loop_start = 1;
    yaw_loop_end = 1;
end

plot_legends = struct('legend_induction', strings(), 'legend_forces', strings(), 'legend_circulation', strings(), 'legend_angles', strings());
azimuthal_angles = linspace(0,2*pi,50);
% for n = 1:length(annulus_array)
%     annuli = annulus_array(n);
    r_R = linspace(0.2, 1, annuli);
    r = r_R * R;
    chord_dist = 3 * (1-r_R) + 1; % chord of blade in [m]
    twist_dist = -14 * (1-r_R); %twist over the radius of the blade

    for y = yaw_loop_start:yaw_loop_end
        yaw = yaw_angles(y) /180 * pi;
        disp('yaw');
        disp(yaw);
        for n = TSR_loop_start:TSR_loop_size %length(TSR_array)
            results = zeros(length(r_R)-1, 10);
            matrix_alpha = zeros(length(r_R)-1,50);
            matrix_inflowangle = zeros(length(r_R)-1,50);
            matrix_fnorm = zeros(length(r_R)-1,50);
            TSR = TSR_array(n);

            omega = TSR*U_0/R;
            area = pi*(R^2-(R*0.2)^2);
            dr = zeros(length(r_R)-1, 1);
            dArea = zeros(length(r_R)-1, 1);

            for i = 1:length(r_R)-1
                chord = interp1(r_R, chord_dist, (r_R(i+1)+r_R(i))/2);
                twist = interp1(r_R, twist_dist, (r_R(i+1)+r_R(i))/2);
                j = i -1;   %index of previous results
                if j  == 0  %cannot go before start of array
                    j = 1;
                end

                [results(i,1), results(i,2), results(i,3), results(i,4), results(i,5), results(i,6), results(i,7), results(i,8)] =... 
                    solveStreamtube(U_0, rho, r_R(i), r_R(i+1), rootradius_R, tipradius_R, omega, R, n_blades, chord, twist, alfaMat, clMat, cdMat, yaw, max_iterations, error_limit, iteratingDesignParameters);
                dr(i) = (r(i+1) - r(i));
                dArea(i)= dr(i)^2 * pi;
            end

            power_coefficient = dr.*results(:,5).*results(:,3)*R*n_blades*omega/(0.5*U_0^3*rho* area); %.* dArea.^-1;
    %         CP = sum(dr.*results(:,5).*results(:,3)*n_blades*R*omega/(0.5*U_0^3*rho) .* dArea.^-1);
            CP = sum(power_coefficient);
            disp('CP');
            disp(CP);
            force_dimension =(.5 * U_0^2 *  rho* R); 
            total_torque = sum(dr.*results(:,5).*results(:,2)*n_blades*R)        %total torque
            total_thrust = sum(results(:,5)*n_blades)            %total thrust
            results(:,4) = results(:,4)./force_dimension;
            results(:,5) = results(:,5)./force_dimension;
            C_Q = dr.*results(:,5).*results(:,2)* n_blades/area;
            circulation_dimension = (pi*U_0^2/(n_blades*omega));
            power_dimension = n_blades*R*omega/(0.5*U_0^3*pi*R^2);


            %% Enthalpy calculations
%             enthalpy_plots(R, results,P_0, rho, r_R, U_0)

            %% twist and chord distribution
%             twist_chord_plots(results)

            

            %% plot legends 
            plot_legends = plotGraphs(results, C_Q, power_coefficient, TSR, yaw, circulation_dimension, plot_legends, variable_TSR, variable_yaw);
        end
    %     figure(1)
    %     plot(plot_legends.legend_induction)
    end
% end


%% polarplot

matrix_alpha = (180/pi)* matrix_alpha;
[azimuthal_angle,r_R] = meshgrid(linspace(0,2*pi,50),results(:,3));
[x,y_1] = pol2cart(azimuthal_angle,r_R);

% figure('Color','white','NumberTitle','off','Name','Angle of attack','Position', [100 100 900 600]);
% polarplot3d(matrix_alpha,'PlotType','surfn','PolarGrid',{4 12},'TickSpacing',15,...
% 'AngularRange',[0 360]*pi/180,'RadialRange',[0.2 1],...
% 'RadLabels',3,'RadLabelLocation',{15 'max'},'RadLabelColor','black','GridColor','#94b8b8','PolarAxisColor','#94b8b8','TickColor','white');
% set(gca,'DataAspectRatio',[1 1 10],'View',[0,90],...
% 'Xlim',[-1.2 1.2],'Xtick',[],...
% 'Ylim',[-1.2 1.2],'Ytick',[],'ycolor','white','xcolor','white');
% title('Angle of attack \alpha [deg]'); colorbar('eastoutside'); caxis([0 12]); colormap(turbo(24))
% 
% matrix_inflowangle = (180/pi)* matrix_inflowangle;
% figure('Color','white','NumberTitle','off','Name','Inflow Angle','Position', [100 100 900 600]);
% polarplot3d(matrix_inflowangle,'PlotType','surfn','PolarGrid',{4 12},'TickSpacing',15,...
% 'AngularRange',[0 360]*pi/180,'RadialRange',[0.2 1],...
% 'RadLabels',3,'RadLabelLocation',{15 'max'},'RadLabelColor','black','GridColor','#94b8b8','PolarAxisColor','#94b8b8','TickColor','white');
% set(gca,'DataAspectRatio',[1 1 10],'View',[0,90],...
% 'Xlim',[-1.2 1.2],'Xtick',[],...
% 'Ylim',[-1.2 1.2],'Ytick',[],'ycolor','white','xcolor','white');
% title('Inflow angle \phi [deg]'); colorbar('eastoutside'); 
% caxis([0 24]); 
% colormap(turbo(24))

matrix_CT = matrix_fnorm .* dr ./ force_dimension;


%determine thrust coefficient for each azimuthal section
% figure('Color','white','NumberTitle','off','Name','Inflow Angle','Position', [100 100 900 600]);
% polarplot3d(matrix_CT,'PlotType','surfn','PolarGrid',{4 12},'TickSpacing',15,...
% 'AngularRange',[0 360]*pi/180,'RadialRange',[0.2 1],...
% 'RadLabels',3,'RadLabelLocation',{15 'max'},'RadLabelColor','black','GridColor','#94b8b8','PolarAxisColor','#94b8b8','TickColor','white');
% set(gca,'DataAspectRatio',[1 1 10],'View',[0,90],...
% 'Xlim',[-1.2 1.2],'Xtick',[],...
% 'Ylim',[-1.2 1.2],'Ytick',[],'ycolor','white','xcolor','white');
% title('C_T'); colorbar('eastoutside'); 
% caxis([0 0.4]); 
% colormap(parula(24))

%determine azial induction
% a_1 = 0.5*(1-sqrt(1-CT));

% figure('Color','white','NumberTitle','off','Name','Inflow Angle','Position', [100 100 900 600]);
% polarplot3d(a_1,'PlotType','surfn','PolarGrid',{4 12},'TickSpacing',15,...
% 'AngularRange',[0 360]*pi/180,'RadialRange',[0.2 1],...
% 'RadLabels',3,'RadLabelLocation',{15 'max'},'RadLabelColor','black','GridColor','#94b8b8','PolarAxisColor','#94b8b8','TickColor','white');
% set(gca,'DataAspectRatio',[1 1 10],'View',[0,90],...
% 'Xlim',[-1.2 1.2],'Xtick',[],...
% 'Ylim',[-1.2 1.2],'Ytick',[],'ycolor','white','xcolor','white');
% title('Axial induction 1'); colorbar('eastoutside'); 
% caxis([0 24]); 
% colormap(parula(24))


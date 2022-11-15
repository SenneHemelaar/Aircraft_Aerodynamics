clear all 
close all
clc

N = 4;
D = 3.04;
tipR = D/2;
rootR = 0.1*tipR;
RPM = 400;
n = RPM/60;
omega = n*pi*2;
rho = 1.225;
sec = [0.2 0.3 0.45 0.6 0.7 0.8 0.9 0.95];

r_steps = length(sec);
dr = [0.1 diff(sec)]*D;

%Import lift data

fid = csvread('CL_alpha.csv',2);

Cl_a.sec_1 = [-4 6; 0.0001 0.0001];
Cl_a.sec_2 = fid(1:2,1:2)';
Cl_a.sec_3 = fid(1:2,3:4)';
Cl_a.sec_4 = fid(1:2,5:6)';
Cl_a.sec_5 = fid(1:2,7:8)';
Cl_a.sec_6 = fid(1:2,9:10)';
Cl_a.sec_7 = fid(1:2,11:12)';
Cl_a.sec_8 = fid(1:2,13:14)';
Cl_a_sec_fn = fieldnames(Cl_a);
%Import L/D data

fid = csvread('L_D_alpha.csv',2);

Cd_a.sec_1 = [-6 6; 0.00025 0.00025];
Cd_a.sec_2 = [-6 6; 0.001 0.001];

a = nonzeros(fid(:,1)');
b = nonzeros(fid(:,2)');
Cd_a.sec_3 = [a';b'];

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
Cd_a.sec_4 = [a';b'];

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
Cd_a.sec_5 = [a';b'];

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
Cd_a.sec_6 = [a';b'];

a = nonzeros(fid(:,9)');
b = nonzeros(fid(:,10)');
Cd_a.sec_7 = [a';b'];

a = nonzeros(fid(:,11)');
b = nonzeros(fid(:,12)');
Cd_a.sec_8 = [a';b'];

Cd_a_sec_fn = fieldnames(Cd_a);

%Import chord length distribution

fid = csvread('Chord_Pitch_radial.csv',2);

chord = fid(:,1:2)';
chord(2,:) = chord(2,:)*D;

%Import pitch distribution

fid = csvread('Chord_Pitch_radial.csv',2);

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
pitch.deg_25 = [a';b'];
pitch.deg_25(2,:) = pitch.deg_25(2,:)*D;

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
pitch.deg_35 = [a';b'];
pitch.deg_35(2,:) = pitch.deg_35(2,:)*D;
% beta = 35 - interp1(pitch.deg_35(1,1:6),pitch.deg_35(2,1:6),.75,'spline');
% pitch.deg_35(2,:) = pitch.deg_35(2,:) + beta;

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
pitch.deg_45 = [a';b'];
pitch.deg_45(2,:) = pitch.deg_45(2,:)*D;
% beta = 45 - interp1(pitch.deg_45(1,1:6),pitch.deg_45(2,1:6),.75,'spline');
% pitch.deg_45(2,:) = pitch.deg_45(2,:) + beta;






fail = 0;
for V_inf=1:50
    fprintf('\n     V_inf = %.2f     \n',V_inf)
    TSR = n*tipR/V_inf;
    %initialise sums
    thrust=0.0;
    torque=0.0;
    %loop over each blade element
    for j=1:r_steps
        
        r = sec(j)*(D/2); % add later -0.5*dr(j);
        r_R = r/tipR;
        theta = interp1(pitch.deg_25(1,:),pitch.deg_25(2,:),sec(j),'spline');
        theta = atan(theta/2/pi/r);
        a_i_0=0.1;
        b_i_0=0.01;
        %set logical variable to control iteration
        finished=false;
        %set iteration count and check flag
        sum=1;
        itercheck=0;
        
        while (~finished)
            %%%====Equations 3 & 4====%%%
            %%% Calculate velocities
            %axial velocity
            Vax=V_inf*(1+a_i_0);
            %disk plane velocity
            Vtan=omega*r*(1-b_i_0);
            
            % Equation 3
            %local velocity at blade
            Vp=sqrt(Vax^2+Vtan^2);
            
            % Equation 4
            %flow angle %which flow angle
            phi=atan2(Vax,Vtan);
            %blade angle of attack 
            alpha=theta-phi;
            
            %%%====Equations 1 & 2 ====%%%
            %%% Get lift and drag coefficients
            % lift coefficient
            if j == 1
                cl = 0;
            elseif j == 2
                cl = 0.045*alpha;
            else
                cl = interp1(Cl_a.(Cl_a_sec_fn{j})(1,:), ... 
                             Cl_a.(Cl_a_sec_fn{j})(2,:), ...
                             rad2deg(alpha));
            end
            
            % drag coefficient
            if j == 1
                cd = 0.4;
            elseif j == 2
                cd = 0.1;
            else 
                cd = 1/(interp1(Cd_a.(Cd_a_sec_fn{j})(1,:), ... 
                                Cd_a.(Cd_a_sec_fn{j})(2,:), ...
                                rad2deg(alpha),'spline')/cl);
            end
            % chord 
            
            chord_j = interp1(chord(1,:), ... 
                              chord(2,:), ...
                              sec(j),'spline');
                          
            % Equation 1 divided by dr
            DTdr = 0.5 * rho * chord_j * Vp^2 * (cl*cos(phi) - cd*sin(phi)) * N;
            
            % Equation 2 divided by dr
            DQdr = 0.5 * rho * chord_j * Vp^2 * (cd*cos(phi) + cl*sin(phi)) * N*r;

            %%%====Equations 5 & 6====%%%
            %inflow and swirl
            a_i_1 = DTdr/(rho * 4 * pi * r * V_inf^2 * (1+a_i_0));
            b_i_1 = DQdr/(rho * 4 * pi * r^3 * V_inf * (1+a_i_0) * omega);
            
            %%% Prandtl tip and root losses %%% 
            
            
            temp1 = -N/2 * ((tipR/tipR)-r_R)/r_R *sqrt( 1+ ((TSR*r_R)^2)/((1-a_i_1)^2));
            F_tip = 2/pi * acos(exp(temp1));
            F_tip(isnan(F_tip)) = 0;
            
            temp2 = N/2 * ((rootR/tipR)-r_R)/r_R *sqrt( 1+ ((TSR*r_R)^2)/((1-a_i_1)^2));
            F_root = 2/pi * acos(exp(temp2));            
            F_root(isnan(F_root)) = 0;
            
            pr_corr = F_tip*F_root;
            
            if pr_corr < 1e-3
                pr_corr = 1e-2;
            end
            %stabilise iteration
            a_i_1 = 0.5 * a_i_0 + 0.5 * a_i_1*pr_corr; 
            b_i_1 = 0.5 * b_i_0 + 0.5 * b_i_1*pr_corr;
            
            %check for convergence
            if (abs(a_i_1-a_i_0)<1.0e-5) && (abs(b_i_1-b_i_0)<1.0e-5)
                finished=true;
            end


            a_i_0=a_i_1;
            b_i_0=b_i_1;
            %increment iteration count
            sum=sum+1;
            %check to see if iteration stuck
            if (sum>20)
                finished=true;
                itercheck=1;
                fail = fail + 1;
            end
            
        end
        thrust=thrust+DTdr*dr(j);
        torque=torque+DQdr*dr(j);

        fprintf('j = %i |  cl = %.2f, cd = %.2f, r =  %.2f, dr =  %.2f,',[j cl cd r dr(j)]);
        fprintf('chord = %.2f, pitch = %.2f, alpha = %.2f, Ftip = %.2f, Froot = %.2f \n',[chord_j,rad2deg(theta),rad2deg(alpha),F_tip,F_root]);
    
    end
    CT(V_inf)=thrust/(rho*n^2*D^4);
    CP(V_inf)=torque/(rho*n^2*D^5);
    J(V_inf)=V_inf/(n*D);
    alpha_lib(V_inf) = alpha;
    if CT(V_inf) < 0
        eff(V_inf) = 0;
    else
        eff(V_inf)=J(V_inf)/2.0/pi*(CT(V_inf)/CP(V_inf));
    end
    
end

%% 
%Import Coefficient of Thrust data
fid = csvread('Ct_J.csv',2);

a = nonzeros(fid(:,1)');
b = nonzeros(fid(:,2)');
Ct_j.deg_25 = [a';b'];

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
Ct_j.deg_35 = [a';b'];

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
Ct_j.deg_45 = [a';b'];

%Import Coefficient of Power data
fid = csvread('C_p_eta_J.csv',2);

a = nonzeros(fid(:,1)');
b = nonzeros(fid(:,2)');
Cp_j.deg_25 = [a';b'];

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
Cp_j.deg_35 = [a';b'];

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
Cp_j.deg_45 = [a';b'];

%Import eta data 
fid = csvread('C_p_eta_J.csv',2);

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
eta_j.deg_25 = [a';b'];

a = nonzeros(fid(:,9)');
b = nonzeros(fid(:,10)');
eta_j.deg_35 = [a';b'];

a = nonzeros(fid(:,11)');
b = nonzeros(fid(:,12)');
eta_j.deg_45 = [a';b'];


% plotting the results
%Jmax=max(J);
%CTmax=max(CT);


%fit curves
f_Cp = fit(Cp_j.deg_25(1,:)',Cp_j.deg_25(2,:)','smoothingspline');
f_Ct = fit(Ct_j.deg_25(1,:)',Ct_j.deg_25(2,:)','smoothingspline');
f_eta = fit(eta_j.deg_25(1,:)',eta_j.deg_25(2,:)','smoothingspline');

f1 = figure(1);
hold on 
plot(J,CT,J,CP)
plot(f_Ct,'m',Ct_j.deg_25(1,:)',Ct_j.deg_25(2,:)')
plot(f_Cp,'c',Cp_j.deg_25(1,:)',Cp_j.deg_25(2,:)')

title('Thrust and Torque Coefficients')
xlabel('Advance Ratio (J)');
ylabel('C_T, C_P');
legend('C_T computed','C_P computed', ...
       'C_T experimental', 'C_T experimental fitted', ...
       'C_P experimental', 'C_P experimental fitted','Location','best');
   
hold off


f2 = figure(2);
hold on
plot(J,eff);
plot(f_eta,'m',eta_j.deg_25(1,:)',eta_j.deg_25(2,:)')
title('Propeller Efficiency');
xlabel('Advance Ratio (J)');
ylabel('\eta');
legend('\eta computed', '\eta experimental', '\eta experimental fitted');

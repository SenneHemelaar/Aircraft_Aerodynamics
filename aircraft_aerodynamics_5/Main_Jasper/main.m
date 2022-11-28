clear all 
close all
clc

N = 4;
D = 3.04;
tipR = D/2;
rootR = 0.2*tipR;
RPM = 400;
n = RPM/60;
omega = n*pi*2;
rho = 1.225;
sec = [0.2 0.3 0.45 0.6 0.7 0.8 0.9 0.95];

r_steps = 9;
sections = linspace(rootR/tipR,tipR/tipR,r_steps);

dr = [diff(sections) 0.05];

%Set if pradtl tip and root correction is used
PrandtlCorrection_deg_25 = true;
PrandtlCorrection_deg_35 = true;
PrandtlCorrection_deg_45 = true;

DiagnosticInfo = true;
%Import lift data

fid = csvread('CL_alpha.csv',2);
a_steps = 10;
a_grid = linspace(-8,8,a_steps);

cl_points = [-4 6; 0.0 0.0];
Cl_a.sec_1(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_1(1,:) = a_grid;

cl_points = fid(1:2,1:2)';
Cl_a.sec_2(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_2(1,:) = a_grid;

cl_points = fid(1:2,3:4)';
Cl_a.sec_3(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_3(1,:) = a_grid;

cl_points = fid(1:2,5:6)';
Cl_a.sec_4(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_4(1,:) = a_grid;

cl_points = fid(1:2,7:8)';
Cl_a.sec_5(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_5(1,:) = a_grid;

cl_points = fid(1:2,9:10)';
Cl_a.sec_6(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_6(1,:) = a_grid;

cl_points = fid(1:2,11:12)';
Cl_a.sec_7(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_7(1,:) = a_grid;

cl_points = fid(1:2,13:14)';
Cl_a.sec_8(2,:) = interp1(cl_points(1,:), cl_points(2,:), a_grid,'linear','extrap');
Cl_a.sec_8(1,:) = a_grid;

Cl_a_sec_fn = fieldnames(Cl_a);

[sec_mesh,a_mesh]  = meshgrid(sec,a_grid);
CL = [Cl_a.sec_1(2,:)' Cl_a.sec_2(2,:)' Cl_a.sec_3(2,:)' Cl_a.sec_4(2,:)' ...
      Cl_a.sec_5(2,:)' Cl_a.sec_6(2,:)' Cl_a.sec_7(2,:)' Cl_a.sec_8(2,:)'];
%Import L/D data

fid = csvread('L_D_alpha.csv',2);

cd_points = [-6 6; 0.4 0.4];
CL_CD = interp1(cd_points(1,:), cd_points(2,:), a_grid,'linear','extrap');
Cd_a.sec_1 = [a_grid; CL_CD ];

cd_points = [-6 6; 0.1 0.1];
CL_CD = interp1(cd_points(1,:), cd_points(2,:), a_grid,'linear','extrap');
Cd_a.sec_2 = [a_grid; CL_CD ];

a = nonzeros(fid(:,1)');
b = nonzeros(fid(:,2)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_3 = [a_grid; 1./(CL_CD./Cl_a.sec_3(2,:)) ];

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_4 = [a_grid; 1./(CL_CD./Cl_a.sec_4(2,:)) ];


a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_5 = [a_grid; 1./(CL_CD./Cl_a.sec_5(2,:)) ];

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_6 = [a_grid; 1./(CL_CD./Cl_a.sec_6(2,:)) ];

a = nonzeros(fid(:,9)');
b = nonzeros(fid(:,10)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_7 = [a_grid; 1./(CL_CD./Cl_a.sec_7(2,:)) ];

a = nonzeros(fid(:,11)');
b = nonzeros(fid(:,12)');
CL_CD = interp1(a, b, a_grid,'pchip','extrap');
Cd_a.sec_8 = [a_grid; 1./(CL_CD./Cl_a.sec_8(2,:)) ];

% %data correction
% Cd_a.sec_4(2,7) = (Cd_a.sec_4(2,6)+Cd_a.sec_4(2,8))*0.5;
% Cd_a.sec_5(2,9) = (Cd_a.sec_5(2,8)+Cd_a.sec_5(2,10))*0.5;
% Cd_a.sec_8(2,3) = (Cd_a.sec_8(2,2)+ Cd_a.sec_8(2,4))*0.5;

Cd_a_sec_fn = fieldnames(Cd_a);
CD = [Cd_a.sec_1(2,:)' Cd_a.sec_2(2,:)' Cd_a.sec_3(2,:)' Cd_a.sec_4(2,:)' ...
      Cd_a.sec_5(2,:)' Cd_a.sec_6(2,:)' Cd_a.sec_7(2,:)' Cd_a.sec_8(2,:)'];
  

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





% Blade setting 25 deg
for V_inf=1:30
    [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_25,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_25,DiagnosticInfo);
    
    %Calculate coefficients
    CT_25(V_inf)=thrust/(rho*n^2*D^4);
    CP_25(V_inf)=power/(rho*n^3*D^5);
    J_25(V_inf)=V_inf/(n*D);
    if CT_25(V_inf) < 0
        eff_25(V_inf) = NaN;
    elseif CT_25(V_inf) == 0
        CT_25(V_inf) = NaN;
        CP_25(V_inf) = NaN;
    else
        eff_25(V_inf)=  J_25(V_inf)*(CT_25(V_inf)/CP_25(V_inf)); %J_25(V_inf)/2.0/pi*(CT_25(V_inf)/CP_25(V_inf));
    end
    
end

% Blade setting 35 deg

for V_inf=1:40
    [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_35,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_35,DiagnosticInfo);    

    CT_35(V_inf)=thrust/(rho*n^2*D^4);
    CP_35(V_inf)=power/(rho*n^3*D^5);
    J_35(V_inf)=V_inf/(n*D);
    if CT_35(V_inf) < 0
        eff_35(V_inf) = NaN;
    elseif CT_35(V_inf) == 0
        CT_35(V_inf) = NaN;
        CP_35(V_inf) = NaN;        
    else
        eff_35(V_inf)=J_35(V_inf)*(CT_35(V_inf)/CP_35(V_inf));
    end
    
end

% Blade setting 45 deg

for V_inf=1:55
        [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);  
    CT_45(V_inf)=thrust/(rho*n^2*D^4);
    CP_45(V_inf)=power/(rho*n^3*D^5);
    J_45(V_inf)=V_inf/(n*D);
    if CT_45(V_inf) < 0
        eff_45(V_inf) = NaN;
    elseif CT_45(V_inf) == 0
        CT_45(V_inf) = NaN;
        CP_45(V_inf) = NaN;  
    else
        eff_45(V_inf)=J_45(V_inf)*(CT_45(V_inf)/CP_45(V_inf));
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

set(groot,'defaultLineLineWidth',2.0)
%fit curves

f_Cp_25 = fit(Cp_j.deg_25(1,:)',Cp_j.deg_25(2,:)','smoothingspline');
f_Ct_25 = fit(Ct_j.deg_25(1,:)',Ct_j.deg_25(2,:)','smoothingspline');
f_eta_25 = fit(eta_j.deg_25(1,:)',eta_j.deg_25(2,:)','smoothingspline');
f_Cp_35 = fit(Cp_j.deg_35(1,:)',Cp_j.deg_35(2,:)','smoothingspline');
f_Ct_35 = fit(Ct_j.deg_35(1,:)',Ct_j.deg_35(2,:)','smoothingspline');
f_eta_35 = fit(eta_j.deg_35(1,:)',eta_j.deg_35(2,:)','smoothingspline');
f_Cp_45 = fit(Cp_j.deg_45(1,:)',Cp_j.deg_45(2,:)','smoothingspline');
f_Ct_45 = fit(Ct_j.deg_45(1,:)',Ct_j.deg_45(2,:)','smoothingspline');
f_eta_45 = fit(eta_j.deg_45(1,:)',eta_j.deg_45(2,:)','smoothingspline');

f1 = figure(1);
hold on 
plot(J_25,CT_25,'m',J_25,CP_25,'c')
p1 = plot(f_Ct_25,Ct_j.deg_25(1,:),Ct_j.deg_25(2,:));
p2 = plot(f_Cp_25,Cp_j.deg_25(1,:),Cp_j.deg_25(2,:));
plot(J_35,CT_35,'Color',[0 0.4470 0.7410])
plot(J_35,CP_35,'Color',[0.8500 0.3250 0.0980])
p3 = plot(f_Ct_35, Ct_j.deg_35(1,:),Ct_j.deg_35(2,:));
p4 = plot(f_Cp_35, Cp_j.deg_35(1,:),Cp_j.deg_35(2,:));
plot(J_45,CT_45,'Color',[0.9290 0.6940 0.1250])
plot(J_45,CP_45,'Color',[0.4660 0.6740 0.1880])
p5 = plot(f_Ct_45, Ct_j.deg_45(1,:),Ct_j.deg_45(2,:));
p6 = plot(f_Cp_45, Cp_j.deg_45(1,:),Cp_j.deg_45(2,:));

p1(2).Color = 'm';
p1(2).LineStyle = '--';
p2(2).Color = 'c';
p2(2).LineStyle = '--';
p3(2).Color = [0 0.4470 0.7410];
p3(2).LineStyle = '--';
p4(2).Color = [0.8500 0.3250 0.0980];
p4(2).LineStyle = '--';
p5(2).Color = [0.9290 0.6940 0.1250];
p5(2).LineStyle = '--';
p6(2).Color = [0.4660 0.6740 0.1880];
p6(2).LineStyle = '--'; 

title('Thrust and Torque Coefficients')
xlabel('Advance Ratio (J)');
ylabel('C_T, C_P');
legend('C_T computed \beta = 25 deg','C_P computed \beta = 25 deg', ...
       'C_T experimental \beta = 25 deg','C_T experimental fitted \beta = 25 deg', ...
       'C_P experimental \beta = 25 deg','C_P experimental fitted \beta = 25 deg', ...
       'C_T computed \beta = 35 deg','C_P computed \beta = 35 deg', ...
       'C_T experimental \beta = 35 deg','C_T experimental fitted \beta = 35 deg', ...
       'C_P experimental \beta = 35 deg','C_P experimental fitted \beta = 35 deg', ...       
       'C_T computed \beta = 45 deg','C_P computed \beta = 45 deg', ...
       'C_T experimental \beta = 45 deg','C_T experimental fitted \beta = 45 deg', ...
       'C_P experimental \beta = 45 deg','C_P experimental fitted \beta = 45 deg', ...
       'Location','best');
hold off


f2 = figure(2);
hold on

plot(J_25,eff_25);
plot(f_eta_25,'m',eta_j.deg_25(1,:)',eta_j.deg_25(2,:)')
plot(J_35,eff_35);
plot(f_eta_35,'m',eta_j.deg_35(1,:)',eta_j.deg_35(2,:)')
plot(J_45,eff_45);
plot(f_eta_45,'m',eta_j.deg_45(1,:)',eta_j.deg_45(2,:)')
title('Propeller Efficiency');
xlabel('Advance Ratio (J)');
ylabel('\eta');
legend('\eta computed \beta = 25 deg', '\eta experimental \beta = 25 deg', '\eta experimental fitted \beta = 25 deg', ...
       '\eta computed \beta = 35 deg', '\eta experimental \beta = 35 deg', '\eta experimental fitted \beta = 35 deg', ...
       '\eta computed \beta = 45 deg', '\eta experimental \beta = 45 deg', '\eta experimental fitted \beta = 45 deg', ...
       'Location','best');

hold off
%%
%distribution
r_R_sec = (sections + 0.5*dr);
V_loading = 35:2:45;
Va = zeros(length(V_loading),length(r_R_sec));
for i = 1:length(V_loading)
   [thrust,power,~, AxialVelocity] = SolveForFreeStreamVelocity(V_loading(i), ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);
    Va(i,:) = AxialVelocity./V_loading(i);
    J = V_loading(i)/(n*D);
    lgnd_str{i} =  strcat('J=',num2str(J,' %.2f'));
end

f3 = figure(3);
hold on 
plot(r_R_sec,Va)

xlabel('Radial location, r/R')
ylabel('V_a / V_\infty')

legend(lgnd_str,'Location','best')


%% Changing the number of blades 
N = 4;
V_inf = 45;

[thrust_baseline,power_baseline,~,~] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);
V_inf = 35;
N = 6;

   
[thrust,power,~, AxialVelocity] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);    


V_inf = 20;
N = 8;
error = 12;
while error > 10
   V_inf = V_inf + 0.5;
   [thrust,power,~, AxialVelocity] = SolveForFreeStreamVelocity(V_inf, ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    CD,CL,sec_mesh,a_mesh, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);    
if isnan(thrust)
    thrust = 1;
end
error = abs(thrust - thrust_baseline);
end





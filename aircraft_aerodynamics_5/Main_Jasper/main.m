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
dr = [diff(sec) 0.05];

%Set if pradtl tip and root correction is used
PrandtlCorrection_deg_25 = false;
PrandtlCorrection_deg_35 = false;
PrandtlCorrection_deg_45 = false;

DiagnosticInfo = false;
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





% Blade setting 25 deg
for V_inf=10:30
    [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf,sec,dr, ...
                                    tipR,rootR,pitch.deg_25,omega,N, n, r_steps, rho, ...
                                    Cl_a, Cl_a_sec_fn,Cd_a,Cd_a_sec_fn, chord, ...
                                    PrandtlCorrection_deg_25,DiagnosticInfo);
    
    %Calculate coefficients
    CT_25(V_inf)=thrust/(rho*n^2*D^4);
    CP_25(V_inf)=power/(rho*n^3*D^5);
    J_25(V_inf)=V_inf/(n*D);
    if CT_25(V_inf) < 0
        eff_25(V_inf) = 0;
    else
        eff_25(V_inf)=  J_25(V_inf)*(CT_25(V_inf)/CP_25(V_inf)); %J_25(V_inf)/2.0/pi*(CT_25(V_inf)/CP_25(V_inf));
    end
    
end

% Blade setting 35 deg

for V_inf=20:40
    [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf,sec,dr, ...
                                    tipR,rootR,pitch.deg_35,omega,N, n, r_steps, rho, ...
                                    Cl_a, Cl_a_sec_fn,Cd_a,Cd_a_sec_fn, chord, ...
                                    PrandtlCorrection_deg_35,DiagnosticInfo);    

    CT_35(V_inf)=thrust/(rho*n^2*D^4);
    CP_35(V_inf)=power/(rho*n^3*D^5);
    J_35(V_inf)=V_inf/(n*D);
    if CT_35(V_inf) < 0
        eff_35(V_inf) = 0;
    else
        eff_35(V_inf)=J_35(V_inf)*(CT_35(V_inf)/CP_35(V_inf));
    end
    
end

% Blade setting 45 deg

for V_inf=35:55
        [thrust, power, torque] = SolveForFreeStreamVelocity(V_inf,sec,dr, ...
                                    tipR,rootR,pitch.deg_45,omega,N, n, r_steps, rho, ...
                                    Cl_a, Cl_a_sec_fn,Cd_a,Cd_a_sec_fn, chord, ...
                                    PrandtlCorrection_deg_45,DiagnosticInfo);   
    CT_45(V_inf)=thrust/(rho*n^2*D^4);
    CP_45(V_inf)=power/(rho*n^3*D^5);
    J_45(V_inf)=V_inf/(n*D);
    if CT_45(V_inf) < 0
        eff_45(V_inf) = 0;
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
plot(J_35,CT_35,'m-.',J_35,CP_35,'c-.')
plot(J_45,CT_45,'m--',J_45,CP_45,'c--')
plot(f_Ct_25,'m',Ct_j.deg_25(1,:)',Ct_j.deg_25(2,:)')
plot(f_Cp_25,'c',Cp_j.deg_25(1,:)',Cp_j.deg_25(2,:)')
plot(f_Ct_35,'m-.',Ct_j.deg_35(1,:)',Ct_j.deg_35(2,:)')
plot(f_Cp_35,'c-.',Cp_j.deg_35(1,:)',Cp_j.deg_35(2,:)')
plot(f_Ct_45,'m--',Ct_j.deg_45(1,:)',Ct_j.deg_45(2,:)')
plot(f_Cp_45,'c--',Cp_j.deg_45(1,:)',Cp_j.deg_45(2,:)')

title('Thrust and Torque Coefficients')
xlabel('Advance Ratio (J)');
ylabel('C_T, C_P');
legend('C_T computed \beta = 25 deg','C_P computed \beta = 25 deg', ...
       'C_T computed \beta = 35 deg','C_P computed \beta = 35 deg', ...
       'C_T computed \beta = 45 deg','C_P computed \beta = 45 deg', ...
       'C_T experimental \beta = 25 deg', 'C_T experimental fitted \beta = 25 deg', ...
       'C_P experimental \beta = 25 deg', 'C_P experimental fitted \beta = 25 deg', ...
       'C_T experimental \beta = 35 deg', 'C_T experimental fitted \beta = 35 deg', ...
       'C_P experimental \beta = 35 deg', 'C_P experimental fitted \beta = 35 deg', ...
       'C_T experimental \beta = 45 deg', 'C_T experimental fitted \beta = 45 deg', ...
       'C_P experimental \beta = 45 deg', 'C_P experimental fitted \beta = 45 deg','Location','best');
   
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

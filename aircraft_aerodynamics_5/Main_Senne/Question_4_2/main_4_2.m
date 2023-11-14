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
for V_inf=1:50
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

%% Plotting Section
close all
 
C = [108/256, 194/256, 74/256;237/256, 104/256, 66/256; 0, 166/256, 214/256; 224/256, 60/256, 49/256];

%%%=============================== CT ==================================%%%

f1 = figure(1);
f1.Position = [200 200 1050 450];
hold on; box on; grid minor;
h=plot(f_Ct_25,'--',Ct_j.deg_25(1,:)',Ct_j.deg_25(2,:)');
set(h,'color',C(1,:));
h(1).Color=C(1,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_Ct_35,'--',Ct_j.deg_35(1,:)',Ct_j.deg_35(2,:)')';
set(h,'color',C(2,:));
h(1).Color=C(2,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_Ct_45,'--',Ct_j.deg_45(1,:)',Ct_j.deg_45(2,:)');
set(h,'color',C(3,:));
h(1).Color=C(3,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;

plot(J_25,CT_25,'-','color',C(1,:),'LineWidth',2.5)
plot(J_35,CT_35,'-','color',C(2,:),'LineWidth',2.5)
plot(J_45,CT_45,'-','color',C(3,:),'LineWidth',2.5)



ylim([-0.05 0.3])
% title('Thrust Coefficients $C_T$', 'Interpreter', 'Latex')
xlabel('Advance Ratio $J$ [-]', 'Interpreter','Latex');
ylabel('$C_T$ [-]', 'Interpreter', 'latex');
legend('$\beta = 25 ^{\circ}$, Experimental Data', ...
       '$\beta = 25 ^{\circ}$, Fitted Data', ...
       '$\beta = 35 ^{\circ}$, Experimental Data', ...
       '$\beta = 35 ^{\circ}$, Fitted Data', ...
       '$\beta = 45 ^{\circ}$, Experimental Data', ...
       '$\beta = 45 ^{\circ}$, Fitted Data', ...
       '$\beta = 25 ^{\circ}$, Blade Element Method', ...
       '$\beta = 35 ^{\circ}$, Blade Element Method', ...
       '$\beta = 45 ^{\circ}$, Blade Element method', ...
       'Interpreter','latex', 'Location','NorthEast');
hold off

%%%=============================== CP ==================================%%%

f2 = figure(2);
f2.Position = [300 300 1050 450];
hold on; box on; grid minor;
h=plot(f_Cp_25,'--',Cp_j.deg_25(1,:)',Cp_j.deg_25(2,:)');
set(h,'color',C(1,:));
h(1).Color=C(1,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_Cp_35,'--',Cp_j.deg_35(1,:)',Cp_j.deg_35(2,:)');
set(h,'color',C(2,:));
h(1).Color=C(2,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_Cp_45,'--',Cp_j.deg_45(1,:)',Cp_j.deg_45(2,:)');
set(h,'color',C(3,:));
h(1).Color=C(3,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;

plot(J_25,CP_25,'-','color',C(1,:),'LineWidth',2.5)
plot(J_35,CP_35,'-','color',C(2,:),'LineWidth',2.5)
plot(J_45,CP_45,'-','color',C(3,:),'LineWidth',2.5)

% title('Torque Coefficients $C_P$', 'Interpreter', 'latex')
xlabel('Advance Ratio $J$ [-]', 'Interpreter', 'latex');
ylabel('$C_P$ [-]', 'Interpreter', 'latex');
legend('$\beta = 25 ^{\circ}$, Experimental Data', ...
       '$\beta = 25 ^{\circ}$, Fitted Data', ...
       '$\beta = 35 ^{\circ}$, Experimental Data', ...
       '$\beta = 35 ^{\circ}$, Fitted Data', ...
       '$\beta = 45 ^{\circ}$, Experimental Data', ...
       '$\beta = 45 ^{\circ}$, Fitted Data', ...
       '$\beta = 25 ^{\circ}$, Blade Element Method', ...
       '$\beta = 35 ^{\circ}$, Blade Element Method', ...
       '$\beta = 45 ^{\circ}$, Blade Element Method', ...
       'Interpreter','latex', 'Location','NorthEast');
hold off

% %%%=============================== eff =================================%%%

f3 = figure(3);
f3.Position = [400 400 1050 450];
hold on; box on; grid minor;
h=plot(f_eta_25,'--',eta_j.deg_25(1,:)',eta_j.deg_25(2,:)');
set(h,'color',C(1,:));
h(1).Color=C(1,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_eta_35,'--',eta_j.deg_35(1,:)',eta_j.deg_35(2,:)');
set(h,'color',C(2,:));
h(1).Color=C(2,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;
h=plot(f_eta_45,'--',eta_j.deg_45(1,:)',eta_j.deg_45(2,:)');
set(h,'color',C(3,:));
h(1).Color=C(3,:);
h(1).MarkerSize=10;
h(2).LineWidth=1.75;

plot(J_25,eff_25,'-','color',C(1,:),'LineWidth',2.5);
plot(J_35,eff_35,'-','color',C(2,:),'LineWidth',2.5);
plot(J_45,eff_45,'-','color',C(3,:),'LineWidth',2.5);




xlim([0 3]); ylim([0.5 0.9])



% title('Propeller Efficiency $\eta$', 'Interpreter', 'latex');
xlabel('Advance Ratio $J$ [-]', 'Interpreter', 'latex');
ylabel('$\eta$ [-]', 'Interpreter', 'latex');


legend('$\beta = 25 ^{\circ}$, Experimental Data', ...
       '$\beta = 25 ^{\circ}$, Fitted Data', ...
       '$\beta = 35 ^{\circ}$, Experimental Data', ...
       '$\beta = 35 ^{\circ}$, Fitted Data', ...
       '$\beta = 45 ^{\circ}$, Experimental Data', ...
       '$\beta = 45 ^{\circ}$, Fitted Data', ...
       '$\beta = 25 ^{\circ}$, Blade Element Method', ...
       '$\beta = 35 ^{\circ}$, Blade Element Method', ...
       '$\beta = 45 ^{\circ}$, Blade Element Method', ...
       'Interpreter','latex', 'Location','best');
   
saveas(f1,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\J_vs_C_T.png")
saveas(f2,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\J_vs_C_P.png")
saveas(f3,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\J_vs_eta.png")

   
   

clear; close all; clc

%% Assignment Variables
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

%% Settings Prandtl tip and root correction
PrandtlCorrection_deg_25 = 0;
PrandtlCorrection_deg_35 = 0;
PrandtlCorrection_deg_45 = 0;

DiagnosticInfo = false;

%% Import lift data

fid = csvread('CL_alpha.csv', 2);

Cl_a.sec_1 = [-4 6; -0.0001 0.0001];
Cl_a.sec_2 = fid(1:2,1:2)';
Cl_a.sec_3 = fid(1:2,3:4)';
Cl_a.sec_4 = fid(1:2,5:6)';
Cl_a.sec_5 = fid(1:2,7:8)';
Cl_a.sec_6 = fid(1:2,9:10)';
Cl_a.sec_7 = fid(1:2,11:12)';
Cl_a.sec_8 = fid(1:2,13:14)';
Cl_a_sec_fn = fieldnames(Cl_a);

%% Import L/D data
fid = csvread('L_D_alpha.csv', 2);

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

%% Import chord length distribution
fid = csvread('Chord_Pitch_radial.csv',2);

chord = fid(:,1:2)';
chord(2,:) = chord(2,:)*D;

%% Import pitch distribution
fid = csvread('Chord_Pitch_radial.csv',2);

a = nonzeros(fid(:,3)');
b = nonzeros(fid(:,4)');
pitch.deg_25 = [a';b'];
pitch.deg_25(2,:) = pitch.deg_25(2,:)*D;

a = nonzeros(fid(:,5)');
b = nonzeros(fid(:,6)');
pitch.deg_35 = [a';b'];
pitch.deg_35(2,:) = pitch.deg_35(2,:)*D;

a = nonzeros(fid(:,7)');
b = nonzeros(fid(:,8)');
pitch.deg_45 = [a';b'];
pitch.deg_45(2,:) = pitch.deg_45(2,:)*D;

%% Loop over free steam velocities: V_inf

V_inf_list =  linspace(10,35,41);

for j = 1:length(V_inf_list)
    
    V_inf = V_inf_list(j);
    
    [thrust, power, torque, Vax, Vtan, r_R, thrust_list, torque_list, a_list, b_list]...
     = SolveForFreeStreamVelocity(V_inf, sec, dr, ...
     tipR, rootR, pitch.deg_35, omega, N, n, r_steps, rho, ...
     Cl_a, Cl_a_sec_fn, Cd_a, Cd_a_sec_fn, chord, ...
     PrandtlCorrection_deg_35, DiagnosticInfo);
    
    % Post Processing
    results.V_inf(j)          = V_inf;
    results.a{j}              = a_list;
    results.b{j}              = b_list;
    results.Vax{j}            = Vax;
    results.Vax_over_Vinf{j}  = Vax/V_inf;
    results.Vtan{j}           = Vtan;
    results.Vtan_over_Vinf{j} = Vtan/V_inf;
    results.CT(j)             = thrust/(rho*n^2*D^4);
    results.CP(j)             = power/(rho*n^3*D^5);
    results.J(j)              = V_inf/(n*D);
    results.thrust_list{j}    = thrust_list;
    results.torque_list{j}    = torque_list;
    if results.CT(j) < 0
        results.eff(j)        = 0;
    else
        results.eff(j)        = results.J(j)*(results.CT(j)/results.CP(j)); 
    end
    
end

%% Plotting Section
close all

C = [108/256, 194/256, 74/256; ...
     237/256, 104/256, 66/256; ...
     0, 166/256, 214/256;      ...
     224/256, 60/256, 49/256];


V_inf_selected = [20, 25, 30, 35];
V_inf_ix1  = find(results.V_inf==V_inf_selected(1));
V_inf_ix2  = find(results.V_inf==V_inf_selected(2));
V_inf_ix3  = find(results.V_inf==V_inf_selected(3));
V_inf_ix4  = find(results.V_inf==V_inf_selected(4));
J_selected = [results.J(V_inf_ix1) results.J(V_inf_ix2) results.J(V_inf_ix3) results.J(V_inf_ix4)];

f1 = figure(1);
hold on; grid minor; box on;
f1.Position = [100 100 1050 450];
plot(r_R, results.Vax_over_Vinf{V_inf_ix1},'o-','Color',C(1,:),'LineWidth', 2)
plot(r_R, results.Vax_over_Vinf{V_inf_ix2},'o-','Color',C(2,:),'LineWidth', 2)
plot(r_R, results.Vax_over_Vinf{V_inf_ix3},'o-','Color',C(3,:),'LineWidth', 2)
plot(r_R, results.Vax_over_Vinf{V_inf_ix4},'o-','Color',C(4,:),'LineWidth', 2)
xlabel('$r/R$ [-]','interpreter','latex','fontsize',13)
ylabel('$V_{ax}/V_{\infty}$ [-]','interpreter','latex','fontsize',13)
legend("$V_\infty = "+V_inf_selected(1)+" m/s$, $J = " + round(J_selected(1),2) + "$",...
       "$V_\infty = "+V_inf_selected(2)+" m/s$, $J = " + round(J_selected(2),2) + "$",...
       "$V_\infty = "+V_inf_selected(3)+" m/s$, $J = " + round(J_selected(3),2) + "$",...
       "$V_\infty = "+V_inf_selected(4)+" m/s$, $J = " + round(J_selected(4),2) + "$",...
       'interpreter','latex','fontsize',12,'location','NorthWest')

f2 = figure(2);
hold on; grid minor; box on;
f2.Position = [200 200 1050 450];
plot(r_R, results.thrust_list{V_inf_ix1},'o-','Color',C(1,:),'LineWidth',2)
plot(r_R, results.thrust_list{V_inf_ix2},'o-','Color',C(2,:),'LineWidth',2)
plot(r_R, results.thrust_list{V_inf_ix3},'o-','Color',C(3,:),'LineWidth',2)
plot(r_R, results.thrust_list{V_inf_ix4},'o-','Color',C(4,:),'LineWidth',2)

xlabel('$r/R$ [-]','interpreter','latex','fontsize',13)
ylabel('$T\;[N]$','interpreter','latex','fontsize',13)
legend("$V_\infty = "+V_inf_selected(1)+" m/s$, $J = " + round(J_selected(1),2) + "$",...
       "$V_\infty = "+V_inf_selected(2)+" m/s$, $J = " + round(J_selected(2),2) + "$",...
       "$V_\infty = "+V_inf_selected(3)+" m/s$, $J = " + round(J_selected(3),2) + "$",...
       "$V_\infty = "+V_inf_selected(4)+" m/s$, $J = " + round(J_selected(4),2) + "$",...
       'interpreter','latex','fontsize',12,'location','NorthWest')
   
saveas(f1,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\r_R_vs_v_V.png")
saveas(f2,"C:\Users\ljvdo\Documents\GitHub\Aircraft_Aerodynamics\aircraft_aerodynamics_5\FigurenJasper\r_R_vs_T.png")   
   
   
   

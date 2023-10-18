%% Program Input Settings
clear; close all; clc

list1 = {'Varying Rotor Diameter (D)',...
         'Varying No. of Blades (N)'};
[indx1,tf1] = listdlg('ListString',list1);

list2 = {'Constant CT',...
         'Constant CP'};
[indx2,tf2] = listdlg('ListString',list2);

%% Settings Prandtl tip and root correction
PrandtlCorrection_deg_25 = 0;
PrandtlCorrection_deg_35 = 0;
PrandtlCorrection_deg_45 = 0;

DiagnosticInfo = false;

%%
if indx1 == 1
    D_list    = 2.5:0.05:3.5;
    N         = 4;
    loop_size = length(D_list);
elseif indx1 == 2
    D         = 3.04;
    N_list    = [1 2 3 4 5 6 7 8 9 10];
    loop_size = length(N_list);
end

for k = 1:loop_size

if indx1 == 1
    D = D_list(k);
    disp("D = "+D)
elseif indx1 == 2
    N = N_list(k);
    disp("N = "+N)
end

tipR = D/2;
rootR = 0.1*tipR;
RPM = 400;
n = RPM/60;
omega = n*pi*2;
rho = 1.225;
sec = [0.2 0.3 0.45 0.6 0.7 0.8 0.9 0.95];

r_steps = length(sec);
dr = [diff(sec) 0.05];

%% Import Lift, Drag, Chord & Pitch data

[Cl_a, Cl_a_sec_fn, Cd_a, Cd_a_sec_fn, ...
chord, pitch] = Import_Data(D);

%% Loop over free steam velocities: V_inf

V_inf_list = linspace(15,40,501);

if indx2 == 1
    CT_aim = 0.15;
elseif indx2 == 2
    CP_aim = 0.2;
end

for j = 1:length(V_inf_list)
    
    V_inf = V_inf_list(j);
    
    [thrust, power, torque, Vax, Vtan, r_R, thrust_list, torque_list, a_list, b_list]...
     = SolveForFreeStreamVelocity(V_inf, sec, dr, ...
     tipR, rootR, pitch.deg_35, omega, N, n, r_steps, rho, ...
     Cl_a, Cl_a_sec_fn, Cd_a, Cd_a_sec_fn, chord, ...
     PrandtlCorrection_deg_35, DiagnosticInfo);
    
    CT = thrust/(rho*n^2*D^4);
    CP = power/(rho*n^3*D^5);
    
    if indx2 == 1
        if CT > CT_aim-0.001 && CT < CT_aim+0.001
            % Post Processing
            results.D(k)              = D;
            results.N(k)              = N;
            results.V_inf(k)          = V_inf;
            results.a{k}              = a_list;
            results.b{k}              = b_list;
            results.Vax{k}            = Vax;
            results.Vax_over_Vinf{k}  = Vax/V_inf;
            results.Vtan{k}           = Vtan;
            results.Vtan_over_Vinf{k} = Vtan/V_inf;
            results.CT(k)             = thrust/(rho*n^2*D^4);
            results.CP(k)             = power/(rho*n^3*D^5);
            results.J(k)              = V_inf/(n*D);
            results.thrust_list{k}    = thrust_list;
            results.torque_list{k}    = torque_list;
            if results.CT(k) < 0
                results.eff(k)        = 0;
            else
                results.eff(k)        = results.J(k)*(results.CT(k)/results.CP(k)); 
            end
            disp("V_inf = " + V_inf)
        end
    
    elseif indx2 == 2
    
        if CP > CP_aim-0.001 && CP < CP_aim+0.001
            % Post Processing
            results.D(k)              = D;
            results.N(k)              = N;
            results.V_inf(k)          = V_inf;
            results.a{k}              = a_list;
            results.b{k}              = b_list;
            results.Vax{k}            = Vax;
            results.Vax_over_Vinf{k}  = Vax/V_inf;
            results.Vtan{k}           = Vtan;
            results.Vtan_over_Vinf{k} = Vtan/V_inf;
            results.CT(k)             = thrust/(rho*n^2*D^4);
            results.CP(k)             = power/(rho*n^3*D^5);
            results.J(k)              = V_inf/(n*D);
            results.thrust_list{k}    = thrust_list;
            results.torque_list{k}    = torque_list;
            if results.CT(k) < 0
                results.eff(k)        = 0;
            else
                results.eff(k)        = results.J(k)*(results.CT(k)/results.CP(k)); 
            end
            disp("V_inf = " + V_inf)
        end
    end
        

end
end

%% Save Results

if indx2 == 1
    if indx1 == 1
        save("results_D_CT_"+CT_aim+".mat", 'results')
    elseif indx1 == 2
        save("results_N_CT_"+CT_aim+".mat", 'results')
    end
elseif indx2 == 2
    if indx1 == 1
        save("results_D_CP_"+CP_aim+".mat", 'results')
    elseif indx1 == 2
        save("results_N_CP_"+CP_aim+".mat", 'results')
    end
end
   
   
   

# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 09:32:50 2021

@author: Thomas
"""

import numpy as np
import numpy.linalg as npla

from scipy import optimize

import matplotlib.pyplot as plt
import matplotlib
from functions import solve_blade_element, find_nearest, solve_blade_element_output_rad_prop
import pandas as pd
from get_polar_fit import get_polar_fit 


#input
lift_polar_fit, drag_polar_fit, alphafit = get_polar_fit()
experimental_data = pd.read_csv('designOfOptimumPropellers.csv')

#diameter of the propeller
dia=1.75;
#engine speed in RPM
RPM=2400.;
#standard sea level atmosphere density
rho=1.225;
#RPM --> revs per sec
n=RPM/60.0;
#rps --> rads per sec
omega=n*2.0*np.pi;
#dynamic viscosity
v = 1.460e-5
#pitch distance in meters. pitch = np.tan(np.deg2rad(EXAMPLE_TWIST_0.75C))*2*np.pi*0.75*R
pitch=1.36;
pitch_R = pitch/dia*2
#chord legth of blade assumed constant with radius
chord=0.1016;
chord_R=0.1016/dia*2
#thickness to chord ratio for propeller section (constant with radius)
tonc=0.12*chord;

# settings
n_steps = 20
Tip_R =  1
Root_R =  0.17
rstep_R = (Tip_R-Root_R)/n_steps
V_max= 135
n_blades = 2  #2 it only runs one setting
dia_range =1.5   #at 0 it only runs one setting 


# Notes to myself
# Velocity considered for Re; 200 [m/s]
# Kinematic viscosity at sea level of 1.460e-5 [m2/s]
# L = chord; Thus Re= 1369863

V_range = np.arange(1, V_max)
n_blades_array = np.arange(2, n_blades+1, 2)
diameter_array= np.arange(dia - dia_range, dia+dia_range+0.25, 0.25)


#memory allocation
if len(n_blades_array) >1:
    t = np.zeros((len(n_blades_array),len(V_range))) 
    q = np.zeros((len(n_blades_array),len(V_range)))
    p = np.zeros((len(n_blades_array),len(V_range)))
    J = np.zeros(len(V_range))
    eff = np.zeros((len(n_blades_array),len(V_range)))
    icheck = np.zeros(len(V_range))
elif len(diameter_array) >1:
    t = np.zeros((len(diameter_array),len(V_range))) 
    q = np.zeros((len(diameter_array),len(V_range)))
    p = np.zeros((len(diameter_array),len(V_range)))
    V_i = np.zeros((len(diameter_array),len(V_range)))
    J = np.zeros((len(diameter_array),len(V_range)))
    eff = np.zeros((len(diameter_array),len(V_range)))
    icheck = np.zeros(len(V_range))
else:
    t = np.zeros(len(V_range))
    q = np.zeros(len(V_range))
    p = np.zeros(len(V_range))
    J = np.zeros(len(V_range))
    eff = np.zeros(len(V_range))
    icheck = np.zeros(len(V_range))
    
    


R=dia/2.0; #tip radius
rstep = rstep_R * R
r_R = np.linspace(Root_R, Tip_R, n_steps)


thrust_dist=np.zeros(len(r_R))
torque_dist=np.zeros(len(r_R))


for V in V_range:
    i = V -1
    
    # initialize radial element position
    
    ## Logic to compute data for different numbers of propellers
    if len(n_blades_array) >1:
        # initialize forces
        thrust=np.zeros((len(n_blades_array), len(r_R)))
        torque=np.zeros((len(n_blades_array), len(r_R)))
        
        for j in np.arange(0, len(r_R)):
            rad = r_R[j]*R
            for k in np.arange(0, len(n_blades_array)):
                n_blade = n_blades_array[k]
                thrust[k, j], torque[k, j], itercheck = \
                    solve_blade_element(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, n_blade, v)
        t[:,i]=np.sum(thrust, axis=1)/(rho*n*n*dia*dia*dia*dia)
        q[:,i]=np.sum(torque, axis=1)/(rho*n*n*dia*dia*dia*dia*dia)
        p[:,i] = q[:,i] * 2 * np.pi
        J[i]=V/(n*dia)
        eff[:,i]=J[i]/2.0/np.pi*t[:,i]/q[:,i]
    
    ## Logic to compute data for different diameters
    elif len(diameter_array)>1:
        # initialize forces
        thrust=np.zeros((len(diameter_array), len(r_R)))
        torque=np.zeros((len(diameter_array), len(r_R)))
        Vlocal = np.zeros(len(r_R))
        for k in np.arange(0, len(diameter_array)):
            diameter = diameter_array[k]  # diameter is the changing diameter for comparison, dia the fixed diameter
            R=diameter/2.0; #tip radius
            rstep = rstep_R * R
            pitch = pitch_R * R
            chord = chord_R * R
        
            for j in np.arange(0, len(r_R)):
                rad = r_R[j]*R
                thrust[k, j], torque[k, j], itercheck = \
                    solve_blade_element(rad, R, diameter, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, 2, v)
            t[k,i]=np.sum(thrust[k,:])/(rho*n*n*diameter*diameter*diameter*diameter)
            q[k,i]=np.sum(torque[k,:])/(rho*n*n*diameter*diameter*diameter*diameter*diameter)
            V_i[k,i] = max(Vlocal)
            p[k,i] = q[k,i] * 2 * np.pi
            J[k,i]=V/(n*diameter)
            eff[k,i]=J[k,i]/2.0/np.pi*t[k,i]/q[k,i]


    else:
        thrust=np.zeros(len(r_R))
        torque=np.zeros(len(r_R))
        V0=np.zeros(len(r_R))
        V2=np.zeros(len(r_R))
        phi=np.zeros(len(r_R))
        alpha=np.zeros(len(r_R))
        for j in np.arange(0, len(r_R)):
            rad = r_R[j]*R
            thrust[j], torque[j], itercheck, V0[j], V2[j], phi[j], alpha[j] = \
                    solve_blade_element_output_rad_prop(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, 2, v)
#                thrust[j], torque[j], itercheck = \
#                    solve_blade_element(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, 2, v)
        t[i]=np.sum(thrust)/(rho*n*n*dia*dia*dia*dia)
        q[i]=np.sum(torque)/(rho*n*n*dia*dia*dia*dia*dia)
        p[i]=q[i] * 2 * np.pi
        J[i]=V/(n*dia)
        eff[i]=J[i]/2.0/np.pi*t[i]/q[i]
        ## plot radial distributions
        if V ==30:
            #Thrust cst chord exp
            fig1, ax1 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
            ax2 = ax1[0,0].twinx()
            ax1[0,0].plot(r_R, V0/V, linewidth=2, linestyle="dotted", label="Vaxi")
            ax1[0,0].plot(r_R, V2/V, linewidth=2, linestyle="dashdot", label="Vtan")
            ax2.plot(r_R, np.rad2deg(phi), "-o", label="φ")
            ax2.plot(r_R, np.rad2deg(alpha), "-x", label="α")
            ax1[0,0].set_xlabel(r"r/R $\,\,[-]$")
            ax1[0,0].set_ylabel(r"$V/V_{\infty}[-]$")
            ax2.set_ylabel(r"$angle[^{\circ}]$")
            ax1[0,0].grid(True,which="major",color="#999999",alpha=0.75)
            ax1[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
            ax1[0,0].minorticks_on()
            ax1[0,0].tick_params(which='major', length=10, width=2, direction='inout')
            ax1[0,0].tick_params(which='minor', length=5, width=2, direction='in')
            ax1[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
            ax2.legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
            
            #Force over radius
            fig9, ax9 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
            ax9[0,0].plot(r_R, thrust, linewidth=2, linestyle="dotted", label="Thrust")
            ax9[0,0].plot(r_R, torque, linewidth=2, linestyle="dashdot", label="Torque")
            ax9[0,0].set_xlabel(r"r/R $\,\,[-]$")
            ax9[0,0].set_ylabel(r"$F[N]$")
            ax9[0,0].grid(True,which="major",color="#999999",alpha=0.75)
            ax9[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
            ax9[0,0].minorticks_on()
            ax9[0,0].tick_params(which='major', length=10, width=2, direction='inout')
            ax9[0,0].tick_params(which='minor', length=5, width=2, direction='in')
            ax9[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
        
    icheck[i]=itercheck
    
if len(n_blades_array) >1:
     # same thrust settings 
    thrust_settings = [0.11, 0.07, 0.04]
    t_fixed_thrust = np.zeros((len(thrust_settings), len(n_blades_array)))
    q_fixed_thrust = np.zeros((len(thrust_settings), len(n_blades_array)))
    p_fixed_thrust = np.zeros((len(thrust_settings), len(n_blades_array)))
    J_fixed_thrust = np.zeros((len(thrust_settings), len(n_blades_array)))
    eff_fixed_thrust = np.zeros((len(thrust_settings), len(n_blades_array)))
    
    for i in np.arange(0, len(thrust_settings)):
        fixed_thrust = thrust_settings[i]
        for k in np.arange(0, len(n_blades_array)):
            t_fixed_thrust[i,k] = find_nearest(t[k,:], fixed_thrust)
            idx = np.where(t[k,:] == t_fixed_thrust[i,k])[0][0]
            q_fixed_thrust[i,k] = q[k, idx]
            p_fixed_thrust[i,k] = p[k, idx]
            J_fixed_thrust[i,k] = J[idx]
            eff_fixed_thrust[i,k] = eff[k, idx]
            
    
    # same power settings 
    power_settings = [0.05, 0.031, 0.02]
    t_fixed_power = np.zeros((len(power_settings), len(n_blades_array)))
    q_fixed_power = np.zeros((len(power_settings), len(n_blades_array)))
    p_fixed_power = np.zeros((len(power_settings), len(n_blades_array)))
    J_fixed_power = np.zeros((len(power_settings), len(n_blades_array)))
    eff_fixed_power = np.zeros((len(power_settings), len(n_blades_array)))
    
    for i in np.arange(0, len(power_settings)):
        fixed_power = power_settings[i]
        for k in np.arange(0, len(n_blades_array)):
            p_fixed_power[i,k] = find_nearest(p[k,:], fixed_power)
            idx = np.where(p[k,:] == p_fixed_power[i,k])[0][0]
            q_fixed_power[i,k] = q[k, idx]
            t_fixed_power[i,k] = t[k, idx]
            J_fixed_power[i,k] = J[idx]
            eff_fixed_power[i,k] = eff[k, idx]
elif len(diameter_array)>1:
    # same thrust settings 
    thrust_settings = [0.11, 0.07, 0.04]
    t_fixed_thrust = np.zeros((len(thrust_settings), len(diameter_array)))
    q_fixed_thrust = np.zeros((len(thrust_settings), len(diameter_array)))
    p_fixed_thrust = np.zeros((len(thrust_settings), len(diameter_array)))
    J_fixed_thrust = np.zeros((len(thrust_settings), len(diameter_array)))
    eff_fixed_thrust = np.zeros((len(thrust_settings), len(diameter_array)))
    
    for i in np.arange(0, len(thrust_settings)):
        fixed_thrust = thrust_settings[i]
        for k in np.arange(0, len(diameter_array)):
            t_fixed_thrust[i,k] = find_nearest(t[k,:], fixed_thrust)
            idx = np.where(t[k,:] == t_fixed_thrust[i,k])[0][0]
            q_fixed_thrust[i,k] = q[k, idx]
            p_fixed_thrust[i,k] = p[k, idx]
            J_fixed_thrust[i,k] = J[k, idx]
            eff_fixed_thrust[i,k] = eff[k, idx]
            
    
    # same power settings 
    power_settings = [0.05, 0.031, 0.02]
    t_fixed_power = np.zeros((len(power_settings), len(diameter_array)))
    q_fixed_power = np.zeros((len(power_settings), len(diameter_array)))
    p_fixed_power = np.zeros((len(power_settings), len(diameter_array)))
    J_fixed_power = np.zeros((len(power_settings), len(diameter_array)))
    eff_fixed_power = np.zeros((len(power_settings), len(diameter_array)))
    
    for i in np.arange(0, len(power_settings)):
        fixed_power = power_settings[i]
        for k in np.arange(0, len(diameter_array)):
            p_fixed_power[i,k] = find_nearest(p[k,:], fixed_power)
            idx = np.where(p[k,:] == p_fixed_power[i,k])[0][0]
            q_fixed_power[i,k] = q[k, idx]
            t_fixed_power[i,k] = t[k, idx]
            J_fixed_power[i,k] = J[k, idx]
            eff_fixed_power[i,k] = eff[k, idx]
    
    # =============================================================================
    # # keep thrust constant
    # for J in [0.3, 0.5, 0.8]:
    #     for j in np.arange(0, len(r_R)):
    #         rad = r_R[j]*R
    #         thrust_dist[j], torque_dist[j], itercheck = \
    #                 solve_blade_element(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, 2, v)
    #     for k in np.arange(0, len(diameter_array)):
    #         diameter = diameter_array[k]  # diameter is the changing diameter for comparison, dia the fixed diameter
    #         V_inf = J*n*diameter    #Compute velocity for this diameter
    #         R=diameter/2.0; #tip radius
    #         rstep = rstep_R * R
    #         pitch = pitch_R * R
    #         chord = chord_R * R
    #     
    #         for j in np.arange(0, len(r_R)):
    #             rad = r_R[j]*R
    #             thrust[k, j], torque[k, j], itercheck, pitch[k,j] = \
    #                 solve_blade_element_cst_thrust(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, 2, v, thrust_dist[j])
    #         t[i,k]=np.sum(thrust[k,:])/(rho*n*n*diameter*diameter*diameter*diameter)
    #         q[i,k]=np.sum(torque[k,:])/(rho*n*n*diameter*diameter*diameter*diameter*diameter)
    #         V_i[i,k] = max(Vlocal)
    #         p[i,k] = q[k,i] * 2 * np.pi
    #         J[i,k]=V/(n*diameter)
    #         eff[i,k]=J[k,i]/2.0/np.pi*t[k,i]/q[k,i]
    # =============================================================================
        
#%% Reverse calculating     
J_CT_Ccst = experimental_data['J_CT_Ccst'] #no direct plotting or manipulation, so there we go
t_Ccst =  experimental_data['CT_Ccst']
eff_Ccst =  experimental_data['eff_Ccst']
J_eff_Ccst = experimental_data['J_eff_Ccst']
eff_in_t_coordinates = np.interp(J_CT_Ccst, J_eff_Ccst, eff_Ccst)
q_exp = J_CT_Ccst*t_Ccst/(2*np.pi*eff_in_t_coordinates)

#%% Plotting
#plot settings
## Define text sizes for **SAVED** pictures (texpsize -- text export size)
texpsize= [26,28,30]

## Graphing Parameters
SMALL_SIZE  = texpsize[0]
MEDIUM_SIZE = texpsize[1]
BIGGER_SIZE = texpsize[2]

plt.style.use('grayscale')
plt.rc('font', size=SMALL_SIZE, family='serif')    ## controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)                ## fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)                ## fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)               ## fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)               ## legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)             ## fontsize of the figure title
plt.rc('text', usetex=False)
matplotlib.rcParams['lines.linewidth']  = 1.5
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['axes.facecolor']   = 'white'
matplotlib.rcParams["legend.fancybox"]  = False


if len(diameter_array)>1:
    Jmax=J[-2,np.argmax(q[-2,:]<0)]
    Tmax=max(t[-1,:]);
    qmax=max(q[-1,:]);
    pmax=max(p[-1,:]);
elif len(n_blades_array) >1:
    Tmax=max(t[-1,:]);
    qmax=max(q[-1,:]);
    pmax=max(p[-1,:]);
    Tmax_2blades = max(t[0,:])
    Jmax=J[np.argmax(eff[0,:]<0)+1];
else:
    Tmax=max(t);
    qmax=max(q);
    pmax=max(p);
    Jmax=J[np.argmax(eff<0)+1];


if len(n_blades_array) >1:
    ## Thrust varying number of propellers
    fig1, ax1 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax1[0,0].plot(J, t[0,:], linewidth=2, label="".join(("nb=", str(n_blades_array[0]))))
    ax1[0,0].plot(J, t[3,:], linewidth=2, linestyle="dashed", label="".join(("nb=", str(n_blades_array[3]))))
    ax1[0,0].plot(J, t[-1,:], linewidth=2, linestyle="dashdot", label="".join(("nb=", str(n_blades_array[-1]))))
    ax1[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax1[0,0].set_ylabel(r"$C_T[-]$")
    ax1[0,0].set_xlim(0, Jmax)
    ax1[0,0].set_ylim(0, 1.1*Tmax)
    ax1[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax1[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax1[0,0].minorticks_on()
    ax1[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax1[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax1[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    ## Thrust varying number of propellers
    fig5, ax5 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax5[0,0].plot(J, q[1,:], linewidth=2, label="".join(("nb=", str(n_blades_array[0]))))
    ax5[0,0].plot(J, q[3,:], linewidth=2,  linestyle="dashed", label="".join(("nb=", str(n_blades_array[3]))))
    ax5[0,0].plot(J, q[-1,:], linewidth=2,  linestyle="dashdot", label="".join(("nb=", str(n_blades_array[-1]))))
    ax5[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax5[0,0].set_ylabel(r"$C_q[-]$")
    ax5[0,0].set_xlim(0, Jmax)
    ax5[0,0].set_ylim(0, 1.1*qmax)
    ax5[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax5[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax5[0,0].minorticks_on()
    ax5[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax5[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax5[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')    
    
    fig6, ax6 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax6[0,0].plot(J, p[1,:], linewidth=2, label="".join(("nb=", str(n_blades_array[0]))))
    ax6[0,0].plot(J, p[3,:], linewidth=2,  linestyle="dashed", label="".join(("nb=", str(n_blades_array[3]))))
    ax6[0,0].plot(J, p[-1,:], linewidth=2,  linestyle="dashdot", label="".join(("nb=", str(n_blades_array[-1]))))
    ax6[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax6[0,0].set_ylabel(r"$C_p[-]$")
    ax6[0,0].set_xlim(0, Jmax)
    ax6[0,0].set_ylim(0, 1.1*pmax)
    ax6[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax6[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax6[0,0].minorticks_on()
    ax6[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax6[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax6[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    ## Efficiency constant chord
    fig2, ax2 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax2[0,0].plot(J[:], eff[0,:], linewidth=2, label="2 blades")
    ax2[0,0].plot(J[:], eff[3,:], linestyle="dashed", linewidth=2, label="6 blades")
    ax2[0,0].plot(J[:], eff[-1,:], linestyle="dashdot", linewidth=2, label="12 blades")
    ax2[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax2[0,0].set_ylabel(r"Propeller efficiency $[-]$")
    ax2[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax2[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax2[0,0].minorticks_on()
    ax2[0,0].set_xlim(0, Jmax)
    ax2[0,0].set_ylim(0, 1)
    ax2[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax2[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax2[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    #Thrust cst chord exp
    J_start = np.argmax(J>0.4)
    fig3, ax3 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax4 = ax3[0,0].twinx()
    ax3[0,0].plot(J[J_start:], t[0,J_start:], linewidth=2, label="Ct BEM")
    ax3[0,0].plot(J[J_start:], q[0,J_start:], linewidth=2, label="Cq BEM")
#    ax3[0,0].plot(J[J_start:], p[0,J_start:], linewidth=2, label="Cp BEM")
    ax3[0,0].plot(J_CT_Ccst, t_Ccst,"-^", label="Ct Exp")
    ax3[0,0].plot(J_CT_Ccst, q_exp,"-o", label="Cq Exp")
    ax4.plot(J[J_start:69], eff[0,J_start:69], linewidth=2, label="η BEM")
    ax4.plot(J_eff_Ccst, eff_Ccst, "-x", label="η Exp")
    ax3[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax3[0,0].set_ylabel(r"$[-]$")
    ax4.set_ylabel(r"Propeller efficiency $[-]$")
    ax4.set_ylim(0, 1)
    ax3[0,0].set_xlim(J[J_start], Jmax)
    ax3[0,0].set_ylim(0, 1.1*Tmax_2blades)
    ax3[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax3[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax3[0,0].minorticks_on()
    ax3[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax3[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax3[0,0].legend(loc=2, framealpha=1.0).get_frame().set_edgecolor('k')
    ax4.legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
        
    fig7, ax7 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax7[0,0].plot(n_blades_array, p_fixed_thrust[0,:], linewidth=2, linestyle="dashed", label="CT=0.118")
    ax7[0,0].plot(n_blades_array, p_fixed_thrust[1,:], linewidth=2,  linestyle="dotted", label="CT=0.094")
    ax7[0,0].plot(n_blades_array, p_fixed_thrust[2,:], linewidth=2,  linestyle="dashdot", label="CT=0.037")
    ax7[0,0].set_xlabel(r"Number of blades $\,\,[-]$")
    ax7[0,0].set_ylabel(r"$C_p[-]$")
    #ax7[0,0].set_ylim(0, 1.1*pmax)
    ax7[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax7[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax7[0,0].minorticks_on()
    ax7[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax7[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax7[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k') 
    
    fig8, ax8 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax8[0,0].plot(n_blades_array, t_fixed_power[0,:], linewidth=2, linestyle="dashed", label="CP=0.05")
    ax8[0,0].plot(n_blades_array, t_fixed_power[1,:], linewidth=2,  linestyle="dotted", label="CP=0.03")
    ax8[0,0].plot(n_blades_array, t_fixed_power[2,:], linewidth=2,  linestyle="dashdot", label="CP=0.02")
    ax8[0,0].set_xlabel(r"Number of blades $\,\,[-]$")
    ax8[0,0].set_ylabel(r"$C_T[-]$")
    #ax8[0,0].set_ylim(0, 1.1*pmax)
    ax8[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax8[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax8[0,0].minorticks_on()
    ax8[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax8[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax8[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')  
    
    fig9, ax9 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax9[0,0].plot(n_blades_array, eff_fixed_thrust[0,:], linewidth=2, linestyle="dashed", label="CT=0.118")
    ax9[0,0].plot(n_blades_array, eff_fixed_thrust[1,:], linewidth=2,  linestyle="dotted", label="CT=0.094")
    ax9[0,0].plot(n_blades_array, eff_fixed_thrust[2,:], linewidth=2,  linestyle="dashdot", label="CT=0.037")
    ax9[0,0].set_xlabel(r"Number of blades $\,\,[-]$")
    ax9[0,0].set_ylabel(r"$\eta[-]$")
    ax9[0,0].set_ylim(0, 1)
    ax9[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax9[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax9[0,0].minorticks_on()
    ax9[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax9[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax9[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')   
    
    fig10, ax10 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax10[0,0].plot(n_blades_array, eff_fixed_power[0,:], linewidth=2, linestyle="dashed", label="CP=0.05")
    ax10[0,0].plot(n_blades_array, eff_fixed_power[1,:], linewidth=2,  linestyle="dotted", label="CP=0.03")
    ax10[0,0].plot(n_blades_array, eff_fixed_power[2,:], linewidth=2,  linestyle="dashdot", label="CP=0.02")
    ax10[0,0].set_xlabel(r"Number of blades $\,\,[-]$")
    ax10[0,0].set_ylabel(r"$\eta[-]$")
    ax10[0,0].set_ylim(0, 1)
    ax10[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax10[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax10[0,0].minorticks_on()
    ax10[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax10[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax10[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')   


elif len(diameter_array)>1:
    ## Thrust varying number of propellers
    fig6, ax6 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    #ax6[0,0].plot(J[0,:], t[0,:], "-o", linewidth=2, label="".join(("ø=", str(diameter_array[0]), ' [m]')))
    ax6[0,0].plot(J[1,:], t[1,:], linewidth=2, linestyle="dashed",label="".join(("ø=", str(diameter_array[1]), ' [m]')))
    ax6[0,0].plot(J[2,:], t[2,:], linewidth=2, linestyle="dashdot",label="".join(("ø=", str(diameter_array[2]), ' [m]')))
    ax6[0,0].plot(J[3,:], t[3,:], linewidth=2, linestyle="dotted",label="".join(("ø=", str(diameter_array[3]), ' [m]')))
    ax6[0,0].plot(J[-1,:], t[-1,:], linewidth=2, linestyle="solid", label="".join(("ø=", str(diameter_array[-1]), ' [m]')))
    ax6[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax6[0,0].set_ylabel(r"$C_T[-]$")
    ax6[0,0].set_xlim(0, Jmax)
    ax6[0,0].set_ylim(0, 1.1*Tmax)
    ax6[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax6[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax6[0,0].minorticks_on()
    ax6[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax6[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax6[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
        
    fig5, ax5 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax5[0,0].plot(J[1,:], p[1,:], linewidth=2, linestyle="dashed", label="".join(("ø=", str(diameter_array[1]), '[m]')))
    ax5[0,0].plot(J[2,:], p[2,:], linewidth=2,  linestyle="dotted", label="".join(("ø=", str(diameter_array[2]), '[m]')))
    ax5[0,0].plot(J[3,:], p[3,:], linewidth=2,  linestyle="dashdot", label="".join(("ø=", str(diameter_array[3]), '[m]')))
    ax5[0,0].plot(J[-1,:], p[-1,:], linewidth=2,  linestyle="solid", label="".join(("ø=", str(diameter_array[-1]), '[m]')))
    ax5[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax5[0,0].set_ylabel(r"$C_p[-]$")
    ax5[0,0].set_xlim(0, Jmax)
    ax5[0,0].set_ylim(0, 1.1*pmax)
    ax5[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax5[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax5[0,0].minorticks_on()
    ax5[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax5[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax5[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    fig3, ax3 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax3[0,0].plot(diameter_array, p_fixed_thrust[0,:], linewidth=2, linestyle="dashed", label="CT=0.118")
    ax3[0,0].plot(diameter_array, p_fixed_thrust[1,:], linewidth=2,  linestyle="dotted", label="CT=0.094")
    ax3[0,0].plot(diameter_array, p_fixed_thrust[2,:], linewidth=2,  linestyle="dashdot", label="CT=0.037")
    ax3[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax3[0,0].set_ylabel(r"$C_P[-]$")
    #ax3[0,0].set_xlim(0, Jmax)
    #ax3[0,0].set_ylim(0, 1.1*pmax)
    ax3[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax3[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax3[0,0].minorticks_on()
    ax3[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax3[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax3[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k') 
    
    fig4, ax4 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax4[0,0].plot(diameter_array, p_fixed_power[0,:], linewidth=2, linestyle="dashed", label="CP=0.05")
    ax4[0,0].plot(diameter_array, p_fixed_power[1,:], linewidth=2,  linestyle="dotted", label="CP=0.03")
    ax4[0,0].plot(diameter_array, p_fixed_power[2,:], linewidth=2,  linestyle="dashdot", label="CP=0.02")
    ax4[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax4[0,0].set_ylabel(r"$C_p[-]$")
    #ax4[0,0].set_xlim(0, Jmax)
    #ax4[0,0].set_ylim(0, 1.1*pmax)
    ax4[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax4[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax4[0,0].minorticks_on()
    ax4[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax4[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax4[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')  
    
    fig13, ax13 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax13[0,0].plot(diameter_array, t_fixed_thrust[0,:], linewidth=2, linestyle="dashed", label="CT=0.118")
    ax13[0,0].plot(diameter_array, t_fixed_thrust[1,:], linewidth=2,  linestyle="dotted", label="CT=0.094")
    ax13[0,0].plot(diameter_array, t_fixed_thrust[2,:], linewidth=2,  linestyle="dashdot", label="CT=0.037")
    ax13[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax13[0,0].set_ylabel(r"$C_t[-]$")
    #ax13[0,0].set_xlim(0, Jmax)
    #ax13[0,0].set_ylim(0, 1.1*pmax)
    ax13[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax13[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax13[0,0].minorticks_on()
    ax13[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax13[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax13[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k') 
    
    fig14, ax14 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    #ax14[0,0].plot(diameter_array, t_fixed_power[0,:], linewidth=2, linestyle="dashed", label="CP=0.05")
    ax14[0,0].plot(diameter_array, t_fixed_power[1,:], linewidth=2,  linestyle="dotted", label="CP=0.03")
    ax14[0,0].plot(diameter_array, t_fixed_power[2,:], linewidth=2,  linestyle="dashdot", label="CP=0.02")
    ax14[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax14[0,0].set_ylabel(r"$C_t[-]$")
    #ax14[0,0].set_xlim(0, Jmax)
    #ax14[0,0].set_ylim(0, 1.1*pmax)
    ax14[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax14[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax14[0,0].minorticks_on()
    ax14[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax14[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax14[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')  
    
    ## Thrust varying number of propellers
    fig7, ax7 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    #ax7[0,0].plot(J[0,:], q[0,:], "-o", linewidth=2, label="".join(("ø=", str(diameter_array[0]), ' [m]')))
    ax7[0,0].plot(J[1,:], q[1,:], linewidth=2,  linestyle="dashed", label="".join(("ø=", str(diameter_array[1]), ' [m]')))
    ax7[0,0].plot(J[2,:], q[2,:], linewidth=2,  linestyle="dotted", label="".join(("ø=", str(diameter_array[2]), ' [m]')))
    ax7[0,0].plot(J[3,:], q[3,:], linewidth=2,  linestyle="dashdot", label="".join(("ø=", str(diameter_array[3]), ' [m]')))
    ax7[0,0].plot(J[-1,:], q[-1,:], linewidth=2,  linestyle="solid", label="".join(("ø=", str(diameter_array[-1]), ' [m]')))
    ax7[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax7[0,0].set_ylabel(r"$C_q[-]$")
    ax7[0,0].set_xlim(0, Jmax)
    ax7[0,0].set_ylim(0, 1.1*qmax)
    ax7[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax7[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax7[0,0].minorticks_on()
    ax7[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax7[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax7[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    ## Efficiency constant chord
    fig8, ax8 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax8[0,0].plot(J[1,:np.argmax(eff[1,:]<0)+1], eff[1,:np.argmax(eff[1,:]<0)+1], linestyle="dashed", linewidth=2, label="".join(("ø=", str(diameter_array[1]), ' [m]')))
    ax8[0,0].plot(J[2,:np.argmax(eff[2,:]<0)+1], eff[2,:np.argmax(eff[2,:]<0)+1], linestyle="dotted", linewidth=2, label="".join(("ø=", str(diameter_array[2]), ' [m]')))
    ax8[0,0].plot(J[3,:np.argmax(eff[3,:]<0)+1], eff[3,:np.argmax(eff[3,:]<0)+1], linestyle="dashdot", linewidth=2, label="".join(("ø=", str(diameter_array[3]), ' [m]')))
    ax8[0,0].plot(J[-1,:np.argmax(eff[-1,:]<0)+1], eff[-1,:np.argmax(eff[-1,:]<0)+1], linestyle="solid", linewidth=2, label="".join(("ø=", str(diameter_array[-1]), ' [m]')))
    ax8[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax8[0,0].set_ylabel(r"Propeller efficiency $[-]$")
    ax8[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax8[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax8[0,0].minorticks_on()
    ax8[0,0].set_xlim(0, Jmax)
    ax8[0,0].set_ylim(0, 1)
    ax8[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax8[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax8[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    
    fig1, ax1 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax1[0,0].plot(diameter_array, eff_fixed_thrust[0,:], linewidth=2, linestyle="dashed", label="CT=0.118")
    ax1[0,0].plot(diameter_array, eff_fixed_thrust[1,:], linewidth=2,  linestyle="dotted", label="CT=0.094")
    ax1[0,0].plot(diameter_array, eff_fixed_thrust[2,:], linewidth=2,  linestyle="dashdot", label="CT=0.037")
    ax1[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax1[0,0].set_ylabel(r"$\eta[-]$")
    ax1[0,0].set_ylim(0, 1)
    ax1[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax1[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax1[0,0].minorticks_on()
    ax1[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax1[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax1[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')   
    
    fig2, ax2 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    #ax2[0,0].plot(diameter_array, eff_fixed_power[0,:], linewidth=2, linestyle="dashed", label="CP=0.05")
    ax2[0,0].plot(diameter_array, eff_fixed_power[1,:], linewidth=2,  linestyle="dotted", label="CP=0.03")
    ax2[0,0].plot(diameter_array, eff_fixed_power[2,:], linewidth=2,  linestyle="dashdot", label="CP=0.02")
    ax2[0,0].set_xlabel(r"Diameter $\,\,[m]$")
    ax2[0,0].set_ylabel(r"$\eta[-]$")
    ax2[0,0].set_ylim(0, 1)
    ax2[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax2[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax2[0,0].minorticks_on()
    ax2[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax2[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax2[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')   
else:
    #Thrust cst chord exp
    J_start = np.argmax(J>0.4)
    fig3, ax3 = plt.subplots(1,1,squeeze=False,figsize=(16,9))
    ax4 = ax3[0,0].twinx()
    ax3[0,0].plot(J, t, linewidth=2, label="Ct BEM")
    ax3[0,0].plot(J, q, linewidth=2, label="Cq BEM")
    ax3[0,0].plot(J_CT_Ccst, t_Ccst,"-^", label="Ct Exp")
    ax3[0,0].plot(J_CT_Ccst, q_exp,"-o", label="Cq Exp")
    ax4.plot(J, eff, linewidth=2, label="η BEM")
    ax4.plot(J_eff_Ccst, eff_Ccst, "-x", label="η Exp")
    ax3[0,0].set_xlabel(r"Advance ratio(J) $\,\,[-]$")
    ax3[0,0].set_ylabel(r"$[-]$")
    ax4.set_ylabel(r"Propeller efficiency $[-]$")
    ax4.set_ylim(0, 1)
    ax3[0,0].set_xlim(J[J_start], Jmax)
    ax3[0,0].set_ylim(0, 1.1*Tmax)
    ax3[0,0].grid(True,which="major",color="#999999",alpha=0.75)
    ax3[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.50)
    ax3[0,0].minorticks_on()
    ax3[0,0].tick_params(which='major', length=10, width=2, direction='inout')
    ax3[0,0].tick_params(which='minor', length=5, width=2, direction='in')
    ax4.legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
    ax3[0,0].legend(loc=0, framealpha=1.0).get_frame().set_edgecolor('k')
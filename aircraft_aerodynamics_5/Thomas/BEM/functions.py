# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 09:31:49 2021

@author: thoma
"""

import numpy as np
import numba as nb
    

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def PrandtlTipRootCorrection(r_R, rootradius_R, tipradius_R, TSR, NBlades, axial_induction):
    """
    This function calcualte steh combined tip and root Prandtl correction at agiven radial position 'r_R' (non-dimensioned by rotor radius), 
    given a root and tip radius (also non-dimensioned), a tip speed ratio TSR, the number lf blades NBlades and the axial induction factor
    """
    temp1 = -NBlades/2*(tipradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Ftip = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Ftip[np.isnan(Ftip)] = 0
    temp1 = NBlades/2*(rootradius_R-r_R)/r_R*np.sqrt( 1+ ((TSR*r_R)**2)/((1-axial_induction)**2))
    Froot = np.array(2/np.pi*np.arccos(np.exp(temp1)))
    Froot[np.isnan(Froot)] = 0
    return Froot*Ftip

def solve_blade_element_output_rad_prop(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, n_blades, v):    
    #calculate local blade element setting angle
    theta=np.arctan(pitch/2/np.pi/rad)
    #calculate solidity
    sigma=n_blades*chord/2.0/np.pi/rad
    #guess initial values of inflow and swirl factor
    axi_induction=0.1
    swirl=0.01
    #set logical variable to control iteration
    finished=False
    #set iteration count and check flag
    total_sum=1
    itercheck=0
    while not finished:
        #axial velocity
        V0=V*(1+axi_induction)
        #disk plane velocity
        V2=omega*rad*(1-swirl)
        #flow angle
        phi=np.arctan2(V0,V2)
        #blade angle of attack
        alpha=theta-phi
        #docal velocity at blade
        Vlocal=np.sqrt(V0*V0+V2*V2)
        Relocal = Vlocal * chord/v
        try:
            # lift coefficient
            Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
            #drag coefficient
            Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
        except: 
            print('ERROR: Reynolds number ', Relocal, ' not within range')
            print('For diameter: ', dia, ', with V: ', V, ', and n_blades: ', n_blades)
            if Relocal>35e6:
                Relocal = 354.99e6
                Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
                #drag coefficient
                Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
            else: 
                Relocal = 1.2e4
                Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
                #drag coefficient
                Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
                
        #thrust grading
        DtDr=0.5*rho*Vlocal*Vlocal*n_blades*chord*(Cl*np.cos(phi)-Cd*np.sin(phi))
        #torque grading
        DqDr=0.5*rho*Vlocal*Vlocal*n_blades*chord*rad*(Cd*np.cos(phi)+Cl*np.sin(phi))
        #momentum check on inflow and swirl factors
        tem1=DtDr/(4.0*np.pi*rad*rho*V*V*(1+axi_induction))
        tem2=DqDr/(4.0*np.pi*rad*rad*rad*rho*V*(1+axi_induction)*omega)

        Prandtl_correction = PrandtlTipRootCorrection(rad/R, 0.1, 1, omega*R/V0, n_blades, tem1)
        if (Prandtl_correction < 1e-3): 
           Prandtl_correction = 1e-3
          
        #stabilise iteration
        anew=0.5*(axi_induction+tem1*Prandtl_correction)
        bnew=0.5*(swirl+tem2*Prandtl_correction)

        #check for convergence
        if (abs(anew-axi_induction)<1.0e-5):
            if (abs(bnew-swirl)<1.0e-5):
                finished=True
                
        axi_induction=anew
        swirl=bnew
        #increment iteration count
        total_sum=total_sum+1
        #check to see if iteration stuck
        if (total_sum>500):
            finished=True
            itercheck=1
            print('iteration failed')
        thrust=DtDr*rstep
        torque=DqDr*rstep
    return [thrust, torque, itercheck, V0, V2, phi, alpha]    

def solve_blade_element(rad, R, dia, n, rstep, pitch, chord, omega, rho, V, lift_polar_fit, drag_polar_fit, alphafit, n_blades, v):
    
    
    #calculate local blade element setting angle
    theta=np.arctan(pitch/2/np.pi/rad)
    #calculate solidity
    sigma=n_blades*chord/2.0/np.pi/rad
    #guess initial values of inflow and swirl factor
    axi_induction=0.1
    swirl=0.01
    #set logical variable to control iteration
    finished=False
    #set iteration count and check flag
    total_sum=1
    itercheck=0
    while not finished:
        #axial velocity
        V0=V*(1+axi_induction)
        #disk plane velocity
        V2=omega*rad*(1-swirl)
        #flow angle
        phi=np.arctan2(V0,V2)
        #blade angle of attack
        alpha=theta-phi
        #docal velocity at blade
        Vlocal=np.sqrt(V0*V0+V2*V2)
        Relocal = Vlocal * chord/v
        try:
            # lift coefficient
            Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
            #drag coefficient
            Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
        except: 
            print('ERROR: Reynolds number ', Relocal, ' not within range')
            print('For diameter: ', dia, ', with V: ', V, ', and n_blades: ', n_blades)
            if Relocal>35e6:
                Relocal = 354.99e6
                Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
                #drag coefficient
                Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
            else: 
                Relocal = 1.2e4
                Cl = np.interp(np.rad2deg(alpha), alphafit(Relocal), lift_polar_fit(Relocal))
                #drag coefficient
                Cd = np.interp(np.rad2deg(alpha), alphafit(Relocal), drag_polar_fit(Relocal))
                
        #thrust grading
        DtDr=0.5*rho*Vlocal*Vlocal*n_blades*chord*(Cl*np.cos(phi)-Cd*np.sin(phi))
        #torque grading
        DqDr=0.5*rho*Vlocal*Vlocal*n_blades*chord*rad*(Cd*np.cos(phi)+Cl*np.sin(phi))
        #momentum check on inflow and swirl factors
        tem1=DtDr/(4.0*np.pi*rad*rho*V*V*(1+axi_induction))
        tem2=DqDr/(4.0*np.pi*rad*rad*rad*rho*V*(1+axi_induction)*omega)

        Prandtl_correction = PrandtlTipRootCorrection(rad/R, 0.1, 1, omega*R/V0, n_blades, tem1)
        if (Prandtl_correction < 1e-3): 
           Prandtl_correction = 1e-3
          
        #stabilise iteration
        anew=0.5*(axi_induction+tem1*Prandtl_correction)
        bnew=0.5*(swirl+tem2*Prandtl_correction)

        #check for convergence
        if (abs(anew-axi_induction)<1.0e-5):
            if (abs(bnew-swirl)<1.0e-5):
                finished=True
                
        axi_induction=anew
        swirl=bnew
        #increment iteration count
        total_sum=total_sum+1
        #check to see if iteration stuck
        if (total_sum>500):
            finished=True
            itercheck=1
            print('iteration failed')
        thrust=DtDr*rstep
        torque=DqDr*rstep
    return [thrust, torque, itercheck]    
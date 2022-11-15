# -*- coding: utf-8 -*-
"""
Created on Wed Aug 11 14:59:22 2021

@author: thoma
"""

import pandas as pd
import numpy as np
from scipy.interpolate import griddata, interp1d

def get_polar_fit():
    Re_1e4 = pd.read_csv('polar_files/4415_Re_1e4.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_5e4 = pd.read_csv('polar_files/4415_Re_5e4.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_15e5 = pd.read_csv('polar_files/4415_Re_15e5.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_3e5 = pd.read_csv('polar_files/4415_Re_3e5.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_5e5 = pd.read_csv('polar_files/4415_Re_5e5.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_1e6 = pd.read_csv('polar_files/4415_Re_1e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_13e6 = pd.read_csv('polar_files/4415_Re_13e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_15e6 = pd.read_csv('polar_files/4415_Re_15e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_20e6 = pd.read_csv('polar_files/4415_Re_20e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_23e6 = pd.read_csv('polar_files/4415_Re_23e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_25e6 = pd.read_csv('polar_files/4415_Re_25e6.txt', skipinitialspace=True, delim_whitespace=True) 
    Re_27e6 = pd.read_csv('polar_files/4415_Re_27e6.txt', skipinitialspace=True, delim_whitespace=True)
    Re_29e6 = pd.read_csv('polar_files/4415_Re_29e6.txt', skipinitialspace=True, delim_whitespace=True)  
    Re_32e6 = pd.read_csv('polar_files/4415_Re_32e6.txt', skipinitialspace=True, delim_whitespace=True)
    Re_35e6 = pd.read_csv('polar_files/4415_Re_35e6.txt', skipinitialspace=True, delim_whitespace=True)  
    Re_4e6  = pd.read_csv('polar_files/4415_Re_4e6.txt', skipinitialspace=True, delim_whitespace=True)
    Re_42e6  = pd.read_csv('polar_files/4415_Re_42e6.txt', skipinitialspace=True, delim_whitespace=True)
    
    liftfit = interp1d([1e4, 1.5e5, 3e5, 5e5, 1e6, 1.3e6, 1.5e6, 2e6, 2.3e6, 2.5e6, 2.7e6, 2.9e6, 3.2e6, 3.5e6, 4e6, 42e6], \
                      np.vstack([Re_1e4['CL'], Re_15e5['CL'], Re_3e5['CL'], Re_5e5['CL'], Re_1e6['CL'], \
                                 Re_13e6['CL'], Re_15e6['CL'], Re_20e6['CL'], Re_23e6['CL'], Re_25e6['CL'], \
                                 Re_27e6['CL'], Re_29e6['CL'], Re_32e6['CL'], Re_35e6['CL'], Re_4e6['CL'], \
                                 Re_42e6['CL']]), axis=0)
    dragfit = interp1d([1e4, 1.5e5, 3e5, 5e5, 1e6, 1.3e6, 1.5e6, 2e6, 2.3e6, 2.5e6, 2.7e6, 2.9e6, 3.2e6, 3.5e6, 4e6, 42e6], \
                      np.vstack([Re_1e4['CD'], Re_15e5['CD'], Re_3e5['CD'], Re_5e5['CD'], Re_1e6['CD'], \
                                 Re_13e6['CD'], Re_15e6['CD'], Re_20e6['CD'], Re_23e6['CD'], Re_25e6['CD'], \
                                 Re_27e6['CD'], Re_29e6['CD'], Re_32e6['CD'], Re_35e6['CD'], Re_4e6['CD'], \
                                 Re_42e6['CD']]), axis=0)
    alphafit = interp1d([1e4, 1.5e5, 3e5, 5e5, 1e6, 1.3e6, 1.5e6, 2e6, 2.3e6, 2.5e6, 2.7e6, 2.9e6, 3.2e6, 3.5e6, 4e6, 42e6], \
                      np.vstack([Re_1e4['alpha'], Re_15e5['alpha'], Re_3e5['alpha'], Re_5e5['alpha'], Re_1e6['alpha'], \
                                 Re_13e6['alpha'], Re_15e6['alpha'], Re_20e6['alpha'], Re_23e6['alpha'], Re_25e6['alpha'], \
                                 Re_27e6['alpha'], Re_29e6['alpha'], Re_32e6['alpha'], Re_35e6['alpha'], Re_4e6['alpha'], \
                                 Re_42e6['alpha']]), axis=0)
    return liftfit, dragfit, alphafit
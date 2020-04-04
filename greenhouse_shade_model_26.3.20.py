#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:16:43 2020

@author: bekah
"""


import pandas as pd
import numpy as np
import math as m
from scipy import optimize, stats
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

### location for OPV greenhouse at UA CEAC in Tucson, AZ ###
latitude = 32.28 # OPV greenhouse latitude (deg)
longitude = -110.94 # OPV greenhouse longitude (deg)
timezone = -7 # Tucson, AZ timezone (UTC)
elevation = 718 # OPV greenhouse elevation (m)

### input greenhouse dimensions ###
gh_width = 30.0 # in feet
gh_length = 48.0 # in feet

### input time interval information for shade simulation ### 
year = 2020 # year for simulation
m_1, m_2 = 3, 4 # months for simulation: 1 = January, 12 = December
d_1,d_2 = 1, 2 # days for simulation
hr_1, hr_2 = 9, 16 # hours for simulation: 0 = midnight, 6 = 6AM, 12 = 12PM, 
                   # 18 = 6PM
min_1, min_2 = 30, 31 # minute for simulation: 0 to 60 = entire hour

    
def gh_roof_profile():
    """Get west greenhouse roof profile measurements from excel file"""
    
    # roof measurements from excel file
    path_inputs = '/Users/bekah/Desktop/class/be523/modeling_project/roof_function.xlsx'
    sheet_name_setup = 'roof'
    roof = pd.read_excel(path_inputs, sheetname = sheet_name_setup)

    # make arrays for east and west roof sections
    x_west_measured = roof['x_west'] # x axis = greenhouse width, west side
    z_west_measured = roof['z_west'] # z axis = greenhouse roof heights, west side
    
    return x_west_measured, z_west_measured

def test_func(x, a, b, c, d):
    """Create third order polynomial test function for curve fit function"""
    
    return a + b * x + c * x**2 + d * x**3

def curve_fit():
    """Fit curve to data with third-order polynomial function from test_func()"""
    
    x_west_measured, z_west_measured = gh_roof_profile() # get west roof profile data

    param_west, param_cov_west = optimize.curve_fit(test_func, x_west_measured, z_west_measured) # get parameters for fitted curve

    # fitted roof geometry function: z = 7.29944696 + (1.27415518*x) + (-0.0680139854*x**2) + (0.00152035861*x**3)
    
    return param_west

def mirror(array):
    """Flip values of one array to create mirrored array"""
    
    mirror_array = np.flip(array,axis=0)
    
    return mirror_array


def create_model_gh_roof(gh_width):
    """Create model greenhouse roof profile from fitted curve"""
    
#    global gh_width, N_x # roof profile dependent only on width (symmetrical down length)

    a, b, c, d = curve_fit() # coefficients from curve fit function

    # calculate roof slopes
    N_x = 100 # number of steps for roof curve calculation
    gh_width_half = gh_width/2.0 # one roof pitch side
    dx = gh_width_half/100.0
   
    x_W = np.linspace(0,(N_x)*dx,N_x+1) # array for width
    z_W = np.zeros(N_x+1) # array for height
    z_slopes_W = np.zeros(N_x+1) # array for slopes of roof 
    roof_lengths_W = np.zeros(N_x+1)
    total_length_W = 0

    # calculate roof heights and slopes (zvalues from test function)
    for i in range(0,(N_x+1)): # starting with west end of greenhouse width
        z_W[i] = a + (b * x_W[i]) + (c * x_W[i]**2) + (d * x_W[i]**3) # function and parameters from curve fit 
        dz = z_W[i] - z_W[i-1] # dividing roof into segments of length dz
        z_slopes_W[i] = dz/dx # slope is change in x over change in z
        roof_lengths_W[i] = (dz**2 + dx**2)**0.5 # length of segment calculated from Pythagorean thrm
        if i == 0:
            z_slopes_W[0] = 0 # initial value for calculating slope of roof segments
            roof_lengths_W[0] = 0 # initial value for calculating length of roof segments 
        total_length_W += roof_lengths_W[i] # calculate total length of west roof pitch
            
    z_E, z_slopes_E = mirror(z_W), mirror(z_slopes_W) # create arrays for eastside roof heights and slopes
    z_E_drop_first = z_E[1:] # drop first item in eastside array to avoid doublecounting apex in entire roof array
    z_roof = np.concatenate((z_W, z_E_drop_first), axis=0) # create array for entire roof heights

    x_width = np.linspace(0,gh_width,(N_x * 2)+1) # array for width of greenhouse

    # roof lengths for east side
    roof_length_E = mirror(roof_lengths_W)
    total_length = (2 * total_length_W) # total length of greenhouse roof
    
    plt.figure(1)
    plt.plot(x_width, z_roof, 'b-')
    plt.axis('equal')
    plt.show()
    
    return z_roof, x_width, N_x


def solar_model(latitude, longitude, timezone, elevation, year, i_month, j_day, k_hr, l_min):
    """Calculate solar position based on time specified"""
    
    date=dt.datetime(year,i_month,j_day,k_hr,l_min)                
    
    time = date.timetuple().tm_yday
    time = time + k_hr/24 + l_min/(24*60)
    year = 2020
                    
    lat_r = latitude/180*np.pi
    lon_r = longitude/180*np.pi  
    n = 0
    for i in range(1900,year):
        if i%4 == 0:
            n += 366
        else:
            n+=365
    JulD = n + time + 2415018.5 - (timezone)/24
    LT = time - int(time)
    JC = (JulD - 2451545) / 36525
    x = 46.815 + JC * (0.00059 - JC * 0.001813)
    M_OE = 23 + (26 + (21.448 - JC * x) / 60) / 60
    EEO = 0.016708634 - JC * (0.000042037 + 0.0000001267 * JC)
    GMAS = 357.52911 + JC * (35999.05029 - 0.0001537 * JC)
    GMAS_r = m.radians(GMAS)
    GMLS = (280.46646 + JC * (36000.76983 + JC * 0.0003032))%360
    GMLS_r = m.radians(GMLS)
    Obliq_C = M_OE + 0.00256 * np.cos((125.04 - 1934.136 * JC) / 180 * np.pi)
    Obliq_C_r = m.radians(Obliq_C)
    SEC = np.sin(GMAS_r) * (1.914602 - JC * (0.004817 + 0.000014 * JC)) + np.sin(2 * GMAS_r) * (0.019993 - 0.000101 * JC) + np.sin(3 * GMAS_r) * 0.000289
    STL = GMLS + SEC
    SAL = STL - 0.00569 - 0.00478 * np.sin((125.04 - 1934.136 * JC) / 180 * np.pi)
    SAL_r = m.radians(SAL)
    sin_Delta = np.sin(Obliq_C_r) * np.sin(SAL_r)
    Delta_r = np.arcsin(sin_Delta)     #in radians   
    Var_y = np.tan((Obliq_C / 2) / 180 * np.pi) * np.tan((Obliq_C / 2) / 180 * np.pi)
    EOT_prime = Var_y * np.sin(2 * GMLS_r) - 2 * EEO * np.sin(GMAS_r) + 4 * EEO * Var_y * np.sin(GMAS_r) * np.cos(2 * GMLS_r) - 0.5 * Var_y * Var_y * np.sin(4 * GMLS_r) - 1.25 * EEO * EEO * np.sin(2 * GMAS_r)
    EOT = 4 * EOT_prime / np.pi * 180       
    TST = (LT * 1440 + EOT + 4 * longitude - 60 * timezone)%1440
    if TST / 4 < 0:
        Omega = TST/4+180
    else:
        Omega = TST/4 - 180   
    Omega_r = m.radians(Omega)
     
    cos_Zenith = np.sin(lat_r) * np.sin(Delta_r) + np.cos(lat_r) * np.cos(Delta_r) * np.cos(Omega_r)
    Zenith_r = np.arccos(cos_Zenith)             #in radians
    Aprime_r = np.arccos((np.sin(lat_r) * np.cos(Zenith_r) - np.sin(Delta_r)) / (np.cos(lat_r) * np.sin(Zenith_r)))
    Aprime = Aprime_r / np.pi * 180
    if Omega > 0:
        Azimuth = (Aprime + 180) % 360   #in degrees
    else:
        Azimuth = (540 - Aprime) % 360   #in degrees                
    Azimuth_r = Azimuth / 180 * np.pi
    Elev_angle = (np.pi)/2 - Zenith_r

    
    Azimuth_N = Azimuth_r   # Azimuth with north is the reference going east
    Elev_angle = np.pi/2 - Zenith_r
    Azimuth_S = np.pi + Azimuth_N      # Azimuth with the south is the reference and going west

    return Azimuth_N, Azimuth_S, Elev_angle
     
def shade(latitude, longitude, timezone, elevation, year, m_1, m_2, d_1, d_2, hr_1, hr_2, min_1, min_2):
    """Calculate shade offset using solar position from solar_model function and greenhouse roof values"""

    z_roof, x_width, N_x = create_model_gh_roof(gh_width)

    delta_y = np.zeros((N_x * 2) + 1) # change in shading on length 
    delta_x = np.zeros((N_x * 2) + 1) # change in shading on width
    floor_position_y_S = np.zeros((N_x * 2) + 1) # shade position on length, y axis
    floor_position_x = np.zeros((N_x * 2) + 1) # shade position on width, x axis
    floor_position_y_N = np.zeros((N_x * 2) + 1) # shade position on length, north edge  
    
    for i_month in range(m_1, m_2):
        for j_day in range(d_1, d_2):
            for k_hr in range(hr_1, hr_2):
                for l_min in range(min_1, min_2):
        
                    Azimuth_N, Azimuth_S, Elev_angle = solar_model(latitude, longitude, timezone, elevation, year, i_month, j_day, k_hr, l_min)
            
                    # calculate shade offset
                    for i in range((N_x * 2) + 1):
                        delta_y[i] = (z_roof[i] / np.tan(Elev_angle)) * np.cos(Azimuth_S) # shading offset for length
                        delta_x[i] = (z_roof[i] / np.tan(Elev_angle)) * np.sin(Azimuth_S) # shading offset for width
                        floor_position_x[i] = x_width[i] + delta_x[i] # position on floor from x point + shading offset
                        floor_position_y_S[i] = 24.67 + delta_y[i] # south edge of OPV panels on length of roof
                        floor_position_y_N[i] = 38 + delta_y[i] # north edge of OPV panels on length of roof
                        
                        
                         
                    plt.figure(1)
                    plt.plot(floor_position_x, floor_position_y_S, floor_position_x, floor_position_y_N)
                    
shade(latitude, longitude, timezone, elevation, year, m_1, m_2, d_1, d_2, hr_1, hr_2, min_1, min_2)

def oasis_data():
    """Import UA OASIS station data for radiation, air temperature, and relative humidity"""
    
    df_oasis = pd.read_csv('/Users/bekah/Desktop/opv_data/ivcurve_measurements/env_data/20191015.csv')
    df_oasis['time'] = pd.to_datetime(df_oasis['DATE (MM/DD/YYYY)'] + ' ' + df_oasis['MST'])
    df_oasis.set_index(keys='time', inplace=True)
    
    return df_oasis




### NEXT STEPS ###
                    
# make plot look nicer, heat map for roof pitch elevation
# make greenhouse footprint on same plot
# OPV panels as objects, plot on graph
# join south edge and north edge for total coverage
# transmittance coming through covering
# where are plants? percentage of plants shaded at any given time
# at every point, you're getting a shading value on floor + incident solar irradiance + transmitted solar irradiance

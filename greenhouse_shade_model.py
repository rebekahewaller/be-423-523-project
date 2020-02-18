#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 21:16:43 2020

@author: bekah
"""

# direct radiation - determine solar angle based on hour, day of year, 
# lattitude, longitude for greenhouse & determine amount of direct solar 
# radiation for clear sky day


# DIFFUSE RADIATON - determine amount of diffuse radiation based on hour,
# day of year, latitude, longitude for greenhouse for clear sky day 

# SHADING - determine shading due to OPV panels deployed on top of
# greenhouse roof

# GROWTH - determine growth (biomass & fruit) for lettuce and tomato plants
# based on daily light integral (DLI) under shade and non-shade treatment
# 


# direct radiation

import pandas as pd
import numpy as np
import math as m
from scipy import optimize, stats
import datetime as dt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def location_input():
    """
    Input geographical location, elevation, and timezone for solar model
    """
    
    # location for OPV greenhouse at UA CEAC in Tucson, AZ
    
    latitude = 32.28 # OPV greenhouse latitude (deg)
    longitude = -110.94 # OPV greenhouse longitude (deg)
    timezone = -7 # Tucson, AZ timezone (UTC)
    elevation = 718 # OPV greenhouse elevation (m)
    
    return latitude, longitude, timezone, elevation
    
def time_input():
    """
    Input time data for solar model
    """
    
    year = 2020
    month = 3 # number 
    day = 12 # number in month
    hour = 12 # integer between 9 (= 9:00AM) and 17 (= 4:00PM) ## CHECK THIS
    minute = 0 # float between 0 (= 0 min) to 0.983 = 59 min)
    
    date=dt.datetime(year,month,day)
    time = date.timetuple().tm_yday
    time = time + hour/24 + minute/24/60
    
    return year, time


def solar_model():
    """
    The function calculates for zenith angle and azimuth angle in radians for horizontal surface.
    time: (doy)
    starting_DOY:starting day for the simulation
    year: the year of the simulation
    location: city
    elevation (m):
    longitude (deg):
    latitude (deg):
    timezone: UTC time zone
    """
    
    latitude, longitude, timezone, elevation = location_input()
    year, time = time_input()

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

    
    # calculate incidence angle
    # Beta is equal to angle of tilted surface to horizontal (in radians)
    Beta = 45 # in degrees
    Beta_r = m.radians(Beta)
    
    cos_incidence = np.sin(Delta_r)* np.sin(lat_r) * np.cos(Beta_r) - np.sin(Delta_r) * np.cos(lat_r) * np.sin(Beta_r) * np.cos(Azimuth_r) + np.cos(Delta_r) * np.cos(lat_r) * np.cos(Beta_r) * np.cos(Omega_r) + np.cos(Delta_r) * np.sin(lat_r) * np.sin(Beta_r) * np.cos(Azimuth_r) * np.cos(Omega_r) + np.cos(Delta_r) * np.sin(Beta_r) * np.sin(Azimuth_r) * np.sin(Omega_r)                      
    incidence_ang_r = np.arccos(cos_incidence)
    
    return Delta_r, lat_r, Omega_r, Zenith_r, Azimuth_r, Elev_angle

def greenhouse_orientation():
    """ 
    Input orientation of OPV greenhouse
    """
    
    # NEED TO CHECK THIS WITH COMPASS (OR IPHONE)
    orientation_angle = 90 # angle between east-west line and the length of the greenhouse (0-90 degree)
    orientation_angle = float(orientation_angle)
    
def roof_measurements():
    """
    Read in measurements from scale drawing of OPV greenhouse roof
    """
    
    # roof measurements from excel
    path_inputs = '/Users/bekah/Desktop/class/be523/modeling_project/roof_function.xlsx'
    sheet_name_setup = 'roof'
    roof = pd.read_excel(path_inputs, sheetname = sheet_name_setup)
    
    # make arrays for east and west roof sections
    x_west = roof['x_west']
    z_west = roof['z_west']
    x_east = roof['x_east']
    z_east = roof['z_east']
    
    return x_west, z_west, x_east, z_east
    
def test_func(x, a, b, c, d):
    """
    Create third order polynomial test function for curve fit function
    """
    return a + b * x + c * x**2 + d * x**3

def curve_fit():
    """
    Fit roof measurements to curve with curve fit function
    """
    x_west, z_west, x_east, z_east = roof_measurements()
    
    # find best curve fit for west roof section
    param_west, param_cov_west = optimize.curve_fit(test_func, x_west, z_west)
    print(param_west)
    
    # z = 7.29944696 + (1.27415518*x) + (-0.0680139854*x**2) + (0.00152035861*x**3)
    
    z_west_fitted = test_func(x_west,*param_west)
    
    # mirror curve for east roof section
    z_east_fitted = np.flip(z_west_fitted,axis=0)
    
    # create array for both west and east roof sections
    x_whole_roof = np.concatenate((x_west, x_east), axis=0)
    z_whole_roof = np.concatenate((z_west_fitted, z_east_fitted), axis=0)

    # plot roof
    plt.plot(x_whole_roof, z_whole_roof, '-', color ='blue', label="roof curve") 
    plt.axis('equal') # make axes square
    plt.show() 

    return param_west

def section_coordinates():
    """
    Create array for x, y, z coordinates for 1 ft square sections of greenhouse footprint
    """
    
    gh_width = 30.0 # in feet
    gh_width_west = gh_width/2.0
    N_x = 100
    dx = gh_width_west/100.0
    gh_length = 48 # in feet
    
    xvalues = np.linspace(0,(N_x)*dx,N_x+1) # array for width
    yvalues = np.linspace(0,gh_length,num=gh_length+1) # array for height
    zvalues_west = np.zeros(N_x+1) # array for height
    
    for i in range(0,len(xvalues)):
        zvalues_west[i] = 7.29944696 + (1.27415518*xvalues[i]) + (-0.0680139854*xvalues[i]**2) + (0.00152035861*xvalues[i]**3)
        i += 1
 
    roof_slopes_west = np.zeros(N_x+1)
    roof_lengths = np.zeros(N_x+1)

    total_length_west = 0

    for i in range(1,len(xvalues)):
        dz = zvalues_west[i] - zvalues_west[i-1]
        roof_slopes_west[i] = dz/dx
        roof_lengths[i] = (dz**2 + dx**2)**0.5
        total_length_west += roof_lengths[i]
    
    zvalues_east = np.flip(zvalues_west, axis=0)
    zvalues_west = zvalues_west[:-1]
    zvalues = np.concatenate((zvalues_west, zvalues_east), axis=0)
    
    xx, yy = np.meshgrid(xvalues, yvalues)      
    
    plt.plot(xx, yy, marker='.', color='k', linestyle='none')
    plt.axis('equal')
    plt.show() 

    return roof_slopes_west


def calc_incidence_angle():
    """ Calculate incidence angle for every slope agle on greenhouse roof """
    
    Delta_r, lat_r, Omega_r, Zenith_r, Azimuth_r, Elev_angle = solar_model()
    
    # Beta is equal to angle of tilted surface to horizontal (in radians)
    roof_slopes_west = section_coordinates()
    Beta_r = np.arctan(roof_slopes_west) 
    incidence_angles_west = np.zeros(101)
    
    
    for i in range(0,len(roof_slopes_west)):
        incidence_angles_west[i] = np.arccos(np.sin(Delta_r)* np.sin(lat_r) * np.cos(Beta_r[i]) - np.sin(Delta_r) * np.cos(lat_r) * np.sin(Beta_r[i]) * np.cos(Azimuth_r) + np.cos(Delta_r) * np.cos(lat_r) * np.cos(Beta_r[i]) * np.cos(Omega_r) + np.cos(Delta_r) * np.sin(lat_r) * np.sin(Beta_r[i]) * np.cos(Azimuth_r) * np.cos(Omega_r) + np.cos(Delta_r) * np.sin(Beta_r[i]) * np.sin(Azimuth_r) * np.sin(Omega_r))
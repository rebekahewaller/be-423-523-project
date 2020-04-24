# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:02:31 2020

@author: Nicole Wigtil
"""

#Assume Size of Greenhouse = 24X36 ft 

import numpy as np
import matplotlib.pyplot as plt


    #Average Mature Fruit Size (g of dry weight)

Num_Plants = int(input("Number of Plants/m^2 = ")) #using 3 plants/m^2 in the greenhouse
Nodes_OPV = int(input("Number of Nodes for OPV section (Average per week) = ")) #using 3 nodes per week
Nodes = int(input("Number of Nodes for non-OPV section (Average per week) = ")) #using 3 nodes per week
Days = int(input("Number of Days of Growth = "))
Temp_OPV = int(input("Temperature for OPV section (C) = "))
Temp = int(input("Temperature for the non-OPV section is (C) = "))

GH_Width = int(input("Width of the greenhouse (in feet)= ")) #allows user to choose for any greenhouse size
GH_Length = int(input("Length of the greenhouse (in feet)= "))


GH = GH_Width * GH_Length * 0.092903 #converts from ft^2 to m^2
dt = 1 #timestep in days
N_t = int(Days/dt)
t = np.linspace(0,(N_t+1)*dt, N_t+1)
avg_weight = 250 #grams


#function for the number of mature fruit in the OPV section(N/m^2)
def Num_Mature_Fruit_OPV (num_node, rho): #inputs are the number of nodes per week, total time period, and plant density
    N = (num_node/7) #number of nodes per week divided by 7 to make number of nodes per day
    M = np.zeros(N_t+1) 
    for n in range(1, N_t + 1):
        M[n] = M[n-1]+ (rho * (N/4) * 3 * dt) #3 represents the number of mature fruit for each node produces, divide by 4 every fourth node is mature fruit growth
    return M;

Mature_Fruit_Per_SQM_OPV = Num_Mature_Fruit_OPV(Nodes_OPV, Num_Plants)
#print("Number of mature fruit per meter squared=", Mature_Fruit_Per_SQM)

plt.figure(1)
plt.plot(t, Mature_Fruit_Per_SQM_OPV, 'r')
plt.legend(['Mature Fruit in OPV (N/m^2)'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Number of Tomatoes per m^2')

#function for the number of mature fruit (N/m^2)
def Num_Mature_Fruit (num_node, rho): #inputs are the number of nodes per week, total time period, and plant density
    N = (num_node/7) #number of nodes per week divided by 7 to make number of nodes per day
    M = np.zeros(N_t+1) 
    for n in range(1, N_t + 1):
        M[n] = M[n-1] + (rho * (N/4) * 3 * dt) #3 represents the number of mature fruit for each node produces, divide by 4 every fourth node is mature fruit growth
    return M;

Mature_Fruit_Per_SQM = Num_Mature_Fruit(Nodes, Num_Plants)
#print("Number of mature fruit per meter squared=", Mature_Fruit_Per_SQM)

plt.figure(2)
plt.plot(t, Mature_Fruit_Per_SQM, 'r')
plt.legend(['Mature Fruit (N/m^2)'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Number of Tomatoes per m^2')

    
#function for mature fruit weight in the OPV section (g/m^2)   
def Fruit_Weight_Mat_OPV (T, rho):
    F_W = np.zeros(N_t+1)
    if(T >= 0 and T < 20): # the weight is dependent on the temperature input
        for n in range(1, N_t+1):
            F_W[n] = (avg_weight + 0.15 * F_W[n-1]) * rho
    elif(T >= 20 and T <= 25):
        for n in range(1, N_t+1):
            F_W[n] = (avg_weight) * rho
    elif(T > 25 and T < 38):
        for n in range (1, N_t+1):
            F_W[n] = (avg_weight - 0.15 * F_W[n-1]) * rho
    else:
        for n in range (1, N_t+1):
            F_W[n] = 0
    return(F_W);
    
Mature_Fruit_Weight_OPV = Fruit_Weight_Mat_OPV(Temp_OPV, Num_Plants)

plt.figure(3)
plt.plot(t, Mature_Fruit_Weight_OPV, 'b')
plt.legend(['Mature Fruit Weight in OPV (g/m^2)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight per m^2')


#function for mature fruit weight in non-OPV section (g/m^2)   
def Fruit_Weight_Mat (T, rho):
    F_W = np.zeros(N_t+1)
    if(T >= 0 and T < 20): # the weight is dependent on the temperature input
        for n in range(1, N_t+1):
            F_W[n] = (avg_weight + 0.15 * F_W[n-1]) * rho
    elif(T >= 20 and T <= 25):
        for n in range(1, N_t+1):
            F_W[n] = (avg_weight) * rho
    elif(T > 25 and T < 38):
        for n in range (1, N_t+1):
            F_W[n] = (avg_weight - 0.15 * F_W[n-1]) * rho
    else:
        for n in range (1, N_t+1):
            F_W[n] = 0
    return(F_W);
    
Mature_Fruit_Weight = Fruit_Weight_Mat(Temp, Num_Plants)

plt.figure(4)
plt.plot(t, Mature_Fruit_Weight, 'b')
plt.legend(['Mature Fruit Weight (g/m^2)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight per m^2')


#function for the average mature fruit size in grams of dry weight for the OPV        
def Avg_Dry_Weight_OPV (mat_weight, rho):
    D_W = (mat_weight / rho) * 0.7      #dry weight is estimated to be 70% of fresh weight
    return(D_W);
    
Avg_Fruit_Size_Dry_OPV = Avg_Dry_Weight_OPV(Mature_Fruit_Weight_OPV, Num_Plants)

plt.figure(5)
plt.plot(t, Avg_Fruit_Size_Dry_OPV, 'g')
plt.legend(['Average Mature Fruit Size OPV (g dry weight)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight (g)')

#function for the average mature fruit size in grams of dry weight with no OPV        
def Avg_Dry_Weight (mat_weight, rho):
    D_W = (mat_weight / rho) * 0.7      #dry weight is estimated to be 70% of fresh weight
    return(D_W);
    
Avg_Fruit_Size_Dry = Avg_Dry_Weight(Mature_Fruit_Weight, Num_Plants)
    
plt.figure(6)
plt.plot(t, Avg_Fruit_Size_Dry, 'g')
plt.legend(['Average Mature Fruit Size (g dry weight)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight (g)')    
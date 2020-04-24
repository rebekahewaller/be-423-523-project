
"""
Created on Sun Mar 29 11:52:10 2020

@author: Caroline Schulte
"""
# Assumptions: Max DLI plants can handle is 60 mol/m^2/day
# Best DLI for tomatoes is 30 mol/m^2/day
# DLI ideal is between 22 mol/m^2/day and 30
# Assume lowest DLI of 5 


import numpy as np
import matplotlib.pyplot as plt

Temp_OPV = int(input("Temperature for OPV section (C) = "))
Temp = int(input("Temperature for the non-OPV section is (C) = "))
Days = int(input("Number of Days of Growth = "))

DLI_OPV_Input= int(input("DLI for OPV section (mol/m^2/day) is = "))
DLI_Input = int(input("DLI for non-OPV section (mol/m^2/day) is = "))


#print("temperature input for OPV section in degrees C is ", Temp_OPV)
#print("temperature input in degrees C for the non OPV section is ", Temp)

dt = 1 # time step (days)
N_t = int(Days/dt)
t = np.linspace(0,(N_t+1)*dt, N_t+1)



# OPV - covered 
def DLI_OPV (DLI): # Function to calculate the effect of DLI on the growth of the plants given optimal DLI is 30 mol/m^2/day and minimum is 20 mol/m^2/day
    if(DLI >= 22 and DLI <= 30): 
        DLI_effect = 1.0 + 0.05625 * (DLI - 30) # (0.55-1) / (22-30) = x = 0.05625
    elif(DLI > 30 and DLI < 45): 
         DLI_effect = 1.0 - 0.066 * (DLI - 30) # 1 / (45 - 30) = x = 0.066
    elif(DLI < 22 and DLI > 5):
        DLI_effect = 1.0 + 0.04 * (DLI - 30) # - 1 / (5-30) = x = 0.04
    else:
        DLI_effect = 0
    #print(F_T)
    return(DLI_effect);
    
DLI_OPV_Effect = DLI_OPV(DLI_OPV_Input)


def F_N_OPV(T): # Function to modify node development rate based on temperature
# This equation assumes that the ideal temp is 28 C for tomato growth
    
    if(T > 20 and T <= 22): 
        F_T = 1.0 + 0.225 * (T - 22)  # (0.55 - 1) / (20-22) = x = 0.225
    elif(T > 22 and T < 35):
        F_T = 1.0 - 0.0455 * (T - 22) # 1 / (35-22) = x = 0.0769
    elif (T < 20 and T > 10 ):
        F_T = 1.0 + 0.0833 * (T-22) # -1 / (10 - 22) = x = 0.0833
    else:
        F_T = 0
    #print(F_T)
    return(F_T);
    
F_T_OPV = F_N_OPV(Temp_OPV) # The output from the function to modify node development rate 

#print("Function to modify node development for OPV is", F_T_OPV)


def dNdt_OPV(FN, DLI): # The function for change in node development versus time using the output from the first function
  Nm = (3/7) # Constant for maximum rate of node appearance (node/day) given 3 nodes develop per week
  dN_dt = Nm * FN * DLI
  return (dN_dt);

dN_dt_OPV = dNdt_OPV(F_T_OPV, DLI_OPV_Effect) # Calling node development function
#print ("Rate of node development for OPV = ", dN_dt_OPV)


def Nodes_developed_OPV(FN, DLI): # Nodes developed given F_T and how many weeks have passed
    N_m = (3/7) # constant. 3 nodes per 1 foot of growth
    N = np.zeros(N_t+1)
    for n in range(1, N_t + 1):
        N[n] = N[n-1]+ (N_m*FN*DLI)*dt
    return N

Nodes_OPV = Nodes_developed_OPV(F_T_OPV, DLI_OPV_Effect) #input the number of days
#print ("Number of nodes developed for OPV is = ", Nodes_OPV)


def Flowers_developed_OPV(Num_nodes):
    T = Num_nodes/3 # Number of trusses
    F = T * 3 # Number of flowers
    return(F, T)
    
Flowers_OPV, Trusses_OPV = Flowers_developed_OPV(Nodes_OPV)
#print("flowers and trusses developed for OPV are", Flowers_OPV, Trusses_OPV)


def Head_Growth_OPV(Num_nodes):
    Growth = (Num_nodes/3) # 3 nodes for every 1 foot of growth
    return(Growth)
    
Head_Growth_OPV_Total = Head_Growth_OPV(Nodes_OPV)
#print("The plant has grown this many feet under the OPV section= ", Head_Growth_OPV_Total)

# Graphs for OPV Section
plt.figure(1)
plt.plot(t,Head_Growth_OPV_Total,'r')
plt.legend(['OPV Growth of Head'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Growth (ft)')

plt.figure(2)
plt.plot(t,Trusses_OPV,'g', t, Flowers_OPV, 'y', t, Nodes_OPV, 'b')
plt.legend(['OPV Trusses', 'OPV Flowers', 'OPV Nodes'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Number')





# Non-Covered

dt = 1 # time step (days)
N_t = int(Days/dt)
t = np.linspace(0,(N_t+1)*dt, N_t+1)

def DLI (DLI): # Function to calculate the effect of DLI on the growth of the plants given optimal DLI is 30 mol/m^2/day and minimum is 20 mol/m^2/day
   if(DLI >= 22 and DLI <= 30): 
        DLI_effect = 1.0 + 0.05625 * (DLI - 30) # (0.55-1) / (22-30) = x = 0.05625
   elif(DLI > 30 and DLI < 45):
         DLI_effect = 1.0 - 0.066 * (DLI - 30) # 1 / (45 - 30) = x = 0.066
   elif(DLI < 22 and DLI > 5):
        DLI_effect = 1.0 + 0.04 * (DLI - 30) # - 1 / (5-30) = x = 0.04
   else:
        DLI_effect = 0
    #print(F_T)
   return(DLI_effect);
    
DLI_Effect = DLI(DLI_Input)

def F_N(T): # Function to modify node development rate based on temperature
# This equation assumes that the ideal temp is 28 C for tomato growth
    
    if(T > 20 and T <= 22): 
        F_T = 1.0 + 0.225 * (T - 22)  # (0.55 - 1) / (20-22) = x = 0.225
    elif(T > 22 and T < 35):
        F_T = 1.0 - 0.0455 * (T - 22) # 1 / (35-22) = x = 0.0769
    elif (T < 20 and T > 10 ):
        F_T = 1.0 + 0.0833 * (T-22) # -1 / (10 - 22) = x = 0.0833
    else:
        F_T = 0
    #print(F_T)
    return(F_T);
    
F_T = F_N(Temp) # The output from the function to modify node development rate 

#print("Function to modify node development for non-OPV is", F_T)


def dNdt(FN, DLI): # The function for change in node development versus time using the output from the first function
  Nm = (4.5/7) # Constant for maximum rate of node appearance (node/day) given 4.5 nodes develop per week because 1.5 ft growth per week
  dN_dt = Nm * FN * DLI
  return (dN_dt);

dN_dt = dNdt(F_T, DLI_Effect) # Calling node development function using the output from the first function
#print ("Rate of node development for non-OPV = ", dN_dt)


def Nodes_developed(FN,DLI): # Nodes developed given F_T and how many weeks have passed
    N_m = (4.5/7) # constant. 3 nodes per 1 foot of growth
    N = np.zeros(N_t+1)
    for n in range(1, N_t + 1):
        N[n] = N[n-1]+ (N_m*FN*DLI)*dt
    return N

Nodes = Nodes_developed(F_T, DLI_Effect) #input the number of days
#print ("Number of nodes developed for non-OPV is = ", Nodes)


def Flowers_developed(Num_nodes):
    T = Num_nodes/3 # Number of trusses
    F = T * 3 # Number of flowers
    return(F,T) 
    
Flowers, Trusses = Flowers_developed(Nodes)
#print("flowers and trusses developed for non-OPV are", Flowers, Trusses)


def Head_Growth(Num_nodes):
    Growth = (Num_nodes/3) # 3 nodes for every 1 foot of growth
    return(Growth)
    
Head_Growth_Total = Head_Growth(Nodes)
#print("The plant has grown this many feet = ", Head_Growth_Total)

# Graphs for Non-OPV      
plt.figure(3)
plt.plot(t,Head_Growth_Total,'r')
plt.legend(['Growth of Head'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Growth (ft)')

plt.figure(4)
plt.plot(t,Trusses,'g', t, Flowers, 'y', t, Nodes, 'b')
plt.legend(['Trusses', 'Flowers', 'Nodes'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Number')
        




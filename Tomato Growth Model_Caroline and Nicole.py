# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 21:02:31 2020

@author: Nicole Wigtil

@author: Caroline Schulte
"""

# Assumptions: Max DLI plants can handle is 60 mol/m^2/day
# Best DLI for tomatoes is 30 mol/m^2/day
# DLI ideal is between 22 mol/m^2/day and 30
# Assume lowest DLI of 5 


import pandas as pd # New Library to make dataframes for data input, output and printing multiple variables
import numpy as np
import matplotlib.pyplot as plt
#Assume Size of Greenhouse = 24X36 ft 


# This is for OPV
excel_file = 'opv_greenhouse_env_data.xlsx'
cols = [0, 1, 11] # Weeks are 0, Temp is 1, DLI is 2
Inputs_OPV = pd.read_excel(excel_file, sheet_name = 0, index_col = None, usecols=cols)

print("Data is", Inputs_OPV)
print("There are this many weeks: ", len(Inputs_OPV))

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
    

#print("Function to modify node development for OPV is", F_T_OPV)


def Nodes_Developed_OPV(FN, DLI): # The function for change in node development versus time using the output from the first function
  Nm = (3) # Constant for maximum rate of node appearance (node/week) given 3 nodes develop per week
  N = Nm * FN * DLI # This is over one week so multiply by 7 for 7 days
  return (N);

#print ("Rate of node development for OPV = ", dN_dt_OPV)


def Flowers_developed_OPV(Num_nodes):
    T = Num_nodes/3 # Number of trusses - one truss per every three nodes
    F = T * 3 # Number of flowers - 3 flowers for every truss (should be same as number of nodes)
    return(F, T)
    
#print("flowers and trusses developed for OPV are", Flowers_OPV, Trusses_OPV)


def Head_Growth_OPV(Num_nodes):
    Growth = (Num_nodes/3) # 3 nodes for every 1 foot of growth
    return(Growth)
    
#function for the number of mature fruit in the OPV section(N/m^2)
def Num_Mature_Fruit_OPV (num_node, rho): #inputs are the number of nodes per week, total time period, and plant density
    N = (num_node/7) #number of nodes per week divided by 7 to make number of nodes per day
    M = np.zeros(N_t+1) 
    for n in range(1, N_t + 1):
        M[n] = M[n-1]+ (rho * (N/4) * 3 * dt) #3 represents the number of mature fruit for each node produces, divide by 4 every fourth node is mature fruit growth
    return M;

Mature_Fruit_Per_SQM_OPV = Num_Mature_Fruit_OPV(Nodes_OPV, Num_Plants)
#print("Number of mature fruit per meter squared=", Mature_Fruit_Per_SQM)

plt.figure(8)
plt.plot(t, Mature_Fruit_Per_SQM_OPV, 'r')
plt.legend(['Mature Fruit in OPV (N/m^2)'], loc = 'upper left')
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

plt.figure(4)
plt.plot(t, Mature_Fruit_Weight_OPV, 'b')
plt.legend(['Mature Fruit Weight in OPV (g/m^2)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight per m^2')

#function for the average mature fruit size in grams of dry weight for the OPV        
def Avg_Dry_Weight_OPV (mat_weight, rho):
    D_W = (mat_weight / rho) * 0.7      #dry weight is estimated to be 70% of fresh weight
    return(D_W);
    
Avg_Fruit_Size_Dry_OPV = Avg_Dry_Weight_OPV(Mature_Fruit_Weight_OPV, Num_Plants)

plt.figure(6)
plt.plot(t, Avg_Fruit_Size_Dry_OPV, 'g')
plt.legend(['Average Mature Fruit Size OPV (g dry weight)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight (g)')


i = 0
DLI_OPV_Effect = np.zeros(len(Inputs_OPV))
F_T_OPV = np.zeros(len(Inputs_OPV))
OPV_Nodes_Developed = np.zeros(len(Inputs_OPV))
Flowers_OPV = np.zeros(len(Inputs_OPV))
Trusses_OPV = np.zeros(len(Inputs_OPV))
OPV_Head_Growth = np.zeros(len(Inputs_OPV))
Week = np.zeros(len(Inputs_OPV))

Total_OPV_Flowers = np.zeros(len(Inputs_OPV))
Total_OPV_Trusses = np.zeros(len(Inputs_OPV))
Total_OPV_Nodes= np.zeros(len(Inputs_OPV))
Total_OPV_Head_Growth = np.zeros(len(Inputs_OPV))

while (i < len(Inputs_OPV)):
    # Import the values from the excel document
    Temp_OPV = Inputs_OPV.iloc[i,1] # Temperature values are in column 1 
    DLI_OPV_Input = Inputs_OPV.iloc[i,2] # DLI values are in column 11

    
    DLI_OPV_Effect[i] = DLI_OPV(DLI_OPV_Input) # Calculates the value for effect of DLI on the growth of the plant



    F_T_OPV[i] = F_N_OPV(Temp_OPV) # The output from the function to modify node development rate 
    
  
    OPV_Nodes_Developed[i] = Nodes_Developed_OPV(F_T_OPV[i], DLI_OPV_Effect[i]) # How many nodes developed that week based on on DLI and T inputs
    #print("Nodes", OPV_Nodes_Developed[i])
    Total_OPV_Nodes[i] = OPV_Nodes_Developed[i] + Total_OPV_Nodes[i-1] # How many nodes have developed over all the weeks
    #print("Total Nodes", Total_Nodes_Developed[i])
    
    Flowers_OPV[i], Trusses_OPV[i] = Flowers_developed_OPV(OPV_Nodes_Developed[i]) # How many flowers developed that week
    Total_OPV_Flowers[i] = Flowers_OPV[i] + Total_OPV_Flowers[i-1]
    Total_OPV_Trusses[i] = Trusses_OPV[i] + Total_OPV_Trusses[i-1]
    
    
    OPV_Head_Growth[i] = Head_Growth_OPV(OPV_Nodes_Developed[i]) # How much the head grew that week
    Total_OPV_Head_Growth[i] = OPV_Head_Growth[i] + Total_OPV_Head_Growth[i-1]
    
    Week[i] = i
    
    i = i + 1
    
    
print("The plant has grown this many feet under the OPV section= ", Total_OPV_Head_Growth[len(Inputs_OPV)-1])
print("The plant has grown this many trusses under the OPV section= ", Total_OPV_Trusses[len(Inputs_OPV)-1])
print("The plant has grown this many flowers under the OPV section= ", Total_OPV_Flowers[len(Inputs_OPV)-1])
print("The plant has grown this many nodes under the OPV section= ", Total_OPV_Nodes[len(Inputs_OPV)-1])




plt.figure(1)
plt.plot(Week, OPV_Head_Growth, label = 'Weekly OPV Head Growth (ft)')
plt.plot(Week, Trusses_OPV, label = 'Weekly OPV Trusses')
plt.plot(Week, Flowers_OPV, label = 'Weekly OPV Flowers')
plt.plot(Week, OPV_Nodes_Developed, label = 'Weekly OPV Nodes')
plt.xlabel('Weeks of Growth'); plt.ylabel('Growth Amount')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.show()

plt.figure(2)
plt.plot(Week, Total_OPV_Head_Growth, label = 'Total OPV Head Growth (ft)')
plt.plot(Week, Total_OPV_Trusses, label = 'Total OPV Trusses')
plt.plot(Week, Total_OPV_Flowers, label = 'Total OPV Flowers')
plt.plot(Week, Total_OPV_Nodes, label = 'Total OPV Nodes')
plt.xlabel('Weeks of Growth'); plt.ylabel('Growth Amount')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.show()









# This is for Non-OPV
cols = [0, 2, 12] # Weeks are 0, Temp is 1, DLI is 2
Inputs_NonOPV = pd.read_excel(excel_file, sheet_name = 0, index_col = None, usecols=cols)

print("Data is", Inputs_NonOPV)
print("There are this many weeks: ", len(Inputs_NonOPV))

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


def Nodes_Dev(FN, DLI): # The function for change in node development versus time using the output from the first function
  Nm = (4.5) # Constant for maximum rate of node appearance (node) given 4.5 nodes develop per week because 1.5 ft growth per week
  N = Nm * FN * DLI
  return (N);


def Flowers_developed(Num_nodes):
    T = Num_nodes/3 # Number of trusses
    F = T * 3 # Number of flowers
    return(F,T) 
    


def Head_Growth(Num_nodes):
    Growth = (Num_nodes/3) # 3 nodes for every 1 foot of growth
    return(Growth)
    
#function for the number of mature fruit (N/m^2)
def Num_Mature_Fruit (num_node, rho): #inputs are the number of nodes per week, total time period, and plant density
    N = (num_node/7) #number of nodes per week divided by 7 to make number of nodes per day
    M = np.zeros(N_t+1) 
    for n in range(1, N_t + 1):
        M[n] = M[n-1] + (rho * (N/4) * 3 * dt) #3 represents the number of mature fruit for each node produces, divide by 4 every fourth node is mature fruit growth
    return M;

Mature_Fruit_Per_SQM = Num_Mature_Fruit(Nodes, Num_Plants)
#print("Number of mature fruit per meter squared=", Mature_Fruit_Per_SQM)

plt.figure(3)
plt.plot(t, Mature_Fruit_Per_SQM, 'r')
plt.legend(['Mature Fruit (N/m^2)'], loc = 'upper left')
plt.xlabel('t (days)'); plt.ylabel('Number of Tomatoes per m^2')

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

plt.figure(5)
plt.plot(t, Mature_Fruit_Weight, 'b')
plt.legend(['Mature Fruit Weight (g/m^2)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight per m^2')

#function for the average mature fruit size in grams of dry weight with no OPV        
def Avg_Dry_Weight (mat_weight, rho):
    D_W = (mat_weight / rho) * 0.7      #dry weight is estimated to be 70% of fresh weight
    return(D_W);
    
Avg_Fruit_Size_Dry = Avg_Dry_Weight(Mature_Fruit_Weight, Num_Plants)
    
plt.figure(7)
plt.plot(t, Avg_Fruit_Size_Dry, 'g')
plt.legend(['Average Mature Fruit Size (g dry weight)'], loc = 'lower right')
plt.xlabel('t (days)'); plt.ylabel('Average Weight (g)') 
    

i = 0
DLI_Effect = np.zeros(len(Inputs_NonOPV))
F_T = np.zeros(len(Inputs_NonOPV))
Nodes_Developed = np.zeros(len(Inputs_NonOPV))
Flowers = np.zeros(len(Inputs_NonOPV))
Trusses = np.zeros(len(Inputs_NonOPV))
Main_Growth = np.zeros(len(Inputs_NonOPV))
Week = np.zeros(len(Inputs_NonOPV))

Total_Flowers = np.zeros(len(Inputs_NonOPV))
Total_Trusses = np.zeros(len(Inputs_NonOPV))
Total_Nodes= np.zeros(len(Inputs_NonOPV))
Total_Head_Growth = np.zeros(len(Inputs_NonOPV))

while (i < len(Inputs_NonOPV)):
    # Import the values from the excel document
    Temp = Inputs_NonOPV.iloc[i,1] # Temperature values are in column 2
    DLI_Input = Inputs_NonOPV.iloc[i,2] # DLI values are in column 12

    
    DLI_Effect[i] = DLI(DLI_Input) # Calculates the value for effect of DLI on the growth of the plant


    F_T[i] = F_N(Temp) # The output from the function to modify node development rate 
    
  
    Nodes_Developed[i] = Nodes_Dev(F_T[i], DLI_Effect[i]) # How many nodes developed that week based on on DLI and T inputs
    #print("Nodes", OPV_Nodes_Developed[i])
    
    Total_Nodes[i] = Nodes_Developed[i] + Total_Nodes[i-1] # How many nodes have developed over all the weeks
    #print("Total Nodes", Total_Nodes_Developed[i])
    
    Flowers[i], Trusses[i] = Flowers_developed(Nodes_Developed[i]) # How many flowers developed that week
    Total_Flowers[i] = Flowers[i] + Total_Flowers[i-1]
    Total_Trusses[i] = Trusses[i] + Total_Trusses[i-1]
    
    
    Main_Growth[i] = Head_Growth(Nodes_Developed[i]) # How much the head grew that week
    Total_Head_Growth[i] = Main_Growth[i] + Total_Head_Growth[i-1]
    
    Week[i] = i
    
    i = i + 1
    
    
print("The plant has grown this many feet under the NON-OPV section= ", Total_Head_Growth[len(Inputs_NonOPV)-1])
print("The plant has grown this many trusses under the NON-OPV section= ", Total_Trusses[len(Inputs_NonOPV)-1])
print("The plant has grown this many flowers under the NON-OPV section= ", Total_Flowers[len(Inputs_NonOPV)-1])
print("The plant has grown this many nodes under the NON-OPV section= ", Total_Nodes[len(Inputs_NonOPV)-1])




plt.figure(1)
plt.plot(Week, Main_Growth, label = 'Weekly Head Growth (ft)')
plt.plot(Week, Trusses, label = 'Weekly Trusses')
plt.plot(Week, Flowers, label = 'Weekly Flowers')
plt.plot(Week, Nodes_Developed, label = 'Weekly Nodes')
plt.xlabel('Weeks of Growth'); plt.ylabel('Growth Amount')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.show()

plt.figure(2)
plt.plot(Week, Total_Head_Growth, label = 'Total Head Growth (ft)')
plt.plot(Week, Total_Trusses, label = 'Total Trusses')
plt.plot(Week, Total_Flowers, label = 'Total Flowers')
plt.plot(Week, Total_Nodes, label = 'Total Nodes')
plt.xlabel('Weeks of Growth'); plt.ylabel('Growth Amount')

plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.show()

plt.figure(2)
plt.plot(Week, Total_Head_Growth, label = 'Total Head Growth (ft)')
plt.plot(Week, Total_OPV_Head_Growth, label = 'Total OPV Head Growth (ft)')
plt.xlabel('Weeks of Growth'); plt.ylabel('Growth (ft)')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=2, mode="expand", borderaxespad=0.) 
plt.show()
 
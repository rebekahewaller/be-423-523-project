# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 12:00:11 2020

@author: Caroline Schulte
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 09:20:16 2020

@author: Caroline Schulte
"""
# Assumptions: Max DLI plants can handle is 60 mol/m^2/day
# Best DLI for tomatoes is 30 mol/m^2/day
# DLI ideal is between 22 mol/m^2/day and 30
# Assume lowest DLI of 5 


import pandas as pd # New Library to make dataframes for data input, output and printing multiple variables
import numpy as np
import matplotlib.pyplot as plt


# This is for OPV
excel_file = 'opv_greenhouse_env_data.xlsx'
cols = [0, 1, 11] # Weeks are 0, Temp is 1, DLI is 2
Inputs_OPV = pd.read_excel(excel_file, sheet_name = 0, index_col = None, usecols=cols)

print("Data is", Inputs_OPV)
print("There are this many weeks: ", len(Inputs_OPV))

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
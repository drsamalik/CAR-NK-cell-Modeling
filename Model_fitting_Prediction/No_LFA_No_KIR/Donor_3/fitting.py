from scipy.optimize import least_squares
import numpy as np
import time
def Residue_Fit(x0,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data):
    y0_Gen4 = [elem for i, elem in enumerate(x0) if i != 1]
    yM_NK_Gen4 = Sys_CAR_Gen4.Effector_vs_Lysis(y0_Gen4)
    yD_NK_Gen4 = (Gen4_NK_data - yM_NK_Gen4)
    y0_Gen2 = [elem for i, elem in enumerate(x0) if i != 0]
    yM_NK_Gen2 = Sys_CAR_Gen2.Effector_vs_Lysis(y0_Gen2)
    yD_NK_Gen2 = (Gen2_NK_data - yM_NK_Gen2)
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0_Gen2)
    yD_WT = (WT_NK_data - yM_WT)
    print(x0.tolist())
    cost = sum(yD_NK_Gen4**2 + yD_NK_Gen2**2 + yD_WT**2)
    print('cost', cost)
    return 2.0*np.concatenate((yD_NK_Gen4, yD_NK_Gen2, yD_WT))

def Model_fit(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data,fit=True):
    if fit:
        result = least_squares(Residue_Fit, x0, args=(Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data), bounds=(LB,UB))
        x0 = (result.x).tolist()
    y0 = x0
    y0_Gen4 = [elem for i, elem in enumerate(y0) if i != 1]
    yM_NK_Gen4 = Sys_CAR_Gen4.Effector_vs_Lysis(y0_Gen4)
    yD_NK_Gen4 = (Gen4_NK_data - yM_NK_Gen4)
    y0_Gen2 = [elem for i, elem in enumerate(y0) if i != 0]
    yM_NK_Gen2 = Sys_CAR_Gen2.Effector_vs_Lysis(y0_Gen2)
    yD_NK_Gen2 = (Gen2_NK_data - yM_NK_Gen2)
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0_Gen2)
    yD_WT = (WT_NK_data - yM_WT)
    cost = sum(yD_NK_Gen4**2 + yD_NK_Gen2**2 + yD_WT**2)
    print('cost', cost)
    return (yM_NK_Gen4, yM_NK_Gen2, yM_WT, x0, cost)

def Residue_Pred(x0,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data,y0):
    y0[2],y0[-1] = x0[0],x0[1]
    #y0 = np.insert(y0,3,fixed_x0[0])
    #y0 = np.insert(y0,6,fixed_x0[1])
    y0_Gen4 = [elem for i, elem in enumerate(y0) if i != 1]
    yM_NK_Gen4 = Sys_CAR_Gen4.Effector_vs_Lysis(y0_Gen4)
    yD_NK_Gen4 = (Gen4_NK_data - yM_NK_Gen4)
    y0_Gen2 = [elem for i, elem in enumerate(y0) if i != 0]
    yM_NK_Gen2 = Sys_CAR_Gen2.Effector_vs_Lysis(y0_Gen2)
    yD_NK_Gen2 = (Gen2_NK_data - yM_NK_Gen2)
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0_Gen2)
    yD_WT = (WT_NK_data - yM_WT)
    print(x0.tolist())
    cost = sum(yD_WT**2)
    print('WT cost', cost)
    cost = sum(yD_NK_Gen4**2 + yD_NK_Gen2**2 + yD_WT**2)
    print('totl cost', cost)
    return 2.0*np.concatenate((yD_NK_Gen4, yD_NK_Gen2, yD_WT))#2.0*yD_WT

def Model_Pred(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data,y0,fit=True):
    start_time = time.time()
    if fit:
        result = least_squares(Residue_Pred, x0, args=(Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data,y0), bounds=(LB,UB))
        x0 = (result.x).tolist()
    y0[2],y0[-1] = x0[0],x0[1]
    y0_Gen4 = [elem for i, elem in enumerate(y0) if i != 1]
    yM_NK_Gen4 = Sys_CAR_Gen4.Effector_vs_Lysis(y0_Gen4)
    yD_NK_Gen4 = (Gen4_NK_data - yM_NK_Gen4)
    y0_Gen2 = [elem for i, elem in enumerate(y0) if i != 0]
    yM_NK_Gen2 = Sys_CAR_Gen2.Effector_vs_Lysis(y0_Gen2)
    yD_NK_Gen2 = (Gen2_NK_data - yM_NK_Gen2)
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0_Gen2)
    yD_WT = (WT_NK_data - yM_WT)
    cost = sum(yD_NK_Gen4**2 + yD_NK_Gen2**2 + yD_WT**2)
    print('cost', cost)
    return (yM_NK_Gen4, yM_NK_Gen2, yM_WT, x0, cost)
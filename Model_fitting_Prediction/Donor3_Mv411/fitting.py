from scipy.optimize import least_squares
import numpy as np
def Residue_Fit(x0,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t):
    y0_v5 = [elem for i, elem in enumerate(x0) if i != 1]
    yM_v5_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v5,t)
    yD_v5_NK  = (yM_v5_NK-v5_NK[0])#/v5_NK[1]

    y0_v6 = [elem for i, elem in enumerate(x0) if i != 0]
    yM_v6_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v6,t)
    yD_v6_NK  = (yM_v6_NK-v6_NK[0])#/v6_NK[1]
    y0 = y0_v5
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0,t)
    mean_WT, std_WT = WT_NK[0],WT_NK[1]
    yD_WT  = (yM_WT-mean_WT)#/std_WT
    cost = sum(yD_v5_NK**2+yD_v6_NK**2+ yD_WT**2)
    print(x0.tolist())
    print('cost', cost)
    return 2.0*np.concatenate((yD_v5_NK, yD_v6_NK, yD_WT))

def Model_fit(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,fit=True):
    if fit:
        result = least_squares(Residue_Fit, x0, args=(Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t), bounds=(LB,UB))
        x0 = (result.x).tolist()

    y0_v5 = [elem for i, elem in enumerate(x0) if i != 1]
    yM_v5_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v5,t)
    yD_v5_NK  = (yM_v5_NK-v5_NK[0])#/v5_NK[1]
    y0_v6 = [elem for i, elem in enumerate(x0) if i != 0]
    yM_v6_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v6,t)
    yD_v6_NK  = (yM_v6_NK-v6_NK[0])#/v6_NK[1]
    y0 = y0_v5
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0,t)
    mean_WT, std_WT = WT_NK[0],WT_NK[1]
    yD_WT  = (yM_WT-mean_WT)#/std_WT
    cost = sum(yD_v5_NK**2+yD_v6_NK**2+ yD_WT**2)
    print('final cost', cost)
    return (yM_v5_NK,yM_v6_NK, yM_WT,x0,cost)

def Residue_Pred(x0,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,y0):
    y0[-2],y0[-1] = x0[0],x0[1]
    y0_v5 = [elem for i, elem in enumerate(y0) if i != 1]
    yM_v5_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v5,t)
    yD_v5_NK  = (yM_v5_NK-v5_NK[0])#/v5_NK[1]

    y0_v6 = [elem for i, elem in enumerate(y0) if i != 0]
    yM_v6_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v6,t)
    yD_v6_NK  = (yM_v6_NK-v6_NK[0])#/v6_NK[1]
    y0 = y0_v5
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0,t)
    mean_WT, std_WT = WT_NK[0],WT_NK[1]
    yD_WT  = (yM_WT-mean_WT)#/std_WT
    cost = sum(yD_v5_NK**2+yD_v6_NK**2+ yD_WT**2)
    print(x0.tolist())
    print('cost', cost)
    return 2.0*np.concatenate((yD_v5_NK, yD_v6_NK, yD_WT))

def Model_Pred(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,y0,fit=True):
    if fit:
        result = least_squares(Residue_Pred, x0, args=(Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,y0), bounds=(LB,UB))
        x0 = (result.x).tolist()
    y0[-2],y0[-1] = x0[0],x0[1]
    y0_v5 = [elem for i, elem in enumerate(y0) if i != 1]
    yM_v5_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v5,t)
    yD_v5_NK  = (yM_v5_NK-v5_NK[0])#/v5_NK[1]
    y0_v6 = [elem for i, elem in enumerate(y0) if i != 0]
    yM_v6_NK = Sys_CAR_NK.Effector_vs_Lysis(y0_v6,t)
    yD_v6_NK  = (yM_v6_NK-v6_NK[0])#/v6_NK[1]
    y0 = y0_v5
    yM_WT = Sys_WT_NK.Effector_vs_Lysis(y0,t)
    mean_WT, std_WT = WT_NK[0],WT_NK[1]
    yD_WT  = (yM_WT-mean_WT)#/std_WT
    cost = sum(yD_v5_NK**2+yD_v6_NK**2+ yD_WT**2)
    print('final cost', cost)
    return (yM_v5_NK,yM_v6_NK, yM_WT,x0,cost)
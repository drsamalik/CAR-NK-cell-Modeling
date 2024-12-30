import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
import random
import sys
import os

# Imported the module
from Model_CAR_NK import Model_obj_CAR_NK as obj_CAR_NK
from Model_Wt_NK import Model_obj_WT_NK as obj_WT_NK
from fitting import Residue_Fit, Model_fit
from imp_exp_data import init_data, new_data, state_para, final_para

data_Mv411_48h = pd.read_excel('Average_CD33-CAR_Mv411.xlsx', header=2, sheet_name=1)
Specific_Mv411_48h = data_Mv411_48h.loc[7:,['Unnamed: 3','Unnamed: 5','Unnamed: 7']]
ET_ratio = ['2.5:1','1.25:1','0.625:1','0.313:1']
ET_ratio_num = np.arange(len(ET_ratio))

Sys_CAR_NK = obj_CAR_NK()
Sys_WT_NK = obj_WT_NK()

_,arr = init_data(Specific_Mv411_48h.iloc[:,0])

LB = np.array([1.0e-1,1.0e-1,100,1.0e-1,1.0e-1,
               1.0e-2,
               1.0e-2,1.0e-2,
               1.e-6, 1.e-4])
UB = np.array([1.0,1.0,20.0e+3,1.0,1,
               1.0e+0,
               1.0e+0,1.0e+0,
               1.e-4,1.0e-2])

outdir = sys.argv[1]
if not os.path.exists(outdir):
    os.makedirs(outdir)
out_file = os.path.join(outdir, 'est_par.csv')
totl_file = os.path.join(outdir, 'all_data.csv')


def Start_fitting(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,cost):
    x0_init = x0 + [cost]
    res0 = Model_fit(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t)
    thrs_cost = 30
    if res0[4] < thrs_cost: # if final cost is less than 20
        param = res0[3]+[res0[4]]
        print(res0[3],[res0[4]])
        print(param)
        state_para(WT_NK,v5_NK,v6_NK,x0_init,param,totl_file)
        final_para(param, out_file)
    else:
        x0 = res0[3]
        res0 = Model_fit(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t)
        if res0[4] < thrs_cost: # if final cost is less than 20     
            param = res0[3]+[res0[4]]
            state_para(WT_NK,v5_NK,v6_NK,x0_init,param,totl_file)
            final_para(param, out_file)
        elif res0[4] < 100.0:
            x0 = res0[3]
            res0 = Model_fit(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t)
            if res0[4] < thrs_cost: # if final cost is less than 20     
                param = res0[3]+[res0[4]]
                state_para(WT_NK,v5_NK,v6_NK,x0_init,param,totl_file)
                final_para(param, out_file)

Tr = 46500.0
Sys_CAR_NK.Cell_type_R_H(Targets = Tr)
Sys_WT_NK.Cell_type_R_H(Targets = Tr)
t = 48.0
for i in range(20):
    WT_NK = new_data(arr[0])
    v5_NK = new_data(arr[2])
    v6_NK = new_data(arr[3])
    cost = 1e+4
    while cost >= 1e+4:
        x0 = [random.uniform(LB[i], UB[i]) for i in range(len(LB))]
        print(x0)
        res0 = Model_fit(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,fit=False)
        cost = res0[4]
    Start_fitting(x0,LB,UB,Sys_CAR_NK,Sys_WT_NK,WT_NK,v5_NK,v6_NK,t,cost)
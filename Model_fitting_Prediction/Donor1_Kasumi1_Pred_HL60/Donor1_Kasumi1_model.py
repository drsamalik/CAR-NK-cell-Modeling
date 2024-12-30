import numpy as np
import pandas as pd
import random
import sys
import os


# Imported the module
from Model_CAR_NK import Model_obj_CAR_NK as obj_CAR_NK
from Model_Wt_NK import Model_obj_WT_NK as obj_WT_NK
from fitting import Residue_Fit, Model_fit
from imp_exp_data import new_data, state_para, final_para

Sys_CAR_Gen2 = obj_CAR_NK()
Sys_CAR_Gen4 = obj_CAR_NK()
Sys_WT_NK = obj_WT_NK()
ET_ratio = ['10:1','5:1','2.5:1','1.25:1','0.6:1']
ET_ratio_num = np.arange(len(ET_ratio))
data_Kasumi1_4h = pd.read_excel('Donor1_Av_Specific_Lysis_vs_ET.xlsx',dtype=object)
data_Kasumi1_4h = data_Kasumi1_4h.iloc[:,1:]

LB = np.array([1.0e-2,1.0e-2,1,1.0e-2,1.0e-2,
               1.0e-3,
               1.0e-3,
               1.5e-7])
UB = np.array([1.0,1.0,5.0e+3,1.0,1.0,
               1.0e-0,
               1.0e+0,
               1.e-4])



outdir = sys.argv[1]
if not os.path.exists(outdir):
    os.makedirs(outdir)
out_file = os.path.join(outdir, 'est_par.csv')
totl_file = os.path.join(outdir, 'all_data.csv')

def Start_fitting(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data):
    x0_init = x0
    res0 = Model_fit(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data)
    thrs_cost = 150
    if res0[4] < thrs_cost: # if final cost is less than 20
        param = res0[3]+[res0[4]]
        state_para(WT_NK_data,Gen2_NK_data,Gen4_NK_data,x0_init,param,totl_file)
        final_para(param, out_file)
    else:
        x0 = res0[3]
        res0 = Model_fit(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data)
        if res0[4] < thrs_cost: # if final cost is less than 20     
            param = res0[3]+[res0[4]]
            state_para(WT_NK_data,Gen2_NK_data,Gen4_NK_data,x0_init,param,totl_file)
            final_para(param, out_file)
        elif res0[4] < 200.0:
            x0 = res0[3]
            res0 = Model_fit(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data)
            if res0[4] < thrs_cost: # if final cost is less than 20     
                param = res0[3]+[res0[4]]
                state_para(WT_NK_data,Gen2_NK_data,Gen4_NK_data,x0_init,param,totl_file)
                final_para(param, out_file)
def main():
    for i in range(15):
        Sys_CAR_Gen4.Cell_type_R_H(cell_type='Kasumi1')
        Sys_CAR_Gen2.Cell_type_R_H(cell_type='Kasumi1',mult=10)
        Sys_WT_NK.Cell_type_R_H(cell_type='Kasumi1')
        WT_NK_data= new_data(data_Kasumi1_4h[['Mean_WT','SD_WT']])
        Gen2_NK_data = new_data(data_Kasumi1_4h[['Mean_CAR_NK_Gen2','SD_CAR_NK_Gen2']])
        Gen4_NK_data = new_data(data_Kasumi1_4h[['Mean_CAR_NK_Gen4','SD_CAR_NK_Gen4']])
        x0 = [random.uniform(LB[j], UB[j]) for j in range(len(LB))]
        Start_fitting(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,Gen4_NK_data,Gen2_NK_data,WT_NK_data)      
if __name__ == "__main__":
    main()
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
from sklearn.metrics import r2_score

from Model_CAR_NK import Model_obj_CAR_NK as obj_CAR_NK
from Model_Wt_NK import Model_obj_WT_NK as obj_WT_NK
from fitting import Model_fit, Model_Pred
from imp_exp_data import new_data, state_para, final_para

Sys_CAR_Gen2 = obj_CAR_NK()
Sys_CAR_Gen4 = obj_CAR_NK()
Sys_WT_NK = obj_WT_NK()
ET_ratio = ['10:1','5:1','2.5:1','1.25:1','0.6:1']
ET_ratio_num = np.arange(len(ET_ratio))

data_Kasumi1_4h = pd.read_excel('Donor1_Av_Specific_Lysis_vs_ET.xlsx',dtype=object)
data_Kasumi1_4h = data_Kasumi1_4h.iloc[:,1:]
mean_WT_Kasumi1, sd_WT_Kasumi1 = new_data(data_Kasumi1_4h[['Mean_WT','SD_WT']])
mean_Gen2_Kasumi1, sd_Gen2_Kasumi1 = new_data(data_Kasumi1_4h[['Mean_CAR_NK_Gen2','SD_CAR_NK_Gen2']])
mean_Gen4_Kasumi1, sd_Gen4_Kasumi1 = new_data(data_Kasumi1_4h[['Mean_CAR_NK_Gen4','SD_CAR_NK_Gen4']])

LB = np.array([1.0e-2,1.0e-2,500,0.01,1.0e-2,
               1.0e-3,
               1.5e-3,1.0e-3,
               1.5e-7])
UB = np.array([0.20,0.26,5.0e+3,0.1,1.0,
               1.0e-0,
               1.5e-1,1.0e+0,
               1.e-4])

def main_Kasumi1():
    Sys_CAR_Gen4.Cell_type_R_H(cell_type='Kasumi1')
    Sys_CAR_Gen2.Cell_type_R_H(cell_type='Kasumi1',mult=10)
    Sys_WT_NK.Cell_type_R_H(cell_type='Kasumi1')
    x0 = [0.18510612989064024, 0.1668436093796393, 612.6622265256295, 0.03117145473112584, 0.43729631556909043, 0.5787585107835337, 0.001500751222309375, 0.04616750953248207, 1.911281041757332e-05]
    res0 = Model_fit(x0,LB,UB,Sys_CAR_Gen4,Sys_CAR_Gen2,Sys_WT_NK,mean_Gen4_Kasumi1,mean_Gen2_Kasumi1,mean_WT_Kasumi1,fit=False)
    return res0

res0 = main_Kasumi1()
print('saeed')
# plt.figure(figsize=(10,8))
# plt.rcParams['axes.linewidth'] = 8
# ax = plt.gca()
# ax.tick_params(axis='both', which='major', labelsize=19, length=8, width=5)
# ls = [(1,(1,1)),(5, (10, 3)), (0,(2,2)), (0, (3, 1, 1, 1))]
# m_size = 20
# lw = 4.5
# plt.plot(ET_ratio_num, res0[0], marker='o', markersize=m_size, markerfacecolor='white', markeredgewidth=3.5, markeredgecolor='darkcyan',color = 'darkcyan',lw=lw,ls = ls[0])#,label="Fit-Gen4")
# plt.plot(ET_ratio_num, res0[1], marker='o', markersize=m_size, markerfacecolor='white', markeredgewidth=3.5, markeredgecolor='blue',color = 'blue',lw=lw,ls = ls[0])#,label="Fit-Gen2")
# plt.plot(ET_ratio_num, res0[2], marker='o', markersize=m_size, markerfacecolor='white', markeredgewidth=3.5, markeredgecolor='gray',color = 'gray',lw=lw,ls = ls[0])#,label="Fit-WT")

# plt.errorbar(ET_ratio_num, mean_Gen4_Kasumi1, yerr = sd_Gen4_Kasumi1, 
#              fmt='^', markersize=m_size, markerfacecolor='lightgreen', markeredgewidth=3.5, markeredgecolor='darkcyan',
#              elinewidth=3.5, capsize=10, capthick=20,
#              ecolor='darkcyan',alpha=0.99)#,label="data-Gen4")
# plt.errorbar(ET_ratio_num, mean_Gen2_Kasumi1, yerr = sd_Gen2_Kasumi1, 
#              fmt='^', markersize=m_size, markerfacecolor='lightblue', markeredgewidth=3.5, markeredgecolor='blue',
#              elinewidth=3.5, capsize=10, capthick=20,
#              ecolor='blue',alpha=0.99)#,label="data-Gen2")

# plt.errorbar(ET_ratio_num, mean_WT_Kasumi1, yerr = sd_WT_Kasumi1, 
#              fmt='^', markersize=m_size, markerfacecolor='lightgray', markeredgewidth=3.5, markeredgecolor='gray',
#              elinewidth=3.5, capsize=10, capthick=20,
#              ecolor='gray',alpha=0.99)#,label="data-WT")
# t_size = 45
# plt.xticks(ET_ratio_num,ET_ratio,fontname="Arial",fontsize = t_size)
# plt.yticks([10,40,70],fontname="Arial",fontsize = t_size)
# #plt.tight_layout(pad=2.5)
# #plt.legend(bbox_to_anchor=(0.95, 1),fontsize=17, loc='upper right', labelcolor='white')
# plt.savefig('Kasumi1_Donor1_ET_vs_AvSpecificLysis_t_4h.pdf')
# plt.show()
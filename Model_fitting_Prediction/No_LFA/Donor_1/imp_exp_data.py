import numpy as np
import pandas as pd
import random

def new_data(df):
    mean,sd = df.iloc[:,0],df.iloc[:,1]#/np.sqrt(3)
    a = np.random.normal(mean,sd,(10,5))
    b = [np.random.choice(a[:,i][(a[:,i]>=mean[i]-sd[i])&(a[:,i]<=mean[i]+sd[i])]) for i in range(len(df))]
    return (mean.values).tolist(),(sd.values).tolist() #b #(mean.values).tolist()
def final_para(param, file):
    with open(file,'a') as f:
        f.write(",".join(map(str, param))+"\n")
def state_para(WT_NK_data,Gen2_NK_data,Gen4_NK_data,Init_x0,Final_x0,file):
    with open(file, 'a') as f:
        f.write("#WT_NK:\n" + ",".join(map(str, WT_NK_data)))
        f.write("\n#Gen2_NK:\n" + ",".join(map(str, Gen2_NK_data)))
        f.write("\n#Gen4_NK:\n" + ",".join(map(str, Gen4_NK_data)))
        f.write("\n#Init_x0:\n" + ",".join(map(str, Init_x0)))
        f.write("\n#Final_x0:\n" + ",".join(map(str, Final_x0))+"\n")
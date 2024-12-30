import numpy as np
import pandas as pd
from ast import literal_eval
def init_data(df):
    arr = df.values.reshape(-1,3)
    _lysis = [arr[i:i+4] for i in range(0,len(arr),4)]
    return arr,_lysis

def new_data(x):
    # generate boostrapped data with mean std 
    mean_list = []
    sigma_list = []
    for i in range(len(x)):
        lysis_data = x[i]
        if isinstance(lysis_data[0], str):
            lysis_data = np.delete(lysis_data,0)
        random_elements = lysis_data #np.random.choice(lysis_data, size=len(lysis_data), replace=True)
        mean = np.mean(random_elements)
        sigma = np.std(random_elements)
        '''
        min_ = np.std(lysis_data)
        min_sigma = 1 if min_ < 1 else min_
        if sigma < min_sigma:
            sigma = min_sigma  # Use the minimum sigma value if below threshold
        '''
        mean_list.append(mean)
        sigma_list.append(sigma)
    return (np.flip(np.array(mean_list)), np.flip(np.array(sigma_list)))

def final_para(param, file):
    with open(file,'a') as f:
        f.write(",".join(map(str, param))+"\n")
def state_para(WT_NK_data,Gen2_NK_data,Gen4_NK_data,Init_x0,Final_x0,file):
    with open(file, 'a') as f:
        f.write("#WT_NK:\n" + ",".join(map(str, WT_NK_data[0])))
        f.write("\n" + ",".join(map(str, WT_NK_data[1])))
        f.write("\n#Gen2_NK:\n" + ",".join(map(str, Gen2_NK_data[0])))
        f.write("\n" + ",".join(map(str, Gen2_NK_data[1])))
        f.write("\n#Gen4_NK:\n" + ",".join(map(str, Gen4_NK_data[0])))
        f.write("\n" + ",".join(map(str, Gen4_NK_data[1])))
        f.write("\n#Init_x0:\n" + ",".join(map(str, Init_x0)))
        f.write("\n#Final_x0:\n" + ",".join(map(str, Final_x0))+"\n")
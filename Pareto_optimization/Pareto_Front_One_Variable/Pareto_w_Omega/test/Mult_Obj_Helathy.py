import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import truncnorm

class Main_Model():
    def __init__(self):
        pass
    def Recp_Ligand_exp(self,NK_type,frac=0.6):
        R,L = R_H_Expression(NK_type,frac)
        self.NK_cell = pd.DataFrame(columns=['E_R', 'R1', 'R3', 'R4'])
        Effector = 1.0
        for index1, row1 in R[0].iterrows():
            for index3, row3 in R[1].iterrows():
                for index4, row4 in R[2].iterrows():
                    E_R = (row1['E_R1']*row3['E_R3']*row4['E_R4']) * Effector
                    self.NK_cell.loc[len(self.NK_cell)] = [E_R, row1['R1'], row3['R3'], row4['R4']]
        Target = 1.0
        self.Target_cell = pd.DataFrame(columns=['T_L', 'L1', 'L3', 'L4'])
        for index1, row1 in L[0].iterrows():
            for index3, row3 in L[1].iterrows():
                for index4, row4 in L[2].iterrows():
                    T_L = (row1['T_L1']*row3['T_L3']*row4['T_L4']) * Target
                    self.Target_cell.loc[len(self.Target_cell)] = [T_L, row1['L1'], row3['L3'], row4['L4']]
    def Solve_eqn(self,Tc_No):
        y0 = (self.Target_cell['T_L'].values)*Tc_No
        options = {
                'rtol': 1e-2,
                'atol': 1e-4,
                }
        tspan = [0, 100]
        teval = [0,4.0]
        sol = solve_ivp(odefcnCONS_LD, tspan, y0, method = 'RK45', t_eval = teval, dense_output = True, events = None, **options, 
                        args = (self.CAR_NK, self.Lam_CAR, len(y0), len(self.rho)))
        Target_cell_w_time = pd.DataFrame(sol.y)
        yF = (1-(Target_cell_w_time[Target_cell_w_time.shape[1]-1].sum()) / (Target_cell_w_time[0].sum())) * 100
        return yF
    def Healthy(self, a1,w, alpha4, NK_No, Tc_No):
        self.Recp_Ligand_exp(NK_type='CAR_NK',frac=w)
        self.CAR_NK = self.NK_cell['E_R'].values*NK_No
        self.Lam_CAR, self.rho = Rho_lambda(a1,alpha4,self.NK_cell,self.Target_cell)
        YY_H = self.Solve_eqn(Tc_No)
        print('H',YY_H)
        return YY_H
    
def odefcnCONS_LD(tR, U_H0, T_R0, Lam_CAR, lenH, lenR):        
    dUHdt = np.zeros(lenH)
    for Hi in range(lenH):
        Lam_CAR_NK = np.zeros(lenR)
        for Ri in range(lenR):
            Lam_CAR_NK[Ri] = Lam_CAR[Hi,Ri] * T_R0[Ri]
        dUHdt[Hi] = -np.sum(Lam_CAR_NK)*U_H0[Hi]
    return dUHdt

def Rho_lambda(a1,alpha4,NK_cell,Target_cell, WT = False):
        x0 = np.array([9.99935926e-01, 3.91486705e+03, 8.38786722e-01, 8.29549637e-01,
       1.75149471e-01, 6.74554036e-01, 5.76540111e-03, 1.49985099e-04])
        x = np.array(x0)
        KD1 = 1.46 * (0.6 * 113 * 0.002)
        KD3 = 50.0 * (0.6 * 113 * 0.002)
        KD4 = 36.0 * (0.6 * 113 * 0.002)
        alph1 = a1
        if WT:
            alph1 = 0 
        C2N = x[1]
        alph3 = x[2]
        alph4 = alpha4
        Vc = x[4]
        K1 = x[5]
        K2 = x[6]
        k = x[7]
        Kr = K1/K2
        R_set = NK_cell[['R1','R3', 'R4']]
        L_set = Target_cell[['L1','L3', 'L4']]
        Lambda = np.zeros((len(L_set), len(R_set)))
        for index1, row1 in L_set.iterrows():
            for index2, row2 in R_set.iterrows():
                # C0 complex of R H interaction
                C10 = 0.5 * (row2['R1'] + row1['L1'] + KD1) * (1 - np.sqrt(1 - (4 * row2['R1']*row1['L1'] / ((row2['R1'] + row1['L1'] + KD1) ** 2))))
                C30 = 0.5 * (row2['R3'] + row1['L3'] + KD3) * (1 - np.sqrt(1 - (4 * row2['R3']*row1['L3'] / ((row2['R3'] + row1['L3'] + KD3) ** 2))))
                C40 = 0.5 * (row2['R4'] + row1['L4'] + KD4) * (1 - np.sqrt(1 - (4 * row2['R4']*row1['L4'] / ((row2['R4'] + row1['L4'] + KD4) ** 2))))
                C1N = alph1 * C10
                C3N = alph3 * C30
                C4N = alph4 * C40
                Vr = Vc*(C1N + C2N + C3N)/C4N
                W2 = ((Vr - 1) - K2 * (Kr + Vr) + np.sqrt((Vr - 1 - K2 * (Kr + Vr)) ** 2 + 4 * K2 * (Vr - 1) * Vr)) / (2 * (Vr - 1))
                LamX = k * W2
                Lambda[index1][index2] = LamX
        return Lambda, Lambda.T

def R_H_Expression(NK_cell,frac):
    # Receptor Lygand Interaction
    dat0 = pd.read_excel('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Mono_CAR_NK_CD33.xlsx')
    dat10 = pd.read_excel('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Mono_LFA1_ICAM1.xlsx',sheet_name=0)
    dat11 = pd.read_excel('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Mono_LFA1_ICAM1.xlsx',sheet_name=1)
    dat2 = pd.read_excel('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Mono_KIR2DL1_HLA.xlsx')
    if NK_cell == 'CAR_NK':
        R1_ER1 = CAR_data(frac)#dat0.loc[0:4,['No_of_CAR_NK','prob_Car_NK']].rename(columns= {'No_of_CAR_NK' : 'R1','prob_Car_NK' : 'E_R1'}, inplace = False)
        R3_ER3 = dat10.loc[0:4,['No_of_LFA1','prob_LFA1']].rename(columns= {'No_of_LFA1' : 'R3','prob_LFA1' : 'E_R3'}, inplace = False)
        R4_ER4 = dat2.loc[0:4,['No_of_KIR2DL1','prob_KIR']].rename(columns= {'No_of_KIR2DL1' : 'R4','prob_KIR' : 'E_R4'}, inplace = False)
    elif NK_cell == 'WT_NK':
        R1_ER1 = CAR_data(frac)#dat0.loc[0:4,['No_of_CAR_NK','prob_Car_NK']].rename(columns= {'No_of_CAR_NK' : 'R1','prob_Car_NK' : 'E_R1'}, inplace = False)
        R3_ER3 = dat11.loc[0:4,['No_of_LFA1','prob_LFA1']].rename(columns= {'No_of_LFA1' : 'R3','prob_LFA1' : 'E_R3'}, inplace = False)
        R4_ER4 = dat2.loc[0:4,['No_of_KIR2DL1','prob_KIR']].rename(columns= {'No_of_KIR2DL1' : 'R4','prob_KIR' : 'E_R4'}, inplace = False) 
    L1_TL1 = dat0.loc[0:4,['No_of_CD33','prob_CD33']].rename(columns= {'No_of_CD33' : 'L1','prob_CD33' : 'T_L1'}, inplace = False)
    L3_TL3 = dat10.loc[0:4,['No_of_ICAM','prob_ICAM']].rename(columns= {'No_of_ICAM' : 'L3','prob_ICAM' : 'T_L3'}, inplace = False)
    L4_TL4 = dat2.loc[0:4,['No_of_HLA','prob_HLA']].rename(columns= {'No_of_HLA' : 'L4','prob_HLA' : 'T_L4'}, inplace = False)
    return [R1_ER1, R3_ER3, R4_ER4], [L1_TL1,L3_TL3,L4_TL4]

def adding_prob_dens(df, param, w=0.54, bins=4):
    np.random.seed(0)
    num = df['MFI'].values
    num = np.insert(num,0,0)
    step = 1
    num_indx = (np.linspace(0,len(num)-1,int(len(num)/step))).astype(int)
    nmbr = []
    for i in num_indx:
        nmbr.append(num[i])
    num = np.array(nmbr)
    size = df['Count'].sum()
    mu1, sigma1, mu2, sigma2 = param
    lb1, ub1 = df['MFI'].min(), 1.14
    lb1 = (np.log(lb1) - mu1) / sigma1
    ub1 = (np.log(ub1) - mu1) / sigma1
    lb2, ub2 = 1.14, df['MFI'].max()
    lb2 = (np.log(lb2) - mu2) / sigma2
    data1 = np.exp(truncnorm.rvs(a=lb1, b=ub1, loc=mu1, scale=sigma1,size=int(size * (1-w))))
    data2 = np.exp(truncnorm.rvs(a=lb2, b=ub2, loc=mu2, scale=sigma2,size=int(size * w)))
    if w==1:
        counts, num = np.histogram(data2, bins=bins+1)
        MFI = (num[1:] + num[:-1]) / 2
    else:
        num_wt = np.array([len(data1)])
        MFI_wt = np.array([(data1[0]+data1[-1])/2])
        counts, num = np.histogram(data2, bins=bins)
        MFI = (num[1:] + num[:-1]) / 2
        MFI = np.concatenate((MFI_wt,MFI))
        counts = np.concatenate((num_wt,counts))
    prob = counts/sum(counts)
    Nmbr_mcle = 10**(0.9837*np.log10(MFI) + 3.0402)
    return Nmbr_mcle,prob
def CAR_data(frac=0.54):
    sheet_names = (pd.ExcelFile('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Histograms_for_Pareto.xlsx')).sheet_names
    xls = pd.ExcelFile('/Users/sxa126/Dropbox/CAR_NK_model/CAR_NK_Python/Pareto_Opt/Pure_Mono/data/Histograms_for_Pareto.xlsx')
    df = (xls.parse(sheet_names[6]))[['MFI', 'Count']]
    df = df[(df['MFI']>0) & (df['MFI']<150)]
    param = np.array([0.04131501, 1.05464126, 3.1699362 , 1.50408986])
    Nmbr_mcle,prob = adding_prob_dens(df,param, w=frac)
    dist_df=pd.DataFrame(columns=[['R1','E_R1']])
    dist_df['R1'] = Nmbr_mcle
    dist_df['E_R1'] =prob
    return dist_df
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
import time
from multiprocessing import Pool
class Model_obj_CAR_NK:
    def __init__(self):
        self.KD1 = 1.46 * (0.6 * 113 * 0.002)
        self.WT = 104000.0
    def Cell_type_R_H(self,E=1, cell_type=None,mult=1):
        R,H = R_H_Expression(cell_type,mult)
        self.NK_cell = pd.DataFrame(columns=['T_R', 'R1'])
        Effector = E
        for index1, row1 in R[0].iterrows():
            T_R = (row1['T_R1']) * Effector
            self.NK_cell.loc[len(self.NK_cell)] = [T_R, row1['R1']]
        Target = 10000.0
        self.Target_cell = pd.DataFrame(columns=['U_H', 'H1'])
        for index1, row1 in H[0].iterrows():
            U_H = (row1['U_H1']) * Target
            self.Target_cell.loc[len(self.Target_cell)] = [U_H, row1['H1']]
    def Rho_Lambda(self,x0):
        R_set = self.NK_cell[['R1']]
        H_set = self.Target_cell[['H1']]
        # Unknown para
        # R1 H1
        alph1 = x0[0]
        # Activating
        C2N = x0[1]
        #Vav to pVav
        Vc = x0[2]
        K1 = x0[3]
        K2 = x0[4]
        Kr = K1 / K2
        # cell rates
        self.k = x0[5]
        Lambda = np.zeros((len(H_set), len(R_set)))
        for index1, row1 in H_set.iterrows():
            for index2, row2 in R_set.iterrows():
                # C0 complex of R H interaction
                C10 = 0.5 * (row2['R1'] + row1['H1'] + self.KD1) * (1 - np.sqrt(1 - (4 * row2['R1']*row1['H1'] / ((row2['R1'] + row1['H1'] + self.KD1) ** 2))))
                C1N = alph1 * C10
                Vr = Vc*(C1N + C2N)
                W2 = ((Vr - 1) - K2 * (Kr + Vr) + np.sqrt((Vr - 1 - K2 * (Kr + Vr)) ** 2 + 4 * K2 * (Vr - 1) * Vr)) / (2 * (Vr - 1))
                LamX = self.k * W2
                Lambda[index1][index2] = LamX
        self.Lam = Lambda
        self.rho = Lambda.T # lenR

    def Sol_Lysis(self,tD,y0,y_NK_Cell):
        options = {
                'rtol': 1e-2,
                'atol': 1e-4,
                }
        tspan = [0,100]
        teval = [0,tD]
        sol = solve_ivp(odefcnCONS_LD, tspan, y0, method = 'RK45', t_eval = teval, dense_output = True, events = None, **options, args = (y_NK_Cell, self.Lam, len(y0), len(self.rho)))
        Target_cell_w_time = pd.DataFrame(sol.y) # Column correspond to time axis
        yF = (1-(Target_cell_w_time[Target_cell_w_time.shape[1]-1].sum()) /(Target_cell_w_time[0].sum())) * 100
        return yF 
    def Effector_vs_Lysis(self,x0):
        self.y0 = self.Target_cell['U_H'].values #
        self.y_NK_Cell = self.NK_cell['T_R'].values #
        Effector = [100000,50000,25000,12500,6200]
        self.tD = 4.0
        Specific_lysis = []
        self.Rho_Lambda(x0)
        for i,E in enumerate(Effector):
           yF = self.Sol_Lysis(self.tD,self.y0,self.y_NK_Cell*E)
           Specific_lysis.append(yF)
        return np.array(Specific_lysis)
        
def odefcnCONS_LD(tR, U_H0, T_R0, Lam, lenH, lenR):        
    dUHdt = np.zeros(lenH)
    for Hi in range(lenH):
        Lam_T_R = np.zeros(lenR)
        for Ri in range(lenR):
            Lam_T_R[Ri] = Lam[Hi,Ri] * T_R0[Ri]
        dUHdt[Hi] = -np.sum(Lam_T_R)*U_H0[Hi]
    return dUHdt

def R_H_Expression(cell_type,mult):
    # Receptor Lygand Interaction
    if cell_type == 'Kasumi1':
        dat1 = pd.read_excel('Kasumi1_CAR_NK_CD33.xlsx')
        dat3 = pd.read_excel('Kasumi1_LFA1_ICAM1.xlsx')
        dat4 = pd.read_excel('Kasumi1_KIR2DL1_HLA.xlsx')
    else:
        dat1 = pd.read_excel('HL60_CAR_NK_CD33.xlsx')
        dat3 = pd.read_excel('HL60_LFA1_ICAM1.xlsx')
        dat4 = pd.read_excel('HL60_KIR2DL1_HLA.xlsx')
    R1_TR1 = dat1.loc[0:6,['No_of_CAR_NK','prob_Car_NK']].rename(columns= {'No_of_CAR_NK' : 'R1','prob_Car_NK' : 'T_R1'}, inplace = False)
    R1_TR1['R1'] = R1_TR1['R1']*mult
    H1_UH1 = dat1.loc[0:4,['No_of_CD33','prob_CD33']].rename(columns= {'No_of_CD33' : 'H1','prob_CD33' : 'U_H1'}, inplace = False)
    R3_TR3 = dat3.loc[0:5,['No_LFA1','prob_LFA1']].rename(columns= {'No_LFA1' : 'R3','prob_LFA1' : 'T_R3'}, inplace = False)
    H3_UH3 = dat3.loc[0:4,['No_of_ICAM','prob_ICAM']].rename(columns= {'No_of_ICAM' : 'H3','prob_ICAM' : 'U_H3'}, inplace = False)
    R4_TR4 = dat4.loc[0:4,['No_of_KIR2DL1','prob_KIR']].rename(columns= {'No_of_KIR2DL1' : 'R4','prob_KIR' : 'T_R4'}, inplace = False)
    H4_UH4 = dat4.loc[0:4,['No_MHC1','prob_MHC1']].rename(columns= {'No_MHC1' : 'H4','prob_MHC1' : 'U_H4'}, inplace = False)
    return [R1_TR1, R3_TR3, R4_TR4], [H1_UH1,H3_UH3,H4_UH4]
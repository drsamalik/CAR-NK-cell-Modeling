import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from multiprocessing import Pool
class Model_obj_WT_NK:
    def __init__(self):
        pass
    def Rho_Lambda(self,x0):
        C2N = x0[1]
        #Vav to pVav
        Vc = x0[2]
        K1 = x0[3]
        K2 = x0[4]
        Kr = K1 / K2
        # cell rates
        self.k = x0[5]
        Vr = Vc*(C2N)
        W2 = ((Vr - 1) - K2 * (Kr + Vr) + np.sqrt((Vr - 1 - K2 * (Kr + Vr)) ** 2 + 4 * K2 * (Vr - 1) * Vr)) / (2 * (Vr - 1))
        LamX = self.k * W2
        Lambda = LamX
        print(Lambda)
        self.Lam = Lambda

    def Sol_Lysis(self,tD,y0,y_NK_Cell):
        options = {
                'rtol': 1e-2,
                'atol': 1e-4,
                }
        tspan = [0,100]
        teval = [0,tD]
        sol = solve_ivp(odefcnCONS_LD, tspan, [y0], method = 'RK45', t_eval = teval, dense_output = True, events = None, **options, args = (y_NK_Cell, self.Lam))
        Target_cell_w_time = sol.y[0]
        yF = (1-(Target_cell_w_time[-1])/(Target_cell_w_time[0])) * 100
        return yF
    def Effector_vs_Lysis(self,x0):
        self.y0 = 10000.0
        self.y_NK_Cell = 1.0
        Effector = [100000,50000,25000,12500,6200]
        self.tD = 4.0
        Specific_lysis = []
        self.Rho_Lambda(x0)
        for i,E in enumerate(Effector):
            yF = self.Sol_Lysis(self.tD,self.y0,self.y_NK_Cell*E)
            Specific_lysis.append(yF)
        return np.array(Specific_lysis) 
def odefcnCONS_LD(tR, U_H0, T_R0, Lam):  
    return -Lam*T_R0*U_H0
def R_H_Expression(cell_type = None):
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
    H1_UH1 = dat1.loc[0:4,['No_of_CD33','prob_CD33']].rename(columns= {'No_of_CD33' : 'H1','prob_CD33' : 'U_H1'}, inplace = False)
    R3_TR3 = dat3.loc[0:5,['No_LFA1','prob_LFA1']].rename(columns= {'No_LFA1' : 'R3','prob_LFA1' : 'T_R3'}, inplace = False)
    H3_UH3 = dat3.loc[0:4,['No_of_ICAM','prob_ICAM']].rename(columns= {'No_of_ICAM' : 'H3','prob_ICAM' : 'U_H3'}, inplace = False)
    R4_TR4 = dat4.loc[0:4,['No_of_KIR2DL1','prob_KIR']].rename(columns= {'No_of_KIR2DL1' : 'R4','prob_KIR' : 'T_R4'}, inplace = False)
    H4_UH4 = dat4.loc[0:4,['No_MHC1','prob_MHC1']].rename(columns= {'No_MHC1' : 'H4','prob_MHC1' : 'U_H4'}, inplace = False)
    return [R1_TR1, R3_TR3, R4_TR4], [H1_UH1,H3_UH3,H4_UH4]
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, simps, ode
from math import exp, ceil
import pdb
#----------------------------------------------------------------------------------------------
# from parameters_HRHR import T_emerge,params['T_GS32'], params['T_GS39'], params['T_GS61'], params['T_GS87'], params['init_den'],  , den_frac, S_0, params['dt'], nstepz, no_variables,S_ind,ER_ind,ERS_ind,ESR_ind,ES_ind,IR_ind,IRS_ind,ISR_ind,IS_ind,R_ind,PR_ind,PRS_ind,PSR_ind,PS_ind,Fung1_ind,Fung2_ind # r, k, omega_1, theta_1, omega_2, theta_2, beta, gamma, mu, nu, params['delta_1'], params['delta_2']
# r,k,omega_1,theta_1,omega_1_L,params['theta_1_L'],omega_2,theta_2,omega_2_L,beta,gamma,mu,nu,params['delta_1'],params['delta_2'],alpha_1,alpha_2,alpha_1_C,alpha_2_C,params['init_den'],den_frac,params['innoc_frac'],params['innoc_frac_integral'], params['dt'], nstepz, no_variables, S_ind,ER_ind,ERS_ind,ESR_ind,ES_ind,IR_ind,IRS_ind,ISR_ind,IS_ind,R_ind,PR_ind,PRS_ind,PSR_ind,PS_ind,Fung1_ind,Fung2_ind, sex_prop #,ratio,ratio2 #t_points
# will be errors in here because of params not sorted through
from parameters_HRHR import params
from functions_HRHR import fngcide, growth, Gamma, solv_ode, calculator_d_free, calculator_d_free_old, lifetime_yield, Total_yield,inoculum_value,resist_prop_calculator,primary_calculator
from calculator_HRHR import calculator_ode
#----------------------------------------------------------------------------------------------
def calculator_ode_int(dose_11, dose_12, dose_21, dose_22,res_prop_1,res_prop_2,innoc,res_prop_calculator,phi_rr=None,phi_rs=None,phi_sr=None,phi_ss=None):
    t    = np.arange(params['T_emerge'],params['T_GS32']+1,1)
    t1   = np.arange(params['T_GS32'],params['T_GS39']+1,1)
    t2   = np.arange(params['T_GS39'],params['T_GS61']+1,1)
    t3   = np.arange(params['T_GS61'],params['T_GS87']+1,1)
#----------------------------------------------------------------------------------------------    
    P_RR0 , P_RS0 , P_SR0, P_SS0 = primary_calculator(res_prop_1,res_prop_2,innoc,phi_rr,phi_rs,phi_sr,phi_ss)
    # y  = [S, ERR,ERS,ESR,ESS,IRR,IRS,ISR,ISS,R,PR,PRS,PSR,PS,Fung1,Fung2]
    #S_ind,ER_ind,ERS_ind,ESR_ind,ES_ind,IR_ind,IRS_ind,ISR_ind,IS_ind,R_ind,PR_ind,PRS_ind,PSR_ind,PS_ind,Fung1_ind,Fung2_ind
    y0   = [params['S_0'],0,  0,  0, 0, 0,  0,  0, 0,0, P_RR0 , P_RS0 , P_SR0, P_SS0,    0,   0]
    sol  = odeint(solv_ode_int,y0,t)#,rtol=0.3*10**(-13),atol=10**(-25))#,args = (r,k,omega_L,theta_L,omega_H,theta_H,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,params['delta_1'],params['delta_2']))
#----------------------------------------------------------------------------------------------
    y1    = [sol[- 1,S_ind],    sol[- 1,ER_ind],    sol[- 1,ERS_ind],    sol[- 1,ESR_ind],    sol[- 1,ES_ind],   sol[- 1,IR_ind],  sol[- 1,IRS_ind],  sol[- 1,ISR_ind],  sol[- 1,IS_ind],  sol[- 1,R_ind],   sol[- 1,PR_ind],  sol[- 1,PRS_ind],  sol[- 1,PSR_ind],  sol[- 1,PS_ind],  sol[- 1,Fung1_ind] +  dose_11,  sol[- 1,Fung2_ind] +  dose_12]
    sol1  = odeint(solv_ode_int,y1,t1,rtol=0.3*10**(-13),atol=10**(-25))#,args = (r,k,omega_L,theta_L,omega_H,theta_H,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,params['delta_1'],params['delta_2']))
#----------------------------------------------------------------------------------------------
    y2   = [sol1[- 1,S_ind],    sol1[- 1,ER_ind],    sol1[- 1,ERS_ind],    sol1[- 1,ESR_ind],    sol1[- 1,ES_ind],   sol1[- 1,IR_ind],  sol1[- 1,IRS_ind],  sol1[- 1,ISR_ind],  sol1[- 1,IS_ind],  sol1[- 1,R_ind],   sol1[- 1,PR_ind],  sol1[- 1,PRS_ind],  sol1[- 1,PSR_ind],  sol1[- 1,PS_ind],  sol1[- 1,Fung1_ind] +  dose_21,  sol1[- 1,Fung2_ind] +  dose_22]
    sol2 = odeint(solv_ode_int,y2,t2)#,rtol=0.3*10**(-13),atol=10**(-25))#,args = (r,k,omega_L,theta_L,omega_H,theta_H,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,params['delta_1'],params['delta_2']))
#----------------------------------------------------------------------------------------------
    y3   = [sol2[- 1,S_ind],    sol2[- 1,ER_ind],    sol2[- 1,ERS_ind],    sol2[- 1,ESR_ind],    sol2[- 1,ES_ind],   sol2[- 1,IR_ind],  sol2[- 1,IRS_ind],  sol2[- 1,ISR_ind],  sol2[- 1,IS_ind],  sol2[- 1,R_ind],   sol2[- 1,PR_ind],  sol2[- 1,PRS_ind],  sol2[- 1,PSR_ind],  sol2[- 1,PS_ind],  sol2[- 1,Fung1_ind],  sol2[- 1,Fung2_ind]]
    sol3 = odeint(solv_ode_int,y3,t3)#,rtol=0.3*10**(-13),atol=10**(-25))#,args = (r,k,omega_L,theta_L,omega_H,theta_H,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,params['delta_1'],params['delta_2']))
#----------------------------------------------------------------------------------------------
    tg     = t[:-1]
    t1g    = t1[:-1]
    t2g    = t2[:-1]
    solg   = sol[:-1]
    sol1g  = sol1[:-1]
    sol2g  = sol2[:-1]
#----------------------------------------------------------------------------------------------
    solutiont = np.concatenate((tg,t1g,t2g,t3))
    solution  = np.concatenate((solg,sol1g,sol2g,sol3),axis=0)
    inte = simps(sol3[:,S_ind]+sol3[:,ER_ind]+sol3[:,ERS_ind]+sol3[:,ESR_ind]+sol3[:,ES_ind],t3)
#----------------------------------------------------------------------------------------------
    res_prop_new_1, res_prop_new_2, p_rr, p_rs, p_sr, p_ss = resist_prop_calculator(np.transpose(solution),np.transpose(solutiont),FV_or_int=res_prop_calculator)
    innoc = inoculum_value(np.transpose(solution),np.transpose(solutiont),FV_or_int=res_prop_calculator)
    if res_prop_1 > 0:
        Selection_1    = res_prop_new_1/res_prop_1
    else:
        Selection_1 = 1
    if res_prop_2 > 0:
        Selection_2    = res_prop_new_2/res_prop_2
    else:
        Selection_2 = 1
    return Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, p_rr, p_rs, p_sr, p_ss, inte, solution, solutiont

#----------------------------------------------------------------------------------------------  
def calculator_d_free_old():
    # y    = [S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2]
    y0_df   = [params['S_0'],0,0,  0, 0, 0,  0,  0, 0,0, 0,  0,  0, 0,  0  ,  0  ]
    tdf    = np.arange(params['T_emerge'],params['T_GS61']+1,1)
    soldf  = odeint(solv_ode_int,y0_df,tdf)#,rtol=10**(-10),atol=10**(-20))#,args = (r,k,omega_1,theta_1,omega_2,theta_2,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,delta_H,delta_L))
#----------------------------------------------------------------------------------------------
    y2df   = [soldf[- 1,0],    soldf[- 1,1],    soldf[- 1,2],    soldf[- 1,3],    soldf[- 1,4],   soldf[- 1,5],  soldf[- 1,6],  soldf[- 1,7],   soldf[- 1,8],soldf[- 1,9],soldf[- 1,10],soldf[- 1,11],soldf[- 1,12],soldf[- 1,13],soldf[- 1,14],soldf[- 1,15]]
    t2df   = np.arange(params['T_GS61'],params['T_GS87']+1,1)
    sol2df = odeint(solv_ode_int,y2df,t2df,rtol=10**(-10),atol=10**(-20))
#----------------------------------------------------------------------------------------------
    int_df = simps(sol2df[:,0],t2df)
    return int_df
#----------------------------------------------------------------------------------------------
#def solv(y,t,r,k,omega_1,theta_1,omega_2,theta_2,params['T_GS61'],params['T_GS87'],beta,gamma,mu,nu,params['delta_1'],params['delta_2']):
# for calculator_ode_int - replaced by calculator ode
def solv_ode_int(y,t):
    S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2 = y
    A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R
    dydt = [growth(A,t) - (Gamma(params['T_GS61'],params['T_GS87'],t))*S -  S * (beta/A) * (  
               (fngcide(        omega_1,          params['theta_1'],Fung1)) * (fngcide(        omega_2,          theta_2,Fung2)) * (IS + PS)
             + (fngcide(alpha_1*omega_1,alpha_1_C*params['theta_1'],Fung1)) * (fngcide(        omega_2,          theta_2,Fung2)) * (IRS+PRS)
             + (fngcide(        omega_1,          params['theta_1'],Fung1)) * (fngcide(alpha_2*omega_2,alpha_2_C*theta_2,Fung2)) * (ISR + PSR)
             + (fngcide(alpha_1*omega_1,alpha_1_C*params['theta_1'],Fung1)) * (fngcide(alpha_2*omega_2,alpha_2_C*theta_2,Fung2)) * (IR + PR)  ),
        S*(beta/A) * (fngcide(alpha_1*omega_1,alpha_1_C*params['theta_1'],Fung1)) * (fngcide(alpha_2*omega_2,alpha_2_C*theta_2,Fung2)) * (IR + PR)   - (Gamma(params['T_GS61'],params['T_GS87'],t)) * ER  - gamma * (fngcide(alpha_1*omega_1_L,alpha_1_C*params['theta_1_L'],Fung1))*(fngcide(alpha_2*omega_2_L,alpha_2_C*params['theta_2_L'],Fung2)) * ER,
        S*(beta/A) * (fngcide(alpha_1*omega_1,alpha_1_C*params['theta_1'],Fung1)) * (fngcide(        omega_2,          theta_2,Fung2)) * (IRS + PRS) - (Gamma(params['T_GS61'],params['T_GS87'],t)) * ERS - gamma * (fngcide(alpha_1*omega_1_L,alpha_1_C*params['theta_1_L'],Fung1))*(fngcide(        omega_2_L,          params['theta_2_L'],Fung2)) * ERS,
        S*(beta/A) * (fngcide(        omega_1,          params['theta_1'],Fung1)) * (fngcide(alpha_2*omega_2,alpha_2_C*theta_2,Fung2)) * (ISR + PSR) - (Gamma(params['T_GS61'],params['T_GS87'],t)) * ESR - gamma * (fngcide(        omega_1_L,          params['theta_1_L'],Fung1))*(fngcide(alpha_2*omega_2_L,alpha_2_C*params['theta_2_L'],Fung2)) * ESR,
        S*(beta/A) * (fngcide(        omega_1,          params['theta_1'],Fung1)) * (fngcide(        omega_2,          theta_2,Fung2)) * (IS + PS)   - (Gamma(params['T_GS61'],params['T_GS87'],t)) * ES  - gamma * (fngcide(        omega_1_L,          params['theta_1_L'],Fung1))*(fngcide(        omega_2_L,          params['theta_2_L'],Fung2)) * ES,
        gamma * (fngcide(alpha_1*omega_1_L,alpha_1_C*params['theta_1_L'],Fung1)) * (fngcide(alpha_2*omega_2_L,alpha_2_C*params['theta_2_L'],Fung2)) * ER   -  mu * IR,
        gamma * (fngcide(alpha_1*omega_1_L,alpha_1_C*params['theta_1_L'],Fung1)) * (fngcide(        omega_2_L,          params['theta_2_L'],Fung2)) * ERS  -  mu * IRS,
        gamma * (fngcide(        omega_1_L,          params['theta_1_L'],Fung1)) * (fngcide(alpha_2*omega_2_L,alpha_2_C*params['theta_2_L'],Fung2)) * ESR  -  mu * ISR,
        gamma * (fngcide(        omega_1_L,          params['theta_1_L'],Fung1)) * (fngcide(        omega_2_L,          params['theta_2_L'],Fung2)) * ES   -  mu * IS,
        mu * (IR + IRS + ISR + IS)   +  (Gamma(t)) * (S + ER + ERS + ESR + ES),
        -nu * PR,
        -nu * PRS,
        -nu * PSR,
        -nu * PS,
        -params['delta_1'] * Fung1,
        -params['delta_2'] * Fung2]
    return dydt

#----------------------------------------------------------------------------------------------
def loop_RF_SR_Y_F(dose_11, dose_12, dose_21, dose_22,res_prop_1,res_prop_2,numberofseasons):
    Failure_year = 0
    Y_df = calculator_d_free()
    # print(Y_df,'is disease free yield, new calculator')
    Y_vec   = np.linspace(0,1,numberofseasons)
    R_vec_1 = np.linspace(0,1,numberofseasons+1)
    R_vec_2 = np.linspace(0,1,numberofseasons+1)
    S_vec   = np.linspace(1,numberofseasons,numberofseasons)
    S_vec2  = np.linspace(0,numberofseasons,numberofseasons+1)
    i=1
    rpn1 = res_prop_1
    rpn2 = res_prop_2
    R_vec_1[0] = res_prop_1
    R_vec_2[0] = res_prop_2
    innoc = params['init_den']
    while i<=numberofseasons:
        ####Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, inte, solution, solutiont = calculator_ode_int(dose_11, dose_12, dose_21, dose_22,rpn1,rpn2,innoc,'integrated')
        Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, inte, solution, solutiont = calculator_ode(dose_11, dose_12, dose_21, dose_22,rpn1,rpn2,innoc,params['dt'],nstepz,'integrated')
        #print(solution.shape)
        Yield = 100 *(inte/Y_df)
        #print(Yield,'compare to old ',100*inte/Y_df_old)
        Y_vec[i-1] = Yield
        #print('Year ',i,' yield is ',Yield,'%')
        R_vec_1[i] = res_prop_new_1
        R_vec_2[i] = res_prop_new_2
        rpn1 = res_prop_new_1
        rpn2 = res_prop_new_2
        innoc = inoculum_integrated(solution,solutiont)
        # if i==1:
        #print('Year',i,'has innoc =',innoc,'. Compare to',params['init_den'])
#----------------------------------------------------------------------------------------------
        # print(params['innoc_frac']*(solution[-1,5] + solution[-1,6] + solution[-1,7] + solution[-1,8]),' is disease times frac')
        # print((solution[-1,5] + solution[-1,6] + solution[-1,7] + solution[-1,8]),' is disease')
        # print(solution[-1,5],' is IRR')
        # print(solution[-1,6],' is IRS')
        # print(solution[-1,7],' is ISR')
        # print(solution[-1,8],' is ISS')
        if Y_vec[i-1]<params['Yield_threshold'] and Failure_year==0:
            Failure_year = i
            #print('Failure year is ',i)
        i=i+1
    if min(Y_vec) >params['Yield_threshold']:
        print('Never fails, min yield is',min(Y_vec))
        #Failure_year=-10
    return R_vec_1, R_vec_2, Y_vec, S_vec, S_vec2, Failure_year

#----------------------------------------------------------------------------------------------
def loop_Y_SR_LY(res_prop_1,res_prop_2,numberofdoses,numberofseasons):
    Y_df = calculator_d_free()
    #print(Y_df,'is disease free yield')
    Y_v_m   = np.linspace(0,1,numberofdoses) 
    Y_v_lh  = np.linspace(0,1,numberofdoses) 
    Y_v_hl  = np.linspace(0,1,numberofdoses) 
    S_v1_m  = np.linspace(0,1,numberofdoses) 
    S_v1_lh = np.linspace(0,1,numberofdoses) 
    S_v1_hl = np.linspace(0,1,numberofdoses) 
    S_v2_m  = np.linspace(0,1,numberofdoses) 
    S_v2_lh = np.linspace(0,1,numberofdoses) 
    S_v2_hl = np.linspace(0,1,numberofdoses) 
    LY_v_m  = np.linspace(0,1,numberofdoses) 
    LY_v_lh = np.linspace(0,1,numberofdoses) 
    LY_v_hl = np.linspace(0,1,numberofdoses) 
    D_vec   = np.linspace(0,1,numberofdoses)
    j=0
    i=1
    innoc_m  = params['init_den']
    innoc_lh = params['init_den']
    innoc_hl = params['init_den']
    while i<=numberofdoses:
        ####Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode_int(0.5, 0.5*j, 0.5, 0.5*j,res_prop_1,res_prop_2,innoc_m,'integrated')
        Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode(0.5, 0.5*j, 0.5, 0.5*j,res_prop_1,res_prop_2,innoc_m,params['dt'],nstepz,'integrated')
        Yield_m     = 100 *(inte_m/Y_df)
        Y_v_m[i-1]  = Yield_m
        S_v1_m[i-1]  = Selection_1_m
        S_v2_m[i-1]  = Selection_2_m
#----------------------------------------------------------------------------------------------
        ####Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode_int(1, 0, 0, j,res_prop_1,res_prop_2,innoc_lh,'integrated')
        Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode(1, 0, 0, j,res_prop_1,res_prop_2,innoc_lh,params['dt'],nstepz,'integrated')
        Yield_lh    = 100 *(inte_lh/Y_df)
        Y_v_lh[i-1] = Yield_lh
        S_v1_lh[i-1] = Selection_1_lh
        S_v2_lh[i-1] = Selection_2_lh
#----------------------------------------------------------------------------------------------
        ####Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode_int(0, j, 1, 0,res_prop_1,res_prop_2,innoc_hl,'integrated')
        Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode(0, j, 1, 0,res_prop_1,res_prop_2,innoc_hl,params['dt'],nstepz,'integrated')
        Yield_hl     = 100 *(inte_hl/Y_df)
        Y_v_hl[i-1]  = Yield_hl
        S_v1_hl[i-1]  = Selection_1_hl
        S_v2_hl[i-1]  = Selection_2_hl
#----------------------------------------------------------------------------------------------
        R_vec1_m,  R_vec2_m,  Y_vec_m,  S_vec_m,  S_vec2_m,  F_y_m   = loop_RF_SR_Y_F(.5,.5*j,.5,.5*j,res_prop_1,res_prop_2,numberofseasons) #mixture
        R_vec1_hl, R_vec2_hl, Y_vec_hl, S_vec_hl, S_vec2_hl, F_y_hl  = loop_RF_SR_Y_F(0,1*j,1,0,      res_prop_1,res_prop_2,numberofseasons) #highlow
        R_vec1_lh, R_vec2_lh, Y_vec_lh, S_vec_lh, S_vec2_lh, F_y_lh  = loop_RF_SR_Y_F(1,0,0,1*j,      res_prop_1,res_prop_2,numberofseasons) #lowhigh
        
        LifeYield_m, LifeYield_hl, LifeYield_lh = lifetime_yield(Y_vec_m,Y_vec_hl,Y_vec_lh,F_y_m,F_y_hl,F_y_lh)

        LY_v_m[i-1]   = LifeYield_m       
        LY_v_lh[i-1]  = LifeYield_lh
        LY_v_hl[i-1]  = LifeYield_hl
#----------------------------------------------------------------------------------------------
        print([i,j])
        j=j+1/(numberofdoses-1)
        i=i+1
    return Y_v_m, Y_v_lh, Y_v_hl, S_v1_m, S_v1_lh, S_v1_hl, S_v2_m, S_v2_lh, S_v2_hl, LY_v_m, LY_v_lh, LY_v_hl, D_vec

#----------------------------------------------------------------------------------------------
def loop_Y_LY_matrix(res_prop_1,res_prop_2,numberofdoses,numberofseasons):
    Y_df = calculator_d_free()
    #print(Y_df,'is disease free yield')
    Y_v_m   = np.zeros((numberofdoses,numberofdoses))
    Y_v_lh  = np.zeros((numberofdoses,numberofdoses))
    Y_v_hl  = np.zeros((numberofdoses,numberofdoses))
    LY_v_m  = np.zeros((numberofdoses,numberofdoses))
    LY_v_lh = np.zeros((numberofdoses,numberofdoses))
    LY_v_hl = np.zeros((numberofdoses,numberofdoses))
    #D_vec   = np.linspace(0,1,numberofdoses)
    d_fn1=0
    no_dose_1=1
    d_fn2 =0
    no_dose_2=1
    innoc_m  = params['init_den']
    innoc_lh = params['init_den']
    innoc_hl = params['init_den']
    while no_dose_1<=numberofdoses:
        no_dose_2=1
        d_fn2=0
        while no_dose_2<=numberofdoses:
            ####Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode_int(0.5*d_fn1, 0.5*d_fn2, 0.5*d_fn1, 0.5*d_fn2,res_prop_1,res_prop_2,innoc_m,'integrated')
            Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode(0.5*d_fn1, 0.5*d_fn2, 0.5*d_fn1, 0.5*d_fn2,res_prop_1,res_prop_2,innoc_m,params['dt'],nstepz,'integrated')
            Yield_m     = 100 *(inte_m/Y_df)
            Y_v_m[no_dose_1-1,no_dose_2-1]  = Yield_m
#----------------------------------------------------------------------------------------------
            ####Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode_int(d_fn1, 0, 0, d_fn2,res_prop_1,res_prop_2,innoc_lh,'integrated')
            Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode(d_fn1, 0, 0, d_fn2,res_prop_1,res_prop_2,innoc_lh,params['dt'],nstepz,'integrated')
            Yield_lh    = 100 *(inte_lh/Y_df)
            Y_v_lh[no_dose_1-1,no_dose_2-1] = Yield_lh
#----------------------------------------------------------------------------------------------
            ####Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode_int(0, d_fn2, d_fn1, 0,res_prop_1,res_prop_2,innoc_hl,'integrated')
            Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode(0, d_fn2, d_fn1, 0,res_prop_1,res_prop_2,innoc_hl,params['dt'],nstepz,'integrated')
            Yield_hl     = 100 *(inte_hl/Y_df)
            Y_v_hl[no_dose_1-1,no_dose_2-1]  = Yield_hl
#----------------------------------------------------------------------------------------------
            R_vec1_m,  R_vec2_m,  Y_vec_m,  S_vec_m,  S_vec2_m,  F_y_m   = loop_RF_SR_Y_F(.5*d_fn1,.5*d_fn2, .5*d_fn1,.5*d_fn2,res_prop_1,res_prop_2,numberofseasons) #mixture
            R_vec1_hl, R_vec2_hl, Y_vec_hl, S_vec_hl, S_vec2_hl, F_y_hl  = loop_RF_SR_Y_F(    0   ,   d_fn2,    d_fn1,       0,res_prop_1,res_prop_2,numberofseasons) #highlow
            R_vec1_lh, R_vec2_lh, Y_vec_lh, S_vec_lh, S_vec2_lh, F_y_lh  = loop_RF_SR_Y_F(   d_fn1,       0,        0,   d_fn2,res_prop_1,res_prop_2,numberofseasons) #lowhigh
            
            LifeYield_m, LifeYield_hl, LifeYield_lh = lifetime_yield(Y_vec_m,Y_vec_hl,Y_vec_lh,F_y_m,F_y_hl,F_y_lh)

            LY_v_m[no_dose_1-1,no_dose_2-1]   = LifeYield_m
            LY_v_lh[no_dose_1-1,no_dose_2-1]  = LifeYield_lh
            LY_v_hl[no_dose_1-1,no_dose_2-1]  = LifeYield_hl
#----------------------------------------------------------------------------------------------
            print([no_dose_1,no_dose_2,d_fn1,d_fn2])
            d_fn2=d_fn2+1/(numberofdoses-1)
            no_dose_2=no_dose_2+1
        d_fn1=d_fn1+1/(numberofdoses-1)
        no_dose_1=no_dose_1+1
    return Y_v_m, Y_v_lh, Y_v_hl, LY_v_m, LY_v_lh, LY_v_hl

#----------------------------------------------------------------------------------------------
def loop_TY_matrix(res_prop_1,res_prop_2,numberofdoses,numberofseasons):
    Y_df = calculator_d_free()
    #print(Y_df,'is disease free yield')
    # Y_v_m   = np.zeros((numberofdoses,numberofdoses))
    # Y_v_lh  = np.zeros((numberofdoses,numberofdoses))
    # Y_v_hl  = np.zeros((numberofdoses,numberofdoses))
    TY_v_m  = np.zeros((numberofdoses,numberofdoses))
    TY_v_lh = np.zeros((numberofdoses,numberofdoses))
    TY_v_hl = np.zeros((numberofdoses,numberofdoses))
    #D_vec   = np.linspace(0,1,numberofdoses)
    d_fn1=0
    no_dose_1=1
    d_fn2 =0
    no_dose_2=1
    while no_dose_1<=numberofdoses:
        no_dose_2=1
        d_fn2=0
        while no_dose_2<=numberofdoses:  ## there was commented stuff here from copy paste about yield as in loop_Y_LY_matrix
            R_vec1_m,  R_vec2_m,  Y_vec_m,  S_vec_m,  S_vec2_m,  F_y_m   = loop_RF_SR_Y_F(.5*d_fn1,.5*d_fn2, .5*d_fn1,.5*d_fn2,res_prop_1,res_prop_2,numberofseasons) #mixture
            R_vec1_hl, R_vec2_hl, Y_vec_hl, S_vec_hl, S_vec2_hl, F_y_hl  = loop_RF_SR_Y_F(    0   ,   d_fn2,    d_fn1,       0,res_prop_1,res_prop_2,numberofseasons) #highlow
            R_vec1_lh, R_vec2_lh, Y_vec_lh, S_vec_lh, S_vec2_lh, F_y_lh  = loop_RF_SR_Y_F(   d_fn1,       0,        0,   d_fn2,res_prop_1,res_prop_2,numberofseasons) #lowhigh
            
            TotalYield_m, TotalYield_hl, TotalYield_lh = Total_yield(Y_vec_m,Y_vec_hl,Y_vec_lh,numberofseasons)

            TY_v_m[no_dose_1-1,no_dose_2-1]   = TotalYield_m
            TY_v_lh[no_dose_1-1,no_dose_2-1]  = TotalYield_lh
            TY_v_hl[no_dose_1-1,no_dose_2-1]  = TotalYield_hl
#----------------------------------------------------------------------------------------------
            print([no_dose_1,no_dose_2,d_fn1,d_fn2])
            d_fn2=d_fn2+1/(numberofdoses-1)
            no_dose_2=no_dose_2+1
        d_fn1=d_fn1+1/(numberofdoses-1)
        no_dose_1=no_dose_1+1
    return TY_v_m, TY_v_lh, TY_v_hl

#----------------------------------------------------------------------------------------------
def loop_SR_RF_matrix(res_prop_1,res_prop_2,numberofdoses):
    # Y_df = calculator_d_free()
    #D_vec   = np.linspace(0,1,numberofdoses)
    SR1_m  = np.zeros((numberofdoses,numberofdoses))
    SR1_lh = np.zeros((numberofdoses,numberofdoses))
    SR1_hl = np.zeros((numberofdoses,numberofdoses))
    SR2_m  = np.zeros((numberofdoses,numberofdoses))
    SR2_lh = np.zeros((numberofdoses,numberofdoses))
    SR2_hl = np.zeros((numberofdoses,numberofdoses))
    RF1_m  = np.zeros((numberofdoses,numberofdoses))
    RF1_lh = np.zeros((numberofdoses,numberofdoses))
    RF1_hl = np.zeros((numberofdoses,numberofdoses))
    RF2_m  = np.zeros((numberofdoses,numberofdoses))
    RF2_lh = np.zeros((numberofdoses,numberofdoses))
    RF2_hl = np.zeros((numberofdoses,numberofdoses))
    d_fn1=0
    no_dose_1=1
    d_fn2 =0
    no_dose_2=1
    innoc_m  = params['init_den']
    innoc_lh = params['init_den']
    innoc_hl = params['init_den']
    while no_dose_1<=numberofdoses:
        no_dose_2=1
        d_fn2=0
        while no_dose_2<=numberofdoses:
            #### Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode_int(0.5*d_fn1, 0.5*d_fn2, 0.5*d_fn1, 0.5*d_fn2,res_prop_1,res_prop_2,innoc_m,'integrated')
            Selection_1_m, Selection_2_m, res_prop_new_1_m, res_prop_new_2_m, inte_m, solution_m, solutiont_m = calculator_ode(0.5*d_fn1, 0.5*d_fn2, 0.5*d_fn1, 0.5*d_fn2,res_prop_1,res_prop_2,innoc_m,params['dt'],nstepz,'integrated')
#----------------------------------------------------------------------------------------------
            #### Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode_int(d_fn1, 0, 0, d_fn2,res_prop_1,res_prop_2,innoc_lh,'integrated')
            Selection_1_lh, Selection_2_lh, res_prop_new_1_lh, res_prop_new_2_lh, inte_lh, solution_lh, solutiont_lh = calculator_ode(d_fn1, 0, 0, d_fn2,res_prop_1,res_prop_2,innoc_lh,params['dt'],nstepz,'integrated')
#----------------------------------------------------------------------------------------------
            #### Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode_int(0, d_fn2, d_fn1, 0,res_prop_1,res_prop_2,innoc_hl,'integrated')
            Selection_1_hl, Selection_2_hl, res_prop_new_1_hl, res_prop_new_2_hl, inte_hl, solution_hl, solutiont_hl = calculator_ode(0, d_fn2, d_fn1, 0,res_prop_1,res_prop_2,innoc_hl,params['dt'],nstepz,'integrated')
#----------------------------------------------------------------------------------------------
            SR1_m[no_dose_1-1,no_dose_2-1]   = Selection_1_m
            SR1_lh[no_dose_1-1,no_dose_2-1]  = Selection_1_lh
            SR1_hl[no_dose_1-1,no_dose_2-1]  = Selection_1_hl
            SR2_m[no_dose_1-1,no_dose_2-1]   = Selection_2_m
            SR2_lh[no_dose_1-1,no_dose_2-1]  = Selection_2_lh
            SR2_hl[no_dose_1-1,no_dose_2-1]  = Selection_2_hl
#----------------------------------------------------------------------------------------------
            RF1_m[no_dose_1-1,no_dose_2-1]   = res_prop_new_1_m
            RF1_lh[no_dose_1-1,no_dose_2-1]  = res_prop_new_1_lh
            RF1_hl[no_dose_1-1,no_dose_2-1]  = res_prop_new_1_hl
            RF2_m[no_dose_1-1,no_dose_2-1]   = res_prop_new_2_m
            RF2_lh[no_dose_1-1,no_dose_2-1]  = res_prop_new_2_lh
            RF2_hl[no_dose_1-1,no_dose_2-1]  = res_prop_new_2_hl
#----------------------------------------------------------------------------------------------            
            print([no_dose_1,no_dose_2,d_fn1,d_fn2])
            d_fn2=d_fn2+1/(numberofdoses-1)
            no_dose_2=no_dose_2+1
        d_fn1=d_fn1+1/(numberofdoses-1)
        no_dose_1=no_dose_1+1
    return SR1_m, SR1_lh, SR1_hl, SR2_m, SR2_lh, SR2_hl, RF1_m, RF1_lh, RF1_hl, RF2_m, RF2_lh, RF2_hl





#----------------------------------------------------------------------------------------------
def loop_RF_SR_Y_F_var_dose(dose_11_vec,dose_12_vec,dose_21_vec,dose_22_vec,res_prop_1,res_prop_2):
    Failure_year = 0
    Y_df = calculator_d_free()
    Y_vec   = np.linspace(0,1,len(dose_11_vec))
    R_vec_1 = np.linspace(0,1,len(dose_11_vec)+1)
    R_vec_2 = np.linspace(0,1,len(dose_11_vec)+1)
    I_vec = np.linspace(0,1,len(dose_11_vec)+1)
    S_vec   = np.linspace(1,len(dose_11_vec),len(dose_11_vec))
    S_vec2  = np.linspace(0,len(dose_11_vec),len(dose_11_vec)+1)
    rpn1 = res_prop_1
    rpn2 = res_prop_2
    R_vec_1[0] = res_prop_1
    R_vec_2[0] = res_prop_2
    I_vec[0]   = params['init_den']
    innoc      = params['init_den']
    for i in range(len(dose_11_vec)):
        ####Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, inte, solution, solutiont = calculator_ode_int(dose_11_vec[i],dose_12_vec[i],dose_21_vec[i],dose_22_vec[i],params['T_GS87'],rpn1,rpn2,innoc,'integrated')
        Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, inte, solution, solutiont = calculator_ode(dose_11_vec[i],dose_12_vec[i],dose_21_vec[i],dose_22_vec[i],rpn1,rpn2,innoc,params['dt'],nstepz,'integrated')
        #print(solution.shape)
        Yield = 100 *(inte/Y_df)
        Y_vec[i] = Yield
        #print('Year ',i,' yield is ',Yield,'%')
        R_vec_1[i+1] = res_prop_new_1
        R_vec_2[i+1] = res_prop_new_2
        rpn1 = res_prop_new_1
        rpn2 = res_prop_new_2
        innoc = inoculum_integrated(solution,solutiont)
        I_vec[i+1] = innoc
        # if i==1:
        #print('Year',i,'has innoc =',innoc,'. Compare to',params['init_den'])
#----------------------------------------------------------------------------------------------
        if Y_vec[i]<params['Yield_threshold'] and Failure_year==0:
            Failure_year = i+1
            #print('Failure year is ',i)
    if min(Y_vec) >params['Yield_threshold']:
        print('Never fails, min yield is',min(Y_vec))
        #Failure_year=-10
    return R_vec_1, R_vec_2, Y_vec, I_vec, S_vec, S_vec2, Failure_year
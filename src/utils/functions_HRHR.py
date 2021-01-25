import numpy as np
from math import exp, ceil, floor
from scipy.integrate import odeint, simps, ode
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import pickle
import json
import os
import pdb

from .parameters_HRHR import params, params_dict
#----------------------------------------------------------------------------------------------
def growth(A,t):
    if t>=params.T_emerge:
        grw = params.r*(params.k-A)
        return grw
    else:
        grw=0    
    return grw


#----------------------------------------------------------------------------------------------
def fngcide(omega,theta,fungicide_dose):
    effect = 1 - omega*(1 - exp(-theta*fungicide_dose))
    return effect


#----------------------------------------------------------------------------------------------
def Gamma(t):
    if t>=params.T_GS61:
        Gmm = 0.005*((t-params.T_GS61)/(params.T_GS87-params.T_GS61)) + 0.1*exp(-0.02*(params.T_GS87-t))
        return Gmm
    else:
        Gmm = 0
        return Gmm


#----------------------------------------------------------------------------------------------
# for calculator_ode
def solv_ode(t,y):
    S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2 = y
    A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R
    dydt = [growth(A,t) - (Gamma(t))*S -  S * (params.beta/A) * (  
            (fngcide(params.alpha_1*params.omega_1,params.alpha_1_C*params.theta_1,Fung1)) * (fngcide(params.alpha_2*params.omega_2,params.alpha_2_C*params.theta_2,Fung2)) * (IR + PR)
            + (fngcide(params.alpha_1*params.omega_1,params.alpha_1_C*params.theta_1,Fung1)) * (fngcide(                  params.omega_2,                    params.theta_2,Fung2)) * (IRS+PRS)
            + (fngcide(                  params.omega_1,                    params.theta_1,Fung1)) * (fngcide(params.alpha_2*params.omega_2,params.alpha_2_C*params.theta_2,Fung2)) * (ISR + PSR)
            + (fngcide(                  params.omega_1,                    params.theta_1,Fung1)) * (fngcide(                  params.omega_2,                    params.theta_2,Fung2)) * (IS + PS)  ),
        S*(params.beta/A) * (fngcide(params.alpha_1*params.omega_1,params.alpha_1_C*params.theta_1,Fung1)) * (fngcide(params.alpha_2*params.omega_2,params.alpha_2_C*params.theta_2,Fung2)) * (IR + PR)   - (Gamma(t)) * ER  - params.gamma * (fngcide(params.alpha_1*params.omega_1_L,params.alpha_1_C*params.theta_1_L,Fung1))*(fngcide(params.alpha_2*params.omega_2_L,params.alpha_2_C*params.theta_2_L,Fung2)) * ER,
        S*(params.beta/A) * (fngcide(params.alpha_1*params.omega_1,params.alpha_1_C*params.theta_1,Fung1)) * (fngcide(                  params.omega_2,                    params.theta_2,Fung2)) * (IRS + PRS) - (Gamma(t)) * ERS - params.gamma * (fngcide(params.alpha_1*params.omega_1_L,params.alpha_1_C*params.theta_1_L,Fung1))*(fngcide(                  params.omega_2_L,                    params.theta_2_L,Fung2)) * ERS,
        S*(params.beta/A) * (fngcide(                  params.omega_1,                    params.theta_1,Fung1)) * (fngcide(params.alpha_2*params.omega_2,params.alpha_2_C*params.theta_2,Fung2)) * (ISR + PSR) - (Gamma(t)) * ESR - params.gamma * (fngcide(                  params.omega_1_L,                    params.theta_1_L,Fung1))*(fngcide(params.alpha_2*params.omega_2_L,params.alpha_2_C*params.theta_2_L,Fung2)) * ESR,
        S*(params.beta/A) * (fngcide(                  params.omega_1,                    params.theta_1,Fung1)) * (fngcide(                  params.omega_2,                    params.theta_2,Fung2)) * (IS + PS)   - (Gamma(t)) * ES  - params.gamma * (fngcide(                  params.omega_1_L,                    params.theta_1_L,Fung1))*(fngcide(                  params.omega_2_L,                    params.theta_2_L,Fung2)) * ES,
        params.gamma * (fngcide(params.alpha_1*params.omega_1_L,params.alpha_1_C*params.theta_1_L,Fung1)) * (fngcide(params.alpha_2*params.omega_2_L,params.alpha_2_C*params.theta_2_L,Fung2)) * ER   -  params.mu * IR,
        params.gamma * (fngcide(params.alpha_1*params.omega_1_L,params.alpha_1_C*params.theta_1_L,Fung1)) * (fngcide(                  params.omega_2_L,                    params.theta_2_L,Fung2)) * ERS  -  params.mu * IRS,
        params.gamma * (fngcide(                  params.omega_1_L,                    params.theta_1_L,Fung1)) * (fngcide(params.alpha_2*params.omega_2_L,params.alpha_2_C*params.theta_2_L,Fung2)) * ESR  -  params.mu * ISR,
        params.gamma * (fngcide(                  params.omega_1_L,                    params.theta_1_L,Fung1)) * (fngcide(                  params.omega_2_L,                    params.theta_2_L,Fung2)) * ES   -  params.mu * IS,
        params.mu * (IR + IRS + ISR + IS)   +  (Gamma(t)) * (S + ER + ERS + ESR + ES),
        -params.nu * PR,
        -params.nu * PRS,
        -params.nu * PSR,
        -params.nu * PS,
        -params.delta_1 * Fung1,
        -params.delta_2 * Fung2]
    return dydt


#----------------------------------------------------------------------------------------------
def calculator_d_free(): # calculates disease free yield
    y0   = [params.S_0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    sol  = ode(solv_ode,jac=None).set_integrator('dopri5',nsteps= params.nstepz)#,rtol=10**(-10),atol=10**(-20))
    t0, t1, t2 = (params.T_emerge,params.T_GS61,params.T_GS87)
    n1= 1 + (t1-t0)/params.dt
    n2= 1 + (t2-t1)/params.dt
    c1 = ceil(n1-0.5)
    c2 = ceil(n2-0.5)
    tim1 = np.linspace(t0,t1,c1)
    tim2 = np.linspace(t1,t2,c2)
    yy1  = np.zeros((params.no_variables,len(tim1)))
    yy2  = np.zeros((params.no_variables,len(tim2)))
    i1=0
    i2=0
#----------------------------------------------------------------------------------------------
    sol.set_initial_value(y0,t0)
    for t in tim1[1:]:
        if sol.successful():
            yy1[:,i1] = sol.y
            i1=i1+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
    #----------------------------------------------------------------------------------------------
    y1 = sol.y
    sol.set_initial_value(y1,t1)
    for t in tim2[1:]:
        if sol.successful():
            yy2[:,i2] = sol.y
            i2=i2+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
# #----------------------------------------------------------------------------------------------
    yy2[:,i2] = sol.y
# #----------------------------------------------------------------------------------------------  
    int_df = simps(yy2[0,:],tim2)
#----------------------------------------------------------------------------------------------
    return int_df


#----------------------------------------------------------------------------------------------
def lifetime_yield(Y_vec,F_y):
    i = 1
    j = 0
    while i < F_y:
        j = j+Y_vec[i-1]
        i = i+1
    j = j/100
    return j

#----------------------------------------------------------------------------------------------
def Total_yield(Y_vec,numberofseasons):
    i = 1
    j = 0
    while i<=numberofseasons:
        j = j+Y_vec[i-1]
        i = i+1
    j = j/100
    return j

#----------------------------------------------------------------------------------------------
def inoculum_value(solution,solutiont=None,FV_or_int=None):
    if FV_or_int is None or FV_or_int == 'final_value':
        inn = params.init_den*(params.den_frac + (1-params.den_frac)*params.innoc_frac*(solution[-1,params.IR_ind]+solution[-1,params.IRS_ind]+solution[-1,params.ISR_ind]+solution[-1,params.IS_ind]))
    if FV_or_int=='integrated':
        disease = simps(solution[:,params.IR_ind]+solution[:,params.IRS_ind]+solution[:,params.ISR_ind]+solution[:,params.IS_ind],solutiont)
        inn = params.init_den*(params.den_frac + (1-params.den_frac)*params.innoc_frac_integral * disease)
    return inn

#----------------------------------------------------------------------------------------------
def resist_prop_calculator(solution,solutiont=None,FV_or_int=None):
    if FV_or_int is None or FV_or_int == 'final_value':
        disease = solution[-1,params.IR_ind]+solution[-1,params.IRS_ind]+solution[-1,params.ISR_ind]+solution[-1,params.IS_ind]
        Res_disease_1 = solution[-1,params.IR_ind]+solution[-1,params.IRS_ind]
        Res_disease_2 = solution[-1,params.IR_ind]+solution[-1,params.ISR_ind]
        RP1 = Res_disease_1/disease # phi_1
        RP2 = Res_disease_2/disease # phi_2
        p_rr = solution[-1,params.IR_ind]/disease
        p_rs = solution[-1,params.IRS_ind]/disease
        p_sr = solution[-1,params.ISR_ind]/disease
        p_ss = solution[-1,params.IS_ind]/disease
    if FV_or_int == 'integrated':
        disease = simps(solution[:,params.IR_ind]+solution[:,params.IRS_ind]+solution[:,params.ISR_ind]+solution[:,params.IS_ind],solutiont)
        Res_disease_1 = simps(solution[:,params.IR_ind]+solution[:,params.IRS_ind],solutiont)
        Res_disease_2 = simps(solution[:,params.IR_ind]+solution[:,params.ISR_ind],solutiont)
        RP1 = Res_disease_1/disease # phi_1
        RP2 = Res_disease_2/disease # phi_2
        p_rr = simps(solution[:,params.IR_ind],solutiont)/disease
        p_rs = simps(solution[:,params.IRS_ind],solutiont)/disease
        p_sr = simps(solution[:,params.ISR_ind],solutiont)/disease
        p_ss = simps(solution[:,params.IS_ind],solutiont)/disease
    return RP1, RP2, p_rr, p_rs, p_sr, p_ss

#----------------------------------------------------------------------------------------------
def interpolate(array,input_number,output_number,kind='linear'):
    Interpolate = interp2d(np.linspace(0,1,input_number),np.linspace(0,1,input_number),np.transpose(array),kind=kind)

    Int = np.zeros((output_number,output_number))
    for i in range(output_number):
            for j in range(output_number):
                    Int[i,j] = Interpolate(i/(output_number-1),j/(output_number-1))
    return Interpolate, Int

#----------------------------------------------------------------------------------------------
def primary_calculator(res_prop_1,res_prop_2,p_rr=None,p_rs=None,p_sr=None,p_ss=None):
    if p_rr is None and res_prop_1 is not None:
        P_RR0 =    res_prop_1*   res_prop_2  #used to have innoc  
        P_RS0 =    res_prop_1*(1-res_prop_2)#used to have innoc
        P_SR0 = (1-res_prop_1)*   res_prop_2 #used to have innoc   
        P_SS0 = (1-res_prop_1)*(1-res_prop_2)#used to have innoc
    else:   
        P_RR0 = ( params.sex_prop*   res_prop_1*   res_prop_2   + (1-params.sex_prop)*p_rr) # used to have innoc
        P_RS0 = ( params.sex_prop*   res_prop_1*(1-res_prop_2)  + (1-params.sex_prop)*p_rs) # used to have innoc
        P_SR0 = ( params.sex_prop*(1-res_prop_1)*  res_prop_2   + (1-params.sex_prop)*p_sr) # used to have innoc
        P_SS0 = ( params.sex_prop*(1-res_prop_1)*(1-res_prop_2) + (1-params.sex_prop)*p_ss) # used to have innoc
    return P_RR0 , P_RS0 , P_SR0, P_SS0




#----------------------------------------------------------------------------------------------
def calculator_ode(dose_11,dose_12,dose_21,dose_22, P_RR0, P_RS0, P_SR0, P_SS0):
    res_prop_1 = P_RR0 + P_RS0
    res_prop_2 = P_RR0 + P_SR0
    ##
    y0   = [params.S_0, 0,  0,  0, 0, 0,  0,  0, 0,0, P_RR0,  P_RS0, P_SR0, P_SS0,    0,   0]
    sol  = ode(solv_ode,jac=None).set_integrator('dopri5',nsteps= params.nstepz)#,rtol=0.3*10**(-13),atol=10**(-25))#,y0,t)#,hmin =0.01)#,hmax=0.0025)#.set_integrator('vode', with_jacobian=False)#  method='bdf', args = (r,k,omega_L,theta_L,omega_H,theta_H,T_GS61,params.T_GS87,beta,gamma,mu,nu,delta_H,delta_L))#,atol=10**(-20))
    t0, t1, t2, t3, t4 = (params.T_emerge,params.T_GS32,params.T_GS39,params.T_GS61,params.T_GS87)
    n1 = 1 + (t1-t0)/params.dt
    n2 = 1 + (t2-t1)/params.dt
    n3 = 1 + (t3-t2)/params.dt
    n4 = 1 + (t4-t3)/params.dt
    c1 = ceil(n1-0.5)
    c2 = ceil(n2-0.5)
    c3 = ceil(n3-0.5)
    c4 = ceil(n4-0.5)
    tim1 = np.linspace(t0,t1,c1)
    tim2 = np.linspace(t1,t2,c2)
    tim3 = np.linspace(t2,t3,c3)
    tim4 = np.linspace(t3,t4,c4)
    yy1  = np.zeros((params.no_variables,len(tim1)))
    yy2  = np.zeros((params.no_variables,len(tim2)))
    yy3  = np.zeros((params.no_variables,len(tim3)))
    yy4  = np.zeros((params.no_variables,len(tim4)))
    i1=0
    i2=0
    i3=0
    i4=0
#----------------------------------------------------------------------------------------------
    sol.set_initial_value(y0,t0)
    for t in tim1[1:]:
        if sol.successful():
            yy1[:,i1] = sol.y
            i1=i1+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
    #----------------------------------------------------------------------------------------------
    y1 = sol.y
    y1[params.Fung1_ind]=y1[params.Fung1_ind] + dose_11
    y1[params.Fung2_ind]=y1[params.Fung2_ind] + dose_12
    sol.set_initial_value(y1,t1)
    for t in tim2[1:]:
        if sol.successful():
            yy2[:,i2] = sol.y
            i2=i2+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
#----------------------------------------------------------------------------------------------
    y2 = sol.y #yy2[:,-1]
    y2[params.Fung1_ind]=y2[params.Fung1_ind] + dose_21
    y2[params.Fung2_ind]=y2[params.Fung2_ind] + dose_22
    sol.set_initial_value(y2,t2)
    for t in tim3[1:]:
        if sol.successful():
            yy3[:,i3] = sol.y
            i3=i3+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
#----------------------------------------------------------------------------------------------
    y3 = sol.y #yy3[:,-1]
    sol.set_initial_value(y3,t3)
    for t in tim4[1:]:
        if sol.successful():
            yy4[:,i4] = sol.y
            i4=i4+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
# #----------------------------------------------------------------------------------------------
    yy4[:,i4] = sol.y
# #----------------------------------------------------------------------------------------------
    tg     = tim1[:-1]
    t1g    = tim2[:-1]
    t2g    = tim3[:-1]
    solg   = yy1[:,:-1]
    sol1g  = yy2[:,:-1]
    sol2g  = yy3[:,:-1]
# #----------------------------------------------------------------------------------------------
    solutiont = np.concatenate((tg,t1g,t2g,tim4))
    solutionTranspose  = np.concatenate((solg,sol1g,sol2g,yy4),axis=1)
    solution  = np.transpose(solutionTranspose)
# #----------------------------------------------------------------------------------------------  
    res_prop_new_1, res_prop_new_2, p_rr, p_rs, p_sr, p_ss = resist_prop_calculator(solution,solutiont,FV_or_int=params.res_prop_calc_method) # gives either integral or final value
    innoc = inoculum_value(solution,solutiont,FV_or_int=params.res_prop_calc_method)
    Selection_1 = 1 # defaults
    Selection_2 = 1
    if res_prop_1 > 0:
        Selection_1    = res_prop_new_1/(res_prop_1/params.init_den)
    if res_prop_2 > 0:
        Selection_2    = res_prop_new_2/(res_prop_2/params.init_den)
# #----------------------------------------------------------------------------------------------  
    inte = simps(yy4[params.S_ind,:]+yy4[params.ER_ind,:]+yy4[params.ERS_ind,:]+yy4[params.ESR_ind,:]+yy4[params.ES_ind,:],tim4)
    return Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, innoc, p_rr, p_rs, p_sr, p_ss, inte, solution, solutiont

#----------------------------------------------------------------------------------------------  
def master_loop_one_tactic(dose_11_vec,dose_12_vec,dose_21_vec,dose_22_vec,res_prop_1=None,res_prop_2=None,yield_stopper=0,p_rr=None,p_rs=None,p_sr=None,p_ss=None,YDF=None,within_season_before=True):
    Failure_year = 0
    if YDF is None:
        YDF = calculator_d_free()
    ##
    Yield_vec = np.zeros(len(dose_11_vec))
    Selection_vec_1,Selection_vec_2,Innoc_vec, Res_vec_1, Res_vec_2, PRR, PRS, PSR, PSS = [np.zeros(len(dose_11_vec)+1) for i in range(9)]
    Sol_array = np.zeros((ceil((params.T_GS87-params.T_emerge)/params.dt),params.no_variables,len(dose_11_vec))) # array that has solution for each state variable for each year.
    t_vec = np.zeros(ceil((params.T_GS87-params.T_emerge)/params.dt))
    ##
    if p_rr is None:
        PRR[0]       = res_prop_1*res_prop_2
        PRS[0]       = res_prop_1*(1-res_prop_2)
        PSR[0]       = (1-res_prop_1)*res_prop_2
        PSS[0]       = (1-res_prop_1)*(1-res_prop_2)
        Res_vec_1[0] = res_prop_1
        Res_vec_2[0] = res_prop_2
    else:
        PRR[0]       = p_rr
        PRS[0]       = p_rs
        PSR[0]       = p_sr
        PSS[0]       = p_ss
    Innoc_vec[0]   = params.init_den
    ##
    for i in range(len(dose_11_vec)):
        if i==0 or (i>0 and Yield_vec[i-1]>yield_stopper): # stops the solver after we drop below threshold
            innoc_in = Innoc_vec[i]
            P_RR0 = PRR[i]
            P_RS0 = PRS[i]
            P_SR0 = PSR[i]
            P_SS0 = PSS[i]
            if within_season_before:
                res_prop_1_in = P_RR0 + P_RS0
                res_prop_2_in = P_RR0 + P_SR0
                P_RR0 , P_RS0 , P_SR0, P_SS0 = primary_calculator(res_prop_1_in,res_prop_2_in,P_RR0,P_RS0,P_SR0,P_SS0)
            Selection_1, Selection_2, res_prop_new_1, res_prop_new_2, innoc_out, p_rr_out, p_rs_out, p_sr_out, p_ss_out, inte, solution, solutiont = calculator_ode(dose_11=dose_11_vec[i],dose_12=dose_12_vec[i],dose_21=dose_21_vec[i],dose_22=dose_22_vec[i],P_RR0=P_RR0*innoc_in,P_RS0=P_RS0*innoc_in,P_SR0=P_SR0*innoc_in,P_SS0=P_SS0*innoc_in)
            Yield = 100*(inte/YDF)
            Yield_vec[i] = Yield
            if not within_season_before:
                res_prop_1_in = p_rr_out + p_rs_out
                res_prop_2_in = p_rr_out + p_sr_out
                p_rr_out , p_rs_out , p_sr_out, p_ss_out = primary_calculator(res_prop_1_in,res_prop_2_in,p_rr_out,p_rs_out,p_sr_out,p_ss_out)
            Res_vec_1[i+1] = res_prop_new_1
            Res_vec_2[i+1] = res_prop_new_2
            PRR[i+1]       = p_rr_out
            PRS[i+1]       = p_rs_out
            PSR[i+1]       = p_sr_out
            PSS[i+1]       = p_ss_out
            Innoc_vec[i+1] = innoc_out
            Selection_vec_1[i+1] = Selection_1
            Selection_vec_2[i+1] = Selection_2
            Sol_array[:,:,i] = solution
            t_vec = solutiont
#----------------------------------------------------------------------------------------------
            if Yield_vec[i]<params.Yield_threshold and Yield_vec[0]>params.Yield_threshold and Failure_year==0:
                Failure_year = i+1
    if min(Yield_vec) >params.Yield_threshold:
        Failure_year = -1
    dictionary = {'Res_vec_1': Res_vec_1, 'Res_vec_2': Res_vec_2, 'PRR': PRR, 'PRS': PRS, 'PSR': PSR, 'PSS': PSS, 'Yield_vec': Yield_vec, 'Innoc_vec': Innoc_vec, 'Selection_vec_1': Selection_vec_1, 'Selection_vec_2': Selection_vec_2, 'Failure_year': Failure_year, 'Sol_array': Sol_array, 't_vec': t_vec}
    return dictionary


#----------------------------------------------------------------------------------------------
def master_loop_grid_of_tactics(number_of_doses,number_of_seasons,res_freq_1=None,res_freq_2=None,yield_stopper=0,strategy=params.strategy,p_rr=None,p_rs=None,p_sr=None,p_ss=None,YDF=None,within_season_before=1):
    if YDF is None:
        YDF = calculator_d_free()
    LTY, TY, FY = [np.zeros((number_of_doses,number_of_doses)) for i in range(3)]
    ##
    Yield = np.zeros((number_of_doses,number_of_doses,number_of_seasons))
    Res_array_1,Res_array_2,PRR_array,PRS_array,PSR_array,PSS_array,Selection_array_1,Selection_array_2, Innoc_array  = [np.zeros((number_of_doses,number_of_doses,number_of_seasons+1)) for i in range(9)]
    ##
    for i in range(number_of_doses):
        for j in range(number_of_doses):
            if strategy == 'mix':
                dose_11_vec = 0.5*(i/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_12_vec = 0.5*(j/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_21_vec = 0.5*(i/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_22_vec = 0.5*(j/(number_of_doses-1))*np.ones(number_of_seasons)
            if strategy == 'alt_12':
                dose_11_vec = (i/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_12_vec = np.zeros(number_of_seasons)
                dose_21_vec = np.zeros(number_of_seasons)
                dose_22_vec = (j/(number_of_doses-1))*np.ones(number_of_seasons)
            if strategy == 'alt_21':
                dose_11_vec = np.zeros(number_of_seasons)
                dose_12_vec = (j/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_21_vec = (i/(number_of_doses-1))*np.ones(number_of_seasons)
                dose_22_vec = np.zeros(number_of_seasons)
            one_tact_output =  master_loop_one_tactic(dose_11_vec=dose_11_vec,dose_12_vec=dose_12_vec,dose_21_vec=dose_21_vec,dose_22_vec=dose_22_vec,res_prop_1= res_freq_1,res_prop_2= res_freq_2,yield_stopper= yield_stopper,p_rr=p_rr,p_rs=p_rs,p_sr=p_sr,p_ss=p_ss,YDF=YDF,within_season_before=within_season_before)
            LTY[i,j]  = lifetime_yield(one_tact_output['Yield_vec'],one_tact_output['Failure_year'])
            TY[i,j]   = Total_yield(one_tact_output['Yield_vec'],number_of_seasons)
            FY[i,j]   = one_tact_output['Failure_year']

            attr = {'Yield_vec': Yield,'Res_vec_1':Res_array_1,'Res_vec_2':Res_array_2,'PRR':PRR_array,'PRS':PRS_array,'PSR':PSR_array,'PSS':PSS_array,'Selection_vec_1':Selection_array_1,'Selection_vec_2':Selection_array_2,'Innoc_vec':Innoc_array}
            for key in attr.keys():
                attr[key][i,j,:] = one_tact_output[key] # update these variables

    T_vec      = one_tact_output['t_vec']
    dictionary = {'LTY': LTY, 'TY': TY, 'FY': FY, 'Yield': Yield, 'Res_array_1': Res_array_1, 'Res_array_2': Res_array_2, 'PRR_array': PRR_array, 'PRS_array': PRS_array, 'PSR_array': PSR_array, 'PSS_array': PSS_array, 'Selection_array_1': Selection_array_1, 'Selection_array_2': Selection_array_2, 'Innoc_array': Innoc_array,'t_vec': T_vec}
    return dictionary


#----------------------------------------------------------------------------------------------
def Z_metric(M1,M2):
    Z = np.zeros((M1.shape[0],M1.shape[1]))
    for i in range(M1.shape[0]):
        for j in range(M1.shape[1]):
            Z[i,j] = M1[i,j]/(M1[i,j] + M2[i,j])
    return Z


#----------------------------------------------------------------------------------------------
def Dose_tuplet_extractor(Selection_array_1,Selection_array_2,Res_array_1,Res_array_2,Yield,i_vec,j_vec,n_doses,separate = None):
    cmap = plt.get_cmap('jet')
    k = 0
    if separate == 'iterate':            
        R_tup1, R_tup2, SR_tup1, SR_tup2, Y_tup, L_tup, C_tup = [[None]*(len(i_vec)*len(j_vec)) for kk in range(7)]
        for i in i_vec:
            for j in j_vec:
                l = k/(len(i_vec)*len(j_vec))
                ii = floor(i*(n_doses-1))
                jj = floor(j*(n_doses-1))
                SR_tup1[k] = Selection_array_1[ii,jj,1:]
                SR_tup2[k] = Selection_array_2[ii,jj,1:]
                R_tup1[k] = Res_array_1[ii,jj,:]
                R_tup2[k] = Res_array_2[ii,jj,:]
                Y_tup[k] =Yield[ii,jj,:]
                L_tup[k] = "Dose %s and %s" % (round(ii/(n_doses-1),4),round(jj/(n_doses-1),4))
                C_tup[k] = cmap(l)
                k = k+1     
    else:
        R_tup1, R_tup2, SR_tup1, SR_tup2, Y_tup, L_tup, C_tup = [[None]*(len(i_vec)) for kk in range(7)]
        for i in range(len(i_vec)):
            l = k/(len(i_vec))
            ii = floor(i_vec[i]*(n_doses-1))
            jj = floor(j_vec[i]*(n_doses-1))
            SR_tup1[k] = Selection_array_1[ii,jj,1:]
            SR_tup2[k] = Selection_array_2[ii,jj,1:]
            R_tup1[k] = Res_array_1[ii,jj,:]
            R_tup2[k] = Res_array_2[ii,jj,:]
            Y_tup[k] =Yield[ii,jj,:]
            L_tup[k] = "Dose %s and %s" % (ii/(n_doses-1),jj/(n_doses-1))
            C_tup[k] = cmap(l)
            k = k+1     
    return SR_tup1,SR_tup2,R_tup1,R_tup2,Y_tup,L_tup,C_tup






    # from Run2_...
#----------------------------------------------------------------------------------------------
def object_dump(file_name,object_to_dump,object_type=None):
    # check if file path exists - if not create
    outdir =  os.path.dirname(file_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir,exist_ok=True) 
    #default pickle
    if object_type is None:
        object_type='pickle'
    
    if object_type == 'pickle':
        with open(file_name, 'wb') as handle:
            pickle.dump(object_to_dump,handle,protocol=pickle.HIGHEST_PROTOCOL) # protocol?
    elif object_type=='json':
        with open(file_name, 'w') as handle:
            json.dump(object_to_dump,handle)
    return None


#----------------------------------------------------------------------------------------------
def object_open(file_name,object_type=None):
    #default pickle
    if object_type is None:
        object_type='pickle'
    
    if object_type=='pickle':
        object_to_load = pickle.load(open(file_name, 'rb'))
    if object_type=='json':
        object_to_load = json.load(open(file_name))
    return object_to_load


#----------------------------------------------------------------------------------------------
def cluster_chunk(i,asex_dictionary,param_string_recursion,YDF=None):
    if YDF is None:
        YDF = calculator_d_free()
    asex_dictionary['phi_rr_val'] = asex_dictionary['phi_vec_rr'][i]
    rec_string  = params.pickle_path + 'rec_logged' + param_string_recursion + ',phi_rr_val=' + str(round(asex_dictionary['phi_rr_val'],2)) + '.pickle'
    #----------------------------------------------------------------------------------------------
    n_p = asex_dictionary['phi_vec'].shape[0]
    n_d = asex_dictionary['n_d']
    prr2, prs2, psr2, Yield = [2*np.ones((n_p,n_p,n_d,n_d)) for ii in range(4)]
    for j in range(n_p):
        for k in range(n_p):
            prr = 10**(asex_dictionary['phi_rr_val'])
            prs = 10**(asex_dictionary['phi_vec'][j])
            psr = 10**(asex_dictionary['phi_vec'][k])
            pss = 1 - prr - prs - psr
            if pss>0:# 
                output = master_loop_grid_of_tactics(n_d,1,p_rr=prr,p_rs=prs,p_sr=psr,p_ss=pss,YDF=YDF,within_season_before=0)
                prr2[j,k,:,:]  = output['PRR_array'][:,:,1] # only one season
                prs2[j,k,:,:]  = output['PRS_array'][:,:,1] # only one season
                psr2[j,k,:,:]  = output['PSR_array'][:,:,1] # only one season
                Yield[j,k,:,:] = output['Yield'][:,:,0]     # only one season
    #----------------------------------------------------------------------------------------------
    dictionary = {'prr2': prr2, 'prs2': prs2, 'psr2': psr2, 'Yield': Yield}
    ##
    rec_dict_to_dump = {**dictionary, **asex_dictionary, **params_dict}
    
    object_dump(rec_string,rec_dict_to_dump)
    
    return None
#----------------------------------------------------------------------------------------------
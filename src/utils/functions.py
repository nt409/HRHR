import numpy as np
from math import exp, ceil, floor
from scipy.integrate import odeint, simps, ode
import matplotlib.pyplot as plt
import pickle
import json
import os
import pdb

from .params import PARAMS, params_dict

# * TOC
# Utility functions


#----------------------------------------------------------------------------------------------
# Utility functions

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



# End of utility functions





#----------------------------------------------------------------------------------------------
class Simulator:
    def __init__(self):
        pass

    def growth(self, A, t):
        if t>=PARAMS.T_emerge:
            grw = PARAMS.r*(PARAMS.k-A)
            return grw
        else:
            return 0


    #----------------------------------------------------------------------------------------------
    def fcide(self, omega, theta, dose):
        effect = 1 - omega*(1 - exp(-theta*dose))
        return effect


    #----------------------------------------------------------------------------------------------
    def senescence(self, t):
        if t>=PARAMS.T_GS61:
            out = 0.005*((t-PARAMS.T_GS61)/(PARAMS.T_GS87-PARAMS.T_GS61)) + 0.1*exp(-0.02*(PARAMS.T_GS87-t))
            return out
        else:
            return 0


    #----------------------------------------------------------------------------------------------
    def ode_system(self, t, y):

        S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2 = y

        A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R

        dydt = [self.growth(A,t) - (self.senescence(t))*S -  S * (PARAMS.beta/A) * (  
                (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (IR + PR)
                + (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(                  PARAMS.omega_2,                    PARAMS.theta_2,Fung2)) * (IRS+PRS)
                + (self.fcide(                  PARAMS.omega_1,                    PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (ISR + PSR)
                + (self.fcide(                  PARAMS.omega_1,                    PARAMS.theta_1,Fung1)) * (self.fcide(                  PARAMS.omega_2,                    PARAMS.theta_2,Fung2)) * (IS + PS)  ),
            
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (IR + PR)   - (self.senescence(t)) * ER  - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1))*(self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ER,
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(                  PARAMS.omega_2,                    PARAMS.theta_2,Fung2)) * (IRS + PRS) - (self.senescence(t)) * ERS - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1))*(self.fcide(                  PARAMS.omega_2_L,                    PARAMS.theta_2_L,Fung2)) * ERS,
            S*(PARAMS.beta/A) * (self.fcide(                  PARAMS.omega_1,                    PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (ISR + PSR) - (self.senescence(t)) * ESR - PARAMS.gamma * (self.fcide(                  PARAMS.omega_1_L,                    PARAMS.theta_1_L,Fung1))*(self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ESR,
            S*(PARAMS.beta/A) * (self.fcide(                  PARAMS.omega_1,                    PARAMS.theta_1,Fung1)) * (self.fcide(                  PARAMS.omega_2,                    PARAMS.theta_2,Fung2)) * (IS + PS)   - (self.senescence(t)) * ES  - PARAMS.gamma * (self.fcide(                  PARAMS.omega_1_L,                    PARAMS.theta_1_L,Fung1))*(self.fcide(                  PARAMS.omega_2_L,                    PARAMS.theta_2_L,Fung2)) * ES,
            
            PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ER   -  PARAMS.mu * IR,
            PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1)) * (self.fcide(               PARAMS.omega_2_L,                 PARAMS.theta_2_L,Fung2)) * ERS  -  PARAMS.mu * IRS,
            PARAMS.gamma * (self.fcide(               PARAMS.omega_1_L,                 PARAMS.theta_1_L,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ESR  -  PARAMS.mu * ISR,
            PARAMS.gamma * (self.fcide(               PARAMS.omega_1_L,                 PARAMS.theta_1_L,Fung1)) * (self.fcide(               PARAMS.omega_2_L,                 PARAMS.theta_2_L,Fung2)) * ES   -  PARAMS.mu * IS,
            
            PARAMS.mu * (IR + IRS + ISR + IS)   +  (self.senescence(t)) * (S + ER + ERS + ESR + ES),
            
            -PARAMS.nu * PR,
            -PARAMS.nu * PRS,
            -PARAMS.nu * PSR,
            -PARAMS.nu * PS,
            
            -PARAMS.delta_1 * Fung1,
            -PARAMS.delta_2 * Fung2
            ]

        return dydt


    #----------------------------------------------------------------------------------------------
    def calculator_d_free(self): # calculates disease free yield
        y0   = [PARAMS.S_0] + [0]*(PARAMS.no_variables-1)

        sol  = ode(self.ode_system,jac=None).set_integrator('dopri5',nsteps= PARAMS.nstepz)
        #,rtol=10**(-10),atol=10**(-20))
        
        t0, t1, t2 = (PARAMS.T_emerge, PARAMS.T_GS61, PARAMS.T_GS87)
        
        n1= 1 + (t1-t0)/PARAMS.dt
        n2= 1 + (t2-t1)/PARAMS.dt
        
        c1 = ceil(n1-0.5)
        c2 = ceil(n2-0.5)
        
        tim1 = np.linspace(t0,t1,c1)
        tim2 = np.linspace(t1,t2,c2)
        
        yy1  = np.zeros((PARAMS.no_variables,len(tim1)))
        yy2  = np.zeros((PARAMS.no_variables,len(tim2)))
        
        i1=0
        i2=0
        
        #----------------------------------------------------------------------------------------------
        sol.set_initial_value(y0,t0)
        for t in tim1[1:]:
            if sol.successful():
                yy1[:,i1] = sol.y
                i1 += 1
                sol.integrate(t)
            else:
                raise RuntimeError('ode solver unsuccessful')
        
        #----------------------------------------------------------------------------------------------
        y1 = sol.y
        sol.set_initial_value(y1,t1)
        for t in tim2[1:]:
            if sol.successful():
                yy2[:,i2] = sol.y
                i2 += 1
                sol.integrate(t)
            else:
                raise RuntimeError('ode solver unsuccessful')
        
        
        yy2[:,i2] = sol.y
        
        int_df = simps(yy2[0,:],tim2)
        
        return int_df




    
    
    #----------------------------------------------------------------------------------------------
    def resist_prop_calculator(self, solution, solutiont=None, method=None):

        if method is None or method == 'final_value':
            disease = solution[-1,PARAMS.IR_ind]+solution[-1,PARAMS.IRS_ind]+solution[-1,PARAMS.ISR_ind]+solution[-1,PARAMS.IS_ind]
            Res_disease_1 = solution[-1,PARAMS.IR_ind]+solution[-1,PARAMS.IRS_ind]
            Res_disease_2 = solution[-1,PARAMS.IR_ind]+solution[-1,PARAMS.ISR_ind]
            RP1 = Res_disease_1/disease # phi_1
            RP2 = Res_disease_2/disease # phi_2
            p_rr = solution[-1,PARAMS.IR_ind]/disease
            p_rs = solution[-1,PARAMS.IRS_ind]/disease
            p_sr = solution[-1,PARAMS.ISR_ind]/disease
            p_ss = solution[-1,PARAMS.IS_ind]/disease
        if method == 'integrated':
            disease = simps(solution[:,PARAMS.IR_ind]+solution[:,PARAMS.IRS_ind]+solution[:,PARAMS.ISR_ind]+solution[:,PARAMS.IS_ind],solutiont)
            Res_disease_1 = simps(solution[:,PARAMS.IR_ind]+solution[:,PARAMS.IRS_ind],solutiont)
            Res_disease_2 = simps(solution[:,PARAMS.IR_ind]+solution[:,PARAMS.ISR_ind],solutiont)
            RP1 = Res_disease_1/disease # phi_1
            RP2 = Res_disease_2/disease # phi_2
            p_rr = simps(solution[:,PARAMS.IR_ind],solutiont)/disease
            p_rs = simps(solution[:,PARAMS.IRS_ind],solutiont)/disease
            p_sr = simps(solution[:,PARAMS.ISR_ind],solutiont)/disease
            p_ss = simps(solution[:,PARAMS.IS_ind],solutiont)/disease
        return RP1, RP2, p_rr, p_rs, p_sr, p_ss




    #----------------------------------------------------------------------------------------------
    def inoculum_value(self, solution, solutiont=None):
        
        method = PARAMS.res_prop_calc_method
        
        if method is None or method=='final_value':
            inoc = PARAMS.init_den*(PARAMS.den_frac +
                    (1-PARAMS.den_frac)*PARAMS.innoc_frac*
                    (solution[-1,PARAMS.IR_ind]+solution[-1,PARAMS.IRS_ind]+
                        solution[-1,PARAMS.ISR_ind]+solution[-1,PARAMS.IS_ind])
                    )
        
        elif method=='integrated':
            disease = simps(solution[:,PARAMS.IR_ind]+solution[:,PARAMS.IRS_ind] +
                solution[:,PARAMS.ISR_ind]+solution[:,PARAMS.IS_ind],solutiont)
            inoc = PARAMS.init_den*(PARAMS.den_frac + 
                        (1-PARAMS.den_frac)*PARAMS.innoc_frac_integral * disease)
        return inoc


    #----------------------------------------------------------------------------------------------
    def solve_ode(self,
            fung1_doses,
            fung2_doses,
            primary_inoc):

        res_prop_1 = primary_inoc['RR'] + primary_inoc['RS']
        res_prop_2 = primary_inoc['RR'] + primary_inoc['SR']
        ##
        y0   = [PARAMS.S_0] + [0]*9 + [primary_inoc['RR'],  primary_inoc['RS'], primary_inoc['SR'], primary_inoc['SS']] + [0]*2
        sol  = ode(self.ode_system,jac=None).set_integrator('dopri5',nsteps= PARAMS.nstepz)
         
        t0 = PARAMS.T_emerge
        t1 = PARAMS.T_GS32
        t2 = PARAMS.T_GS39
        t3 = PARAMS.T_GS61
        t4 = PARAMS.T_GS87
        
        n1 = 1 + (t1-t0)/PARAMS.dt
        n2 = 1 + (t2-t1)/PARAMS.dt
        n3 = 1 + (t3-t2)/PARAMS.dt
        n4 = 1 + (t4-t3)/PARAMS.dt
        
        c1 = ceil(n1-0.5)
        c2 = ceil(n2-0.5)
        c3 = ceil(n3-0.5)
        c4 = ceil(n4-0.5)

        tim1 = np.linspace(t0,t1,c1)
        tim2 = np.linspace(t1,t2,c2)
        tim3 = np.linspace(t2,t3,c3)
        tim4 = np.linspace(t3,t4,c4)
        yy1  = np.zeros((PARAMS.no_variables,len(tim1)))
        yy2  = np.zeros((PARAMS.no_variables,len(tim2)))
        yy3  = np.zeros((PARAMS.no_variables,len(tim3)))
        yy4  = np.zeros((PARAMS.no_variables,len(tim4)))
        i1=0
        i2=0
        i3=0
        i4=0

        ##
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
        y1[PARAMS.Fung1_ind] = y1[PARAMS.Fung1_ind] + fung1_doses['spray_1']
        y1[PARAMS.Fung2_ind] = y1[PARAMS.Fung2_ind] + fung2_doses['spray_1']
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
        y2[PARAMS.Fung1_ind]=y2[PARAMS.Fung1_ind] + fung1_doses['spray_2']
        y2[PARAMS.Fung2_ind]=y2[PARAMS.Fung2_ind] + fung2_doses['spray_2']
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
        res_prop_new_1, res_prop_new_2, p_rr, p_rs, p_sr, p_ss = self.resist_prop_calculator(solution,solutiont,method=PARAMS.res_prop_calc_method) # gives either integral or final value
        innoc = self.inoculum_value(solution, solutiont)
        Selection_1 = 1 # defaults
        Selection_2 = 1
        if res_prop_1 > 0:
            Selection_1    = res_prop_new_1/(res_prop_1/PARAMS.init_den)
        if res_prop_2 > 0:
            Selection_2    = res_prop_new_2/(res_prop_2/PARAMS.init_den)
        # #----------------------------------------------------------------------------------------------  
        inte = simps(yy4[PARAMS.S_ind,:]+yy4[PARAMS.ER_ind,:]+yy4[PARAMS.ERS_ind,:]+yy4[PARAMS.ESR_ind,:]+yy4[PARAMS.ES_ind,:],tim4)

        props_out = dict(
            RR = p_rr,
            SR = p_sr,
            RS = p_rs,
            SS = p_ss,
        )

        out = [Selection_1,
                Selection_2,
                res_prop_new_1,
                res_prop_new_2,
                innoc,
                props_out,
                inte,
                solution,
                solutiont]
        
        return out
















class RunModel:
    def __init__(self):
        self.simulator = Simulator()

        self.dis_free_yield = None






    #----------------------------------------------------------------------------------------------
    def primary_calculator(self,
                res_prop_1,
                res_prop_2,
                proportions= None,
                is_mixed_sex = PARAMS.mixed_sex):

        sex = dict(
            RR = res_prop_1*res_prop_2,
            RS = res_prop_1*(1-res_prop_2),
            SR = (1-res_prop_1)*res_prop_2,
            SS = (1-res_prop_1)*(1-res_prop_2)
            )
        
        if not is_mixed_sex:
            return sex

        else:
            asex = proportions
            out = {}
            for key in sex.keys():
                out[key] = PARAMS.sex_prop*sex[key] + (1-PARAMS.sex_prop)*asex[key]
            return out


    #----------------------------------------------------------------------------------------------
    def lifetime_yield(self, Y_vec, F_y):
        i = 1
        j = 0
        while i < F_y:
            j = j+Y_vec[i-1]
            i = i+1
        j = j/100
        return j

    #----------------------------------------------------------------------------------------------
    def total_yield(self, Y_vec, numberofseasons):
        i = 1
        j = 0
        while i<=numberofseasons:
            j = j+Y_vec[i-1]
            i = i+1
        j = j/100
        return j




    def master_loop_one_tactic(self,
                            dose_11_vec,
                            dose_12_vec,
                            dose_21_vec,
                            dose_22_vec,
                            res_prop_1=None,
                            res_prop_2=None,
                            yield_stopper=0,
                            p_rr=None,
                            p_rs=None,
                            p_sr=None,
                            p_ss=None,
                            within_season_before=True):
        """
        Run HRHR model for one strategy
        """

        failure_year = 0
        
        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.calculator_d_free()
        ##
        yield_vec = np.zeros(len(dose_11_vec))
        
        selection_vec_1, selection_vec_2, innoc_vec, res_vec_1, res_vec_2 = [np.zeros(len(dose_11_vec)+1) for i in range(5)]
        
        # array that has solution for each state variable for each year.
        sol_array = np.zeros((ceil((PARAMS.T_GS87-PARAMS.T_emerge)/PARAMS.dt), 
                                PARAMS.no_variables,
                                len(dose_11_vec))) 
        
        t_vec = np.zeros(ceil((PARAMS.T_GS87-PARAMS.T_emerge)/PARAMS.dt))


        if p_rr is None:
            primary_inoculum = self.primary_calculator(res_prop_1, res_prop_2, is_mixed_sex=False)

            res_vec_1[0] = res_prop_1
            res_vec_2[0] = res_prop_2

        else:
            primary_inoculum = dict(
                RR = p_rr,
                RS = p_rs,
                SR  = p_sr,
                SS = p_ss
                )
        
        primary_lists = dict(
            RR = [],
            RS = [],
            SR = [],
            SS = []
        )

        for key in primary_lists.keys():
            primary_lists[key].append(primary_inoculum[key])

        innoc_vec[0]   = PARAMS.init_den

        for i in range(len(dose_11_vec)):
            
            # stop the solver after we drop below threshold
            if not (i>0 and yield_vec[i-1]<yield_stopper): 

                innoc_in = innoc_vec[i]
                
                if within_season_before:
                    res_prop_1_in = primary_inoculum['RR'] + primary_inoculum['RS']
                    res_prop_2_in = primary_inoculum['RR'] + primary_inoculum['SR']
                    primary_inoculum = self.primary_calculator(res_prop_1_in, 
                                                    res_prop_2_in,
                                                    primary_inoculum)
                    
                    for key in primary_lists.keys():
                        primary_lists[key][-1] = primary_inoculum[key]
                
                fung1_doses = dict(
                    spray_1 = dose_11_vec[i],
                    spray_2 = dose_12_vec[i]
                    )
                
                fung2_doses = dict(
                    spray_1 = dose_21_vec[i],
                    spray_2 = dose_22_vec[i]
                    )
                
                model_inoc = {}
                for key in primary_inoculum.keys():
                    model_inoc[key] = innoc_in*primary_lists[key][-1]

                out = self.simulator.solve_ode(fung1_doses, fung2_doses, model_inoc)

                [Selection_1,
                 Selection_2,
                 res_prop_new_1,
                 res_prop_new_2,
                 innoc_out,
                 prop_out,
                 inte,
                 solution,
                 solutiont] = out

                Yield = 100*(inte/self.dis_free_yield)
                yield_vec[i] = Yield
                
                if not within_season_before:
                    res_prop_1_in = prop_out['RR'] + prop_out['RS']
                    res_prop_2_in = prop_out['RR'] + prop_out['SR']
                    prop_out = self.primary_calculator(res_prop_1_in,
                                            res_prop_2_in,
                                            prop_out)

                res_vec_1[i+1] = res_prop_new_1
                res_vec_2[i+1] = res_prop_new_2

                for key in primary_lists.keys():                
                    primary_lists[key].append(prop_out[key])

                innoc_vec[i+1] = innoc_out
                selection_vec_1[i+1] = Selection_1
                selection_vec_2[i+1] = Selection_2
                sol_array[:,:,i] = solution
                t_vec = solutiont

                if yield_vec[i]<PARAMS.yield_threshold and yield_vec[0]>PARAMS.yield_threshold and failure_year==0:
                    failure_year = i+1
        
        if min(yield_vec)>PARAMS.yield_threshold:
            failure_year = -1
        
        dictionary = {
                'res_vec_1': res_vec_1,
                'res_vec_2': res_vec_2,
                'primary_lists': primary_lists,
                'yield_vec': yield_vec,
                'innoc_vec': innoc_vec, 
                'selection_vec_1': selection_vec_1, 
                'selection_vec_2': selection_vec_2, 
                'failure_year': failure_year, 
                'sol_array': sol_array, 
                't_vec': t_vec
                }
        
        return dictionary


    #----------------------------------------------------------------------------------------------
    def master_loop_grid_of_tactics(self,
                n_doses,
                n_seasons,
                res_freq_1=None,
                res_freq_2=None,
                yield_stopper=0,
                strategy=PARAMS.strategy,
                p_rr=None,
                p_rs=None,
                p_sr=None,
                p_ss=None,
                within_season_before=True):
        """
        Run across grid
        """
        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.calculator_d_free()
        LTY, TY, FY = [np.zeros((n_doses,n_doses)) for i in range(3)]
        ##
        Yield = np.zeros((n_doses,n_doses,n_seasons))
        Res_array_1,Res_array_2,PRR_array,PRS_array,PSR_array,PSS_array,Selection_array_1,Selection_array_2, Innoc_array  = [np.zeros((n_doses,n_doses,n_seasons+1)) for i in range(9)]
        ##
        for i in range(n_doses):
            for j in range(n_doses):
                if strategy == 'mix':
                    dose_11_vec = 0.5*(i/(n_doses-1))*np.ones(n_seasons)
                    dose_12_vec = 0.5*(j/(n_doses-1))*np.ones(n_seasons)
                    dose_21_vec = 0.5*(i/(n_doses-1))*np.ones(n_seasons)
                    dose_22_vec = 0.5*(j/(n_doses-1))*np.ones(n_seasons)
                if strategy == 'alt_12':
                    dose_11_vec = (i/(n_doses-1))*np.ones(n_seasons)
                    dose_12_vec = np.zeros(n_seasons)
                    dose_21_vec = np.zeros(n_seasons)
                    dose_22_vec = (j/(n_doses-1))*np.ones(n_seasons)
                if strategy == 'alt_21':
                    dose_11_vec = np.zeros(n_seasons)
                    dose_12_vec = (j/(n_doses-1))*np.ones(n_seasons)
                    dose_21_vec = (i/(n_doses-1))*np.ones(n_seasons)
                    dose_22_vec = np.zeros(n_seasons)
                
                one_tact_output =  self.master_loop_one_tactic(dose_11_vec=dose_11_vec,
                                            dose_12_vec=dose_12_vec,
                                            dose_21_vec=dose_21_vec,
                                            dose_22_vec=dose_22_vec,
                                            res_prop_1=res_freq_1,
                                            res_prop_2= res_freq_2,
                                            yield_stopper=yield_stopper,
                                            p_rr=p_rr,
                                            p_rs=p_rs,
                                            p_sr=p_sr,
                                            p_ss=p_ss,
                                            within_season_before=within_season_before)
                
                LTY[i,j]  = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                TY[i,j]   = self.total_yield(one_tact_output['yield_vec'],n_seasons)
                FY[i,j]   = one_tact_output['failure_year']

                attr = {
                    'yield_vec': Yield,
                    'res_vec_1': Res_array_1,
                    'res_vec_2': Res_array_2,
                    'PRR': PRR_array,
                    'PRS': PRS_array,
                    'PSR': PSR_array,
                    'PSS': PSS_array,
                    'selection_vec_1': Selection_array_1,
                    'selection_vec_2': Selection_array_2,
                    'innoc_vec': Innoc_array
                    }

                for key in attr.keys():
                    attr[key][i,j,:] = one_tact_output[key] # update these variables

        T_vec = one_tact_output['t_vec']
        
        dictionary = {'LTY': LTY, 'TY': TY, 'FY': FY, 'Yield': Yield, 'Res_array_1': Res_array_1, 'Res_array_2': Res_array_2, 'PRR_array': PRR_array, 'PRS_array': PRS_array, 'PSR_array': PSR_array, 'PSS_array': PSS_array, 'Selection_array_1': Selection_array_1, 'Selection_array_2': Selection_array_2, 'Innoc_array': Innoc_array,'t_vec': T_vec}

        return dictionary









#----------------------------------------------------------------------------------------------
# def Z_metric(M1,M2):
#     Z = np.zeros((M1.shape[0],M1.shape[1]))
#     for i in range(M1.shape[0]):
#         for j in range(M1.shape[1]):
#             Z[i,j] = M1[i,j]/(M1[i,j] + M2[i,j])
#     return Z


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













    #----------------------------------------------------------------------------------------------
    def cluster_chunk(self, i, asex_dictionary, param_string_recursion):
        
        
        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.calculator_d_free()
        
        asex_dictionary['phi_rr_val'] = asex_dictionary['phi_vec_rr'][i]
        rec_string  = PARAMS.pickle_path + 'rec_logged' + param_string_recursion + ',phi_rr_val=' + str(round(asex_dictionary['phi_rr_val'],2)) + '.pickle'
        
        
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
                if pss>0:
                    output = self.master_loop_grid_of_tactics(n_d,1,p_rr=prr,p_rs=prs,p_sr=psr,p_ss=pss,within_season_before=False)
                    prr2[j,k,:,:]  = output['PRR_array'][:,:,1] # only one season
                    prs2[j,k,:,:]  = output['PRS_array'][:,:,1] # only one season
                    psr2[j,k,:,:]  = output['PSR_array'][:,:,1] # only one season
                    Yield[j,k,:,:] = output['Yield'][:,:,0]     # only one season
        #----------------------------------------------------------------------------------------------
        dictionary = {'prr2': prr2, 'prs2': prs2, 'psr2': psr2, 'Yield': Yield}
        ##
        rec_dict_to_dump = {**dictionary, **asex_dictionary, **params_dict}
        
        object_dump(rec_string, rec_dict_to_dump)
        
        return None

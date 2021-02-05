import numpy as np
from math import exp, ceil, floor, pi, sin, cos
from scipy.integrate import simps, ode
import pickle
import json
import os
# import pdb
import copy
from scipy.optimize import fsolve
import pandas as pd

from .params import PARAMS # , params_dict

# * TOC
# Utility functions
# class Simulator
# class RunModel
# Dyn Prog (?)


#----------------------------------------------------------------------------------------------
# Utility functions

def object_dump(file_name, object_to_dump, object_type=None):
    # check if file path exists - if not create
    outdir =  os.path.dirname(file_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir,exist_ok=True) 
    #default pickle
    if object_type is None:
        object_type='pickle'
    
    if object_type == 'pickle':
        with open(file_name, 'wb') as handle:
            pickle.dump(object_to_dump, handle, protocol=pickle.HIGHEST_PROTOCOL) # protocol?
    elif object_type=='json':
        with open(file_name, 'w') as handle:
            json.dump(object_to_dump, handle)
    return None


#----------------------------------------------------------------------------------------------
def object_open(file_name,object_type=None):
    #default pickle
    if object_type is None:
        object_type='pickle'
    
    if object_type=='pickle':
        object_to_load = pickle.load(open(file_name, 'rb'))
    elif object_type=='json':
        object_to_load = json.load(open(file_name))
    else:
        object_to_load = None

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
                + (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(               PARAMS.omega_2,                 PARAMS.theta_2,Fung2)) * (IRS+PRS)
                + (self.fcide(               PARAMS.omega_1,                 PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (ISR + PSR)
                + (self.fcide(               PARAMS.omega_1,                 PARAMS.theta_1,Fung1)) * (self.fcide(               PARAMS.omega_2,                 PARAMS.theta_2,Fung2)) * (IS + PS)  ),
            
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (IR + PR)   - (self.senescence(t)) * ER  - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1))*(self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ER,
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1,PARAMS.alpha_1_C*PARAMS.theta_1,Fung1)) * (self.fcide(               PARAMS.omega_2,                 PARAMS.theta_2,Fung2)) * (IRS + PRS) - (self.senescence(t)) * ERS - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*PARAMS.omega_1_L,PARAMS.alpha_1_C*PARAMS.theta_1_L,Fung1))*(self.fcide(               PARAMS.omega_2_L,                 PARAMS.theta_2_L,Fung2)) * ERS,
            S*(PARAMS.beta/A) * (self.fcide(               PARAMS.omega_1,                 PARAMS.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*PARAMS.omega_2,PARAMS.alpha_2_C*PARAMS.theta_2,Fung2)) * (ISR + PSR) - (self.senescence(t)) * ESR - PARAMS.gamma * (self.fcide(               PARAMS.omega_1_L,                 PARAMS.theta_1_L,Fung1))*(self.fcide(PARAMS.alpha_2*PARAMS.omega_2_L,PARAMS.alpha_2_C*PARAMS.theta_2_L,Fung2)) * ESR,
            S*(PARAMS.beta/A) * (self.fcide(               PARAMS.omega_1,                 PARAMS.theta_1,Fung1)) * (self.fcide(               PARAMS.omega_2,                 PARAMS.theta_2,Fung2)) * (IS + PS)   - (self.senescence(t)) * ES  - PARAMS.gamma * (self.fcide(               PARAMS.omega_1_L,                 PARAMS.theta_1_L,Fung1))*(self.fcide(               PARAMS.omega_2_L,                 PARAMS.theta_2_L,Fung2)) * ES,
            
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
    def find_disease_free_yield(self):

        y0   = [PARAMS.S_0] + [0]*(PARAMS.no_variables-1)

        sol  = ode(self.ode_system,jac=None).set_integrator('dopri5',nsteps= PARAMS.nstepz)
        
        t0 = PARAMS.T_emerge
        t1 = PARAMS.T_GS61
        t2 = PARAMS.T_GS87
        
        n1= 1 + (t1-t0)/PARAMS.dt
        n2= 1 + (t2-t1)/PARAMS.dt
        
        c1 = ceil(n1-0.5)
        c2 = ceil(n2-0.5)
        
        tim1 = np.linspace(t0,t1,c1)
        tim2 = np.linspace(t1,t2,c2)
        
        yy1  = np.zeros((PARAMS.no_variables, len(tim1)))
        yy2  = np.zeros((PARAMS.no_variables, len(tim2)))
        
        #----------------------------------------------------------------------------------------------
        sol.set_initial_value(y0,t0)
        for ind, t in enumerate(tim1[1:]):
            if sol.successful():
                yy1[:,ind] = sol.y
                sol.integrate(t)
            else:
                raise RuntimeError('ode solver unsuccessful')
        
        #----------------------------------------------------------------------------------------------
        y1 = sol.y
        sol.set_initial_value(y1,t1)
        for ind, t in enumerate(tim2[1:]):
            if sol.successful():
                yy2[:,ind] = sol.y
                sol.integrate(t)
            else:
                raise RuntimeError('ode solver unsuccessful')
        
        yy2[:, ind] = sol.y
        
        df_yield_integral = simps(yy2[0,:],tim2)

        self.dis_free_yield = df_yield_integral
        
        return df_yield_integral




    
    
    #----------------------------------------------------------------------------------------------
    def resist_prop_calculator(self, solution, solutiont=None, method=None):

        res_props_out = {}
        strain_frequencies = {}

        if method is None or method == 'final_value':
            disease = (solution[-1,PARAMS.IR_ind] + solution[-1,PARAMS.IRS_ind]
                         + solution[-1,PARAMS.ISR_ind] + solution[-1,PARAMS.IS_ind])
            
            Res_disease_1 = solution[-1,PARAMS.IR_ind] + solution[-1,PARAMS.IRS_ind]
            Res_disease_2 = solution[-1,PARAMS.IR_ind] + solution[-1,PARAMS.ISR_ind]

            # phi_1, phi_2
            res_props_out['f1'] = Res_disease_1/disease
            res_props_out['f2'] = Res_disease_2/disease
            
            strain_frequencies['RR'] = solution[-1,PARAMS.IR_ind]/disease
            strain_frequencies['RS'] = solution[-1,PARAMS.IRS_ind]/disease
            strain_frequencies['SR'] = solution[-1,PARAMS.ISR_ind]/disease
            strain_frequencies['SS'] = solution[-1,PARAMS.IS_ind]/disease

        if method == 'integrated':
            disease = simps(solution[:,PARAMS.IR_ind] + solution[:,PARAMS.IRS_ind]
                         + solution[:,PARAMS.ISR_ind] + solution[:,PARAMS.IS_ind],solutiont)
            
            Res_disease_1 = simps(solution[:,PARAMS.IR_ind] + solution[:,PARAMS.IRS_ind],solutiont)
            Res_disease_2 = simps(solution[:,PARAMS.IR_ind] + solution[:,PARAMS.ISR_ind],solutiont)
            
            # phi_1, phi_2
            res_props_out['f1'] = Res_disease_1/disease
            res_props_out['f2'] = Res_disease_2/disease
            
            strain_frequencies['RR'] = simps(solution[:,PARAMS.IR_ind],solutiont)/disease
            strain_frequencies['RS'] = simps(solution[:,PARAMS.IRS_ind],solutiont)/disease
            strain_frequencies['SR'] = simps(solution[:,PARAMS.ISR_ind],solutiont)/disease
            strain_frequencies['SS'] = simps(solution[:,PARAMS.IS_ind],solutiont)/disease

        return res_props_out, strain_frequencies




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

        ##
        y0 = [PARAMS.S_0] + [0]*9 + [primary_inoc['RR'],  primary_inoc['RS'], primary_inoc['SR'], primary_inoc['SS']] + [0]*2
        
        sol = ode(self.ode_system,jac=None).set_integrator('dopri5',nsteps= PARAMS.nstepz)

        time_list = [PARAMS.T_emerge,
            PARAMS.T_GS32,
            PARAMS.T_GS39,
            PARAMS.T_GS61,
            PARAMS.T_GS87]
        
        y0_new = None

        spraying_time = [False, True, True, False]
        which_spray = [None, "spray_1", "spray_2", None]
        
        y_list = []
        t_list = []

        sum_ns = 0

        for ii in range(len(time_list)-1):
            n = 1 + (time_list[ii+1]-time_list[ii])/PARAMS.dt
            n = floor(n)


            if ii==len(time_list)-2:
                n = 3 + PARAMS.t_points - sum_ns
            else:
                sum_ns += n

            time_vec = np.linspace(time_list[ii], time_list[ii+1], n)
        
            y_array  = np.zeros((PARAMS.no_variables, len(time_vec)))

            if y0_new is None:
                y0_new = y0
            else:
                y0_new = sol.y
                if spraying_time[ii]:                
                    key = which_spray[ii]
                    
                    y0_new[PARAMS.Fung1_ind] = y0_new[PARAMS.Fung1_ind] + fung1_doses[key]
                    y0_new[PARAMS.Fung2_ind] = y0_new[PARAMS.Fung2_ind] + fung2_doses[key]

            sol.set_initial_value(y0_new, time_vec[0])
            
            for index, t in enumerate(time_vec[1:]):
                if sol.successful():
                    y_array[:,index] = sol.y
                    sol.integrate(t)
                else:
                    raise RuntimeError('ode solver unsuccessful')
            
            if ii!=len(time_list)-2:
                y_list.append(y_array[:,:-1])
                t_list.append(time_vec[:-1])
            else:
                # final one of loop - need to add final time
                # rather than leave it for start condition 
                # of next run through loop
                y_array[:,index+1] = sol.y
                y_list.append(y_array)
                t_list.append(time_vec)

        
        #----------------------------------------------------------------------------------------------
        solutiont = np.concatenate(t_list)
        solutionTranspose  = np.concatenate(y_list,axis=1)
        solution  = np.transpose(solutionTranspose)
        

        # #----------------------------------------------------------------------------------------------  
        final_res_dict, props_out = self.resist_prop_calculator(solution, solutiont, method=PARAMS.res_prop_calc_method) # gives either integral or final value
        
        innoc = self.inoculum_value(solution, solutiont)
        
        initial_res_dict = dict(
            f1 = primary_inoc['RR'] + primary_inoc['RS'],
            f2 = primary_inoc['RR'] + primary_inoc['SR']
            ) 

        # get selection
        selection = dict(
            f1 = 1,
            f2 = 1
            )

        for key in ['f1','f2']:
            if initial_res_dict[key] > 0:
                selection[key]    = final_res_dict[key]/(initial_res_dict[key]/PARAMS.init_den)

        # get integral
        y_yield = y_list[-1]
        t_yield = t_list[-1]
        yield_integral = simps(y_yield[PARAMS.S_ind,:] + y_yield[PARAMS.ER_ind,:] 
                            + y_yield[PARAMS.ERS_ind,:] + y_yield[PARAMS.ESR_ind,:]
                            + y_yield[PARAMS.ES_ind,:],t_yield)

        out = dict(selection=selection,
                final_res_dict=final_res_dict,
                innoc=innoc,
                props_out=props_out,
                yield_integral=yield_integral,
                solution=solution,
                solutiont=solutiont)
        
        return out
















class RunModel:
    def __init__(self):
        self.simulator = Simulator()

        self.dis_free_yield = None

        self.yield_stopper = 0






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




    def master_loop_one_tactic(self, Config):
        """
        Run HRHR model for one strategy
        """
        
        if Config.load_saved:
            filename = Config.config_string
            if os.path.isfile(filename) and "single" in filename:
                loaded_run = pickle.load(open(filename, 'rb'))
                return loaded_run

        f1_doses = Config.fung1_doses
        f2_doses = Config.fung2_doses
        res_props = Config.res_props
        primary_inoculum = Config.primary_inoculum

        n_years = len(f1_doses['spray_1'])

        failure_year = 0
        
        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.find_disease_free_yield()
        
        yield_vec = np.zeros(n_years)
        
        innoc_vec = np.zeros(n_years+1)

        res_vec_dict = dict(
            f1 = [res_props['f1']],
            f2 = [res_props['f2']]
            )
        
        # array that has solution for each state variable for each year.
        sol_array = np.zeros((PARAMS.t_points, 
                                PARAMS.no_variables,
                                n_years)) 
        
        t_vec = np.zeros(PARAMS.t_points)

        if primary_inoculum is None:
            primary_inoculum = self.primary_calculator(res_props['f1'], res_props['f2'], is_mixed_sex=False)
        

        selection_vec_dict = dict(
            f1 = [],
            f2 = []
            )

        primary_lists = {}
        for key in primary_inoculum.keys():
            primary_lists[key] = []
        
        innoc_vec[0]   = PARAMS.init_den

        strain_freqs = primary_inoculum

        for i in range(n_years):
            
            # stop the solver after we drop below threshold
            if not (i>0 and yield_vec[i-1]< self.yield_stopper): 

                innoc_in = innoc_vec[i]
                
                if Config.within_season_before:
                    # sex/asex before each season
                    res_prop_1_start = strain_freqs['RR'] + strain_freqs['RS']
                    res_prop_2_start = strain_freqs['RR'] + strain_freqs['SR']
                    
                    strain_freqs = self.primary_calculator(res_prop_1_start, 
                                                    res_prop_2_start,
                                                    strain_freqs)
                    
                fung1_doses = dict(
                    spray_1 = f1_doses['spray_1'][i],
                    spray_2 = f1_doses['spray_2'][i]
                    )
                
                fung2_doses = dict(
                    spray_1 = f2_doses['spray_1'][i],
                    spray_2 = f2_doses['spray_2'][i]
                    )
                
                model_inoc_in = {}
                for key in primary_inoculum.keys():
                    model_inoc_in[key] = innoc_in*strain_freqs[key]
                    primary_lists[key].append(strain_freqs[key])
                
                output = self.simulator.solve_ode(fung1_doses, fung2_doses, model_inoc_in)

                yield_vec[i] = 100*(output['yield_integral']/self.dis_free_yield)
                
                strain_freqs = output['props_out']
                
                if not Config.within_season_before:
                    # sex/asex after each season
                    res_prop_1_end = strain_freqs['RR'] + strain_freqs['RS']
                    res_prop_2_end = strain_freqs['RR'] + strain_freqs['SR']
                    strain_freqs = self.primary_calculator(res_prop_1_end,
                                            res_prop_2_end,
                                            strain_freqs)

                innoc_vec[i+1] = output['innoc']

                for key in ['f1', 'f2']:
                    selection_vec_dict[key].append(output['selection'][key])
                    res_vec_dict[key].append(output['final_res_dict'][key])

                sol_array[:,:,i] = output['solution']
                t_vec = output['solutiont']

                if yield_vec[i]<PARAMS.yield_threshold and yield_vec[0]>PARAMS.yield_threshold and failure_year==0:
                    failure_year = i+1
        
        if min(yield_vec)>PARAMS.yield_threshold:
            failure_year = -1
        
        model_output = {
                'res_vec_dict': res_vec_dict,
                'primary_lists': primary_lists,
                'yield_vec': yield_vec,
                'innoc_vec': innoc_vec, 
                'selection_vec_dict': selection_vec_dict, 
                'failure_year': failure_year, 
                'sol_array': sol_array, 
                't_vec': t_vec
                }
        
        filename = Config.config_string
        if "single" in filename:
            object_dump(filename, model_output)
        
        return model_output


    #----------------------------------------------------------------------------------------------
    def master_loop_grid_of_tactics(self, ConfigG):
        """
        Run across grid
        """

        if ConfigG.load_saved:
            filename = ConfigG.config_string
            print(filename, 'grid')
            if os.path.isfile(filename):
                loaded_run = pickle.load(open(filename, 'rb'))
                return loaded_run
        

        n_seasons = ConfigG.n_years
        n_doses = ConfigG.n_doses

        ConfRun = copy.copy(ConfigG)

        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.find_disease_free_yield()
        
        LTY, TY, FY = [np.zeros((n_doses,n_doses)) for i in range(3)]
        
        yield_array = np.zeros((n_doses, n_doses, n_seasons))

        res_arrays = {}
        selection_arrays = {}
        for key in ['f1', 'f2']:
            res_arrays[key] = np.zeros((n_doses, n_doses, n_seasons+1))
            selection_arrays[key] = np.zeros((n_doses, n_doses, n_seasons))

        primary_strain_arrays = {}
        for key in ['RR', 'RS', 'SR', 'SS']:
            primary_strain_arrays[key] = np.zeros((n_doses, n_doses, n_seasons))
        
        inoc_array  = np.zeros((n_doses,n_doses,n_seasons+1))
        
        for i in range(n_doses):
            for j in range(n_doses):
                
                if ConfigG.strategy == 'mix':
                    dose_f1s1_vec = 0.5*(i/(n_doses-1))*np.ones(n_seasons)
                    dose_f1s2_vec = 0.5*(i/(n_doses-1))*np.ones(n_seasons)
                    dose_f2s1_vec = 0.5*(j/(n_doses-1))*np.ones(n_seasons)
                    dose_f2s2_vec = 0.5*(j/(n_doses-1))*np.ones(n_seasons)
                
                elif ConfigG.strategy == 'alt_12':
                    dose_f1s1_vec = (i/(n_doses-1))*np.ones(n_seasons)
                    dose_f1s2_vec = np.zeros(n_seasons)
                    dose_f2s1_vec = np.zeros(n_seasons)
                    dose_f2s2_vec = (j/(n_doses-1))*np.ones(n_seasons)
                
                elif ConfigG.strategy == 'alt_21':
                    dose_f1s1_vec = np.zeros(n_seasons)
                    dose_f1s2_vec = (j/(n_doses-1))*np.ones(n_seasons)
                    dose_f2s1_vec = (i/(n_doses-1))*np.ones(n_seasons)
                    dose_f2s2_vec = np.zeros(n_seasons)
                

                ConfRun.fung1_doses = dict(
                    spray_1 = dose_f1s1_vec,
                    spray_2 = dose_f1s2_vec
                    )
          
                ConfRun.fung2_doses = dict(
                    spray_1 = dose_f2s1_vec,
                    spray_2 = dose_f2s2_vec
                    )

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                LTY[i,j] = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                TY[i,j] = self.total_yield(one_tact_output['yield_vec'],n_seasons)
                FY[i,j] = one_tact_output['failure_year']


                attr = {
                    'yield_vec': yield_array,
                    'res_vec_dict': res_arrays,
                    'primary_lists': primary_strain_arrays,
                    'selection_vec_dict': selection_arrays,
                    'innoc_vec': inoc_array
                    }

                # update these variables
                for key in attr.keys():
                    if isinstance(attr[key], dict):
                        for key_ in attr[key].keys():
                            attr[key][key_][i,j,:] = one_tact_output[key][key_]
                    else:
                        attr[key][i,j,:] = one_tact_output[key]

        t_vec = one_tact_output['t_vec']
        
        grid_output = {'LTY': LTY,
                    'TY': TY,
                    'FY': FY,
                    'yield_array': yield_array,
                    'res_arrays': res_arrays,
                    'primary_strain_arrays': primary_strain_arrays,
                    'selection_arrays': selection_arrays,
                    'inoc_array': inoc_array,
                    't_vec': t_vec}
        
        filename = ConfigG.config_string
        object_dump(filename, grid_output)

        return grid_output

    
    
    def constant_effect(self, x, cont_radial, cont_perp):
        out = (1- exp(-PARAMS.theta_1*x)) * (1- exp(-PARAMS.theta_2*cont_radial*x)) - cont_perp
        return out






    def master_loop_dose_space_coordinate(self, ConfigDS):
        """
        Run across grid
        """

        if ConfigDS.load_saved:
            filename = ConfigDS.config_string
            filename = filename.replace("grid", "dose_space")
            print(filename, 'ds')

            if os.path.isfile(filename):
                loaded_run = pickle.load(open(filename, 'rb'))
                return loaded_run
        

        n_seasons = ConfigDS.n_years
        n_doses = ConfigDS.n_doses

        ConfRun = copy.copy(ConfigDS)

        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.find_disease_free_yield()
        
        LTY, TY, FY = [np.zeros((n_doses, n_doses)) for i in range(3)]
        
        yield_array = np.zeros((n_doses, n_doses, n_seasons))

        res_arrays = {}
        selection_arrays = {}
        for key in ['f1', 'f2']:
            res_arrays[key] = np.zeros((n_doses, n_doses, n_seasons+1))
            selection_arrays[key] = np.zeros((n_doses, n_doses, n_seasons))

        primary_strain_arrays = {}
        for key in ['RR', 'RS', 'SR', 'SS']:
            primary_strain_arrays[key] = np.zeros((n_doses, n_doses, n_seasons))
        
        inoc_array  = np.zeros((n_doses,n_doses,n_seasons+1))

        # contour_perp = [2**(xx) for xx in np.linspace(-6, 0, n_doses)]
        contour_perp = np.linspace(0, 1, n_doses+2)[1:-1]
        contours_radial = [2**(n) for n in range(-floor(n_doses/2), floor(n_doses/2)+1, 1)]

        f1_vals = np.zeros((n_doses, n_doses))
        f2_vals = np.zeros((n_doses, n_doses))
        
        for i in range(n_doses):
            for j in range(n_doses):


                f1_val = fsolve(self.constant_effect, args=(contours_radial[j], contour_perp[i]), x0=0.001)[0]
                f2_val = contours_radial[j]*f1_val

                # f1_val = (contour_perp[i]*contours_radial[j])**(0.5)
                # f2_val = (contour_perp[i]/contours_radial[j])**(0.5)

                if f1_val>1 or f2_val>1:
                    LTY[i,j] = None
                    TY[i,j] = None
                    FY[i,j] = None
                    continue

                if ConfigDS.strategy == 'mix':
                    dose_f1s1_vec = 0.5*f1_val*np.ones(n_seasons)
                    dose_f1s2_vec = 0.5*f1_val*np.ones(n_seasons)
                    dose_f2s1_vec = 0.5*f2_val*np.ones(n_seasons)
                    dose_f2s2_vec = 0.5*f2_val*np.ones(n_seasons)
                
                elif ConfigDS.strategy == 'alt_12':
                    dose_f1s1_vec = f1_val*np.ones(n_seasons)
                    dose_f1s2_vec = np.zeros(n_seasons)
                    dose_f2s1_vec = np.zeros(n_seasons)
                    dose_f2s2_vec = f2_val*np.ones(n_seasons)
                
                elif ConfigDS.strategy == 'alt_21':
                    dose_f1s1_vec = np.zeros(n_seasons)
                    dose_f1s2_vec = f1_val*np.ones(n_seasons)
                    dose_f2s1_vec = f2_val*np.ones(n_seasons)
                    dose_f2s2_vec = np.zeros(n_seasons)
                

                ConfRun.fung1_doses = dict(
                    spray_1 = dose_f1s1_vec,
                    spray_2 = dose_f1s2_vec
                    )
          
                ConfRun.fung2_doses = dict(
                    spray_1 = dose_f2s1_vec,
                    spray_2 = dose_f2s2_vec
                    )

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                LTY[i,j] = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                TY[i,j] = self.total_yield(one_tact_output['yield_vec'],n_seasons)
                FY[i,j] = one_tact_output['failure_year']
                f1_vals[i,j] = f1_val
                f2_vals[i,j] = f2_val


                attr = {
                    'yield_vec': yield_array,
                    'res_vec_dict': res_arrays,
                    'primary_lists': primary_strain_arrays,
                    'selection_vec_dict': selection_arrays,
                    'innoc_vec': inoc_array
                    }

                # update these variables
                for key in attr.keys():
                    if isinstance(attr[key], dict):
                        for key_ in attr[key].keys():
                            attr[key][key_][i,j,:] = one_tact_output[key][key_]
                    else:
                        attr[key][i,j,:] = one_tact_output[key]

        t_vec = one_tact_output['t_vec']
        
        grid_output = {'LTY': LTY,
                    'TY': TY,
                    'FY': FY,
                    'f1_vals': f1_vals,
                    'f2_vals': f2_vals,
                    'contour_perp': contour_perp,
                    'contours_radial': contours_radial,
                    'yield_array': yield_array,
                    'res_arrays': res_arrays,
                    'primary_strain_arrays': primary_strain_arrays,
                    'selection_arrays': selection_arrays,
                    'inoc_array': inoc_array,
                    't_vec': t_vec}
        
        filename = ConfigDS.config_string
        filename = filename.replace("grid", "dose_space")
        object_dump(filename, grid_output)

        return grid_output



    def master_loop_radial(self, ConfigDS):
        """
        Run radially in dose space
        """

        if ConfigDS.load_saved:
            filename = ConfigDS.config_string
            filename = filename.replace("grid", "radial")
            print(filename, 'rad')
            if os.path.isfile(filename):
                loaded_run = pickle.load(open(filename, 'rb'))
                return loaded_run
        

        n_seasons = ConfigDS.n_years
        n_doses = ConfigDS.n_doses

        ConfRun = copy.copy(ConfigDS)

        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.find_disease_free_yield()
                
        angles = np.linspace(0, pi/2, n_doses)
        radius = np.linspace(0, 2**(0.5), n_doses+1)[1:]

        row_list = []
        
        for i in range(n_doses):
            for j in range(n_doses):

                f1_val = radius[j]*cos(angles[i])
                f2_val = radius[j]*sin(angles[i])

                if f1_val>1 or f2_val>1:
                    continue

                if ConfigDS.strategy == 'mix':
                    dose_f1s1_vec = 0.5*f1_val*np.ones(n_seasons)
                    dose_f1s2_vec = 0.5*f1_val*np.ones(n_seasons)
                    dose_f2s1_vec = 0.5*f2_val*np.ones(n_seasons)
                    dose_f2s2_vec = 0.5*f2_val*np.ones(n_seasons)
                
                elif ConfigDS.strategy == 'alt_12':
                    dose_f1s1_vec = f1_val*np.ones(n_seasons)
                    dose_f1s2_vec = np.zeros(n_seasons)
                    dose_f2s1_vec = np.zeros(n_seasons)
                    dose_f2s2_vec = f2_val*np.ones(n_seasons)
                
                elif ConfigDS.strategy == 'alt_21':
                    dose_f1s1_vec = np.zeros(n_seasons)
                    dose_f1s2_vec = f1_val*np.ones(n_seasons)
                    dose_f2s1_vec = f2_val*np.ones(n_seasons)
                    dose_f2s2_vec = np.zeros(n_seasons)
                

                ConfRun.fung1_doses = dict(
                    spray_1 = dose_f1s1_vec,
                    spray_2 = dose_f1s2_vec
                    )
          
                ConfRun.fung2_doses = dict(
                    spray_1 = dose_f2s1_vec,
                    spray_2 = dose_f2s2_vec
                    )

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                lty = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                ty = self.total_yield(one_tact_output['yield_vec'],n_seasons)
                fy = one_tact_output['failure_year']

                row_list.append(dict(d1=f1_val,
                            d2=f2_val,
                            LTY=lty,
                            TY=ty,
                            FY=fy,
                            angle=angles[i],
                            radius=radius[j],
                            ))
                
        df_out = pd.DataFrame(row_list)

        
        filename = ConfigDS.config_string
        filename = filename.replace("grid", "radial")
        object_dump(filename, df_out)

        return df_out







# * Dyn Prog?

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
def cluster_chunk(i, asex_dictionary, param_string_recursion):
    
    
    # if self.dis_free_yield is None:
        # self.dis_free_yield = self.simulator.find_disease_free_yield()
    
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

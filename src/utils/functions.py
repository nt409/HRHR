import itertools
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
from tqdm import tqdm

from .params import PARAMS
from runHRHR.config_classes import SingleConfig

# * TOC
# Utility functions
# Changing doses fns
# Changing fcide fns
# class Simulator
# class RunModel


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


# def object_open(file_name, object_type=None):
#     #default pickle
#     if object_type is None:
#         object_type='pickle'
    
#     if object_type=='pickle':
#         object_to_load = pickle.load(open(file_name, 'rb'))
#     elif object_type=='json':
#         object_to_load = json.load(open(file_name))
#     else:
#         object_to_load = None

#     return object_to_load


# End of utility functions


# * changing doses fns

def get_SR_by_doses(doses, freqs):
    outputs = {}
    for dose, rf in itertools.product(doses, freqs):
        ConfigSingleRun = SingleConfig(1, rf, rf, dose, dose, dose, dose)
        output = RunModel().master_loop_one_tactic(ConfigSingleRun)
        outputs[f"dose={dose},rf={rf}"] = output


    conf_str = ConfigSingleRun.config_string_img
    str_freqs = [str(round(f,2)) for f in freqs]
    str_doses = [str(round(d,2)) for d in doses]

    middle_string = ("=" + ",_".join(str_freqs) +
                "_doses=" + ",_".join(str_doses))
    middle_string = middle_string.replace(".", ",")

    conf_str = ("=".join(conf_str.split("=")[0:2]) + 
            middle_string + conf_str.split("=")[-1])

    
    return outputs, conf_str

# End of changing doses fns



# * changing fcide fns

def get_cf_filename(Config, vars, names):
    c_string = Config.config_string.split(".pickle")[0]
    c_string = c_string.replace("single/", "changing_fcide/")

    out = c_string
    my_folder = ""
    for var, name in zip(vars,names):
        my_str = f"_{name}={round(var[0],2)},{round(var[-1],2)},{round(len(var),2)}"
        out += my_str
        my_folder += name + "_"
    
    out = out.replace(".", ",")
    out = out.replace(",,/", "../")

    out_file = out + ".pickle"

    out_img = out + ".png"
    out_img = out_img.replace("saved_runs/changing_fcide/", 
                    f"figures/changing_fcide/{my_folder[:-1]}/")
    
    return out_file, out_img




def process_changing_fcide(df, xname, yname, zname, wname=None):

    x = df[xname].unique()
    y = df[yname].unique()

    z = -2*np.ones((len(y), len(x)))

    x_inds = range(len(x))
    y_inds = range(len(y))

    for i, j in itertools.product(x_inds, y_inds):
        zz = df[(df[xname]==x[i]) & (df[yname]==y[j])]
        filtered = zz[zname]
        if len(filtered):
            if wname is None:
                z[j, i] = float(filtered)
            else:
                if max(zz[wname])>0:
                    w_filt = zz[zz[wname]==max(zz[wname])]
                    best_ones = w_filt[zname]
                    print(w_filt)
                    use = np.mean(best_ones)
                    z[j, i] = float(use)
                else:
                    z[j, i] = None
    
    return x, y, z





def changing_fcide_dose_curve(doses, curvatures, rf, NY):

    ConfigSingleRun = SingleConfig(NY, rf, rf, 0, 0, 0, 0)
    
    filename, filename_img = get_cf_filename(ConfigSingleRun,
                                [doses,
                                curvatures],
                                ["doses",
                                "curv"])
    
    if ConfigSingleRun.load_saved and os.path.isfile(filename):
        print("loading df")
        df = pickle.load(open(filename, 'rb'))
    else:
        print("running to find df")


        rows = []
        for dose, curve in itertools.product(doses, curvatures):
            ConfigSingleRun = SingleConfig(NY, rf, rf, dose, dose, dose, dose)
            
            ConfigSingleRun.load_saved = False
            
            fungicide_params = dict(
                omega_1 = PARAMS.omega_1,
                omega_2 = PARAMS.omega_2,
                theta_1 = curve,
                theta_2 = curve,
            )
            output = RunModel(fungicide_params).master_loop_one_tactic(ConfigSingleRun)
            
            FY = output['failure_year']
            rows.append(dict(dose=dose, curve=curve, failure_year=FY))
        
        df = pd.DataFrame(rows)

        object_dump(filename, df)

    x, y, z = process_changing_fcide(df, 'dose', 'curve', 'failure_year')
    
    return x, y, z, filename_img


def changing_fcide_curve_asymp(curvatures, asymps, rf, NY):
    
    ConfigSingleRun = SingleConfig(NY, rf, rf, 1, 1, 1, 1)
    
    filename, filename_img = get_cf_filename(ConfigSingleRun,
                                [curvatures,
                                asymps],
                                ["curv",
                                "asymp"])
    
    if ConfigSingleRun.load_saved and os.path.isfile(filename):
        print("loading df")
        df = pickle.load(open(filename, 'rb'))
    else:
        print("running to find df")
        
        ConfigSingleRun.load_saved = False

        rows = []
        for curve, asymp in itertools.product(curvatures, asymps):
            
            fungicide_params = dict(
                omega_1 = asymp,
                omega_2 = asymp,
                theta_1 = curve,
                theta_2 = curve,
            )
            output = RunModel(fungicide_params).master_loop_one_tactic(ConfigSingleRun)
            FY = output['failure_year']
            rows.append(dict(curve=curve, asymp=asymp, failure_year=FY))

        df = pd.DataFrame(rows)

        object_dump(filename, df)

    x, y, z = process_changing_fcide(df, 'curve', 'asymp', 'failure_year')
    
    return x, y, z, filename_img



def changing_fcide_sexp_asymp_curv(sex_props,
                            asymps, curvatures, rf, NY):
    
    ConfigSingleRun = SingleConfig(NY, rf, rf, 1, 1, 1, 1)

    filename, filename_img = get_cf_filename(ConfigSingleRun,
                                [sex_props,
                                asymps,
                                curvatures],
                                ["sex-p",
                                "asymp",
                                "curv"])
    
    if ConfigSingleRun.load_saved and os.path.isfile(filename):
        print("loading df")
        df = pickle.load(open(filename, 'rb'))
    else:
        print("running to find df")
        rows = []
        ConfigSingleRun.load_saved = False

        for sex_p, curve, asymp in tqdm(itertools.product(sex_props, curvatures, asymps)):

            ConfigSingleRun.sex_prop = sex_p
            
            fungicide_params = dict(
                omega_1 = asymp,
                omega_2 = asymp,
                theta_1 = curve,
                theta_2 = curve,
            )
            output = RunModel(fungicide_params).master_loop_one_tactic(ConfigSingleRun)
            FY = output['failure_year']
            rows.append(dict(sex_p=sex_p, curve=curve, asymp=asymp, failure_year=FY))

        df = pd.DataFrame(rows)

        object_dump(filename, df)
    
    x, y, z = process_changing_fcide(df, 'sex_p', 'asymp', 'curve', wname="failure_year")
    
    return x, y, z, filename_img

# End of changing fcide fns



#----------------------------------------------------------------------------------------------
class Simulator:
    def __init__(self, fungicide_params):
        if fungicide_params is None:
            self.omega_1 = PARAMS.omega_1
            self.omega_2 = PARAMS.omega_2
            
            self.theta_1 = PARAMS.theta_1
            self.theta_2 = PARAMS.theta_2

        else:
            self.omega_1 = fungicide_params['omega_1']
            self.omega_2 = fungicide_params['omega_2']
            
            self.theta_1 = fungicide_params['theta_1']
            self.theta_2 = fungicide_params['theta_2']
    
    @staticmethod
    def growth(A, t):
        if t>=PARAMS.T_emerge:
            grw = PARAMS.r*(PARAMS.k-A)
            return grw
        else:
            return 0


    #----------------------------------------------------------------------------------------------
    @staticmethod
    def fcide(omega, theta, dose):
        effect = 1 - omega*(1 - exp(-theta*dose))
        return effect


    #----------------------------------------------------------------------------------------------
    @staticmethod
    def senescence(t):
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
                  (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * (IR + PR)
                + (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * (IRS + PRS)
                + (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * (ISR + PSR)
                + (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * (IS + PS)),
            
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * (IR + PR)   - (self.senescence(t)) * ER  - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1))*(self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * ER,
            S*(PARAMS.beta/A) * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * (IRS + PRS) - (self.senescence(t)) * ERS - PARAMS.gamma * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1))*(self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * ERS,
            S*(PARAMS.beta/A) * (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * (ISR + PSR) - (self.senescence(t)) * ESR - PARAMS.gamma * (self.fcide(               self.omega_1,                 self.theta_1,Fung1))*(self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * ESR,
            S*(PARAMS.beta/A) * (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * (IS + PS)   - (self.senescence(t)) * ES  - PARAMS.gamma * (self.fcide(               self.omega_1,                 self.theta_1,Fung1))*(self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * ES,
            
            PARAMS.gamma * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * ER   -  PARAMS.mu * IR,
            PARAMS.gamma * (self.fcide(PARAMS.alpha_1*self.omega_1,PARAMS.alpha_1_C*self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * ERS  -  PARAMS.mu * IRS,
            PARAMS.gamma * (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(PARAMS.alpha_2*self.omega_2,PARAMS.alpha_2_C*self.theta_2,Fung2)) * ESR  -  PARAMS.mu * ISR,
            PARAMS.gamma * (self.fcide(               self.omega_1,                 self.theta_1,Fung1)) * (self.fcide(               self.omega_2,                 self.theta_2,Fung2)) * ES   -  PARAMS.mu * IS,
            
            PARAMS.mu * (IR + IRS + ISR + IS)   +  (self.senescence(t)) * (S + ER + ERS + ESR + ES),
            
            - PARAMS.nu * PR,
            - PARAMS.nu * PRS,
            - PARAMS.nu * PSR,
            - PARAMS.nu * PS,
            
            -PARAMS.delta_1 * Fung1,
            -PARAMS.delta_2 * Fung2
            ]

        return dydt


    #----------------------------------------------------------------------------------------------
    def find_disease_free_yield(self):

        y0   = [PARAMS.S_0] + [0]*(PARAMS.no_variables-1)

        sol  = ode(self.ode_system, jac=None).set_integrator('dopri5', nsteps=PARAMS.nstepz)
        
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
    @staticmethod
    def resist_prop_calculator(solution, solutiont=None, method=None):

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
    @staticmethod
    def inoculum_value(solution, solutiont=None):
        
        method = PARAMS.res_prop_calc_method
        
        if method is None or method=='final_value':
            inoc = PARAMS.init_den*((1-PARAMS.last_year_prop) +
                    PARAMS.last_year_prop*PARAMS.inoc_frac*
                    (solution[-1,PARAMS.IR_ind]+solution[-1,PARAMS.IRS_ind]+
                        solution[-1,PARAMS.ISR_ind]+solution[-1,PARAMS.IS_ind])
                    )
        
        elif method=='integrated':
            disease = simps(solution[:,PARAMS.IR_ind]+solution[:,PARAMS.IRS_ind] +
                solution[:,PARAMS.ISR_ind]+solution[:,PARAMS.IS_ind],solutiont)
            inoc = PARAMS.init_den*((1-PARAMS.last_year_prop) + 
                        PARAMS.last_year_prop*PARAMS.inoc_frac_integral * disease)
        
        else:
            raise Exception("inoculum method incorrect")

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
        
        inoc = self.inoculum_value(solution, solutiont)
        
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
                selection[key] = final_res_dict[key]/(initial_res_dict[key]/PARAMS.init_den)

        # get integral
        y_yield = y_list[-1]
        t_yield = t_list[-1]
        yield_integral = simps(y_yield[PARAMS.S_ind,:] + y_yield[PARAMS.ER_ind,:] 
                            + y_yield[PARAMS.ERS_ind,:] + y_yield[PARAMS.ESR_ind,:]
                            + y_yield[PARAMS.ES_ind,:],t_yield)

        out = dict(selection=selection,
                final_res_dict=final_res_dict,
                inoc=inoc,
                props_out=props_out,
                yield_integral=yield_integral,
                solution=solution,
                solutiont=solutiont)
        
        return out
















class RunModel:
    def __init__(self, fungicide_params=None):

        self.simulator = Simulator(fungicide_params)

        self.dis_free_yield = None

        self.yield_stopper = 0






    #----------------------------------------------------------------------------------------------
    @staticmethod
    def primary_calculator(
                Config,
                res_prop_1,
                res_prop_2,
                proportions=None,
                ):
        
        # is_mixed_sex = Config.is_mixed_sex
        sex_prop = Config.sex_prop

        sex = dict(
            RR = res_prop_1*res_prop_2,
            RS = res_prop_1*(1-res_prop_2),
            SR = (1-res_prop_1)*res_prop_2,
            SS = (1-res_prop_1)*(1-res_prop_2)
            )

        if proportions is None:
            return sex

        else:
            asex = proportions
            out = {}
            for key in sex.keys():
                out[key] = sex_prop*sex[key] + (1 - sex_prop)*asex[key]
            
            return out



    @staticmethod
    def lifetime_yield(Y_vec, F_y):
        return sum(Y_vec[:(F_y+1)])/100


    @staticmethod
    def total_yield(Y_vec):
        return sum(Y_vec)/100


    @staticmethod
    def get_doses(my_strategy, f1_val, f2_val, n_doses, n_seasons):
        if n_doses is not None:
            conc_f1 = f1_val/(n_doses-1)
            conc_f2 = f2_val/(n_doses-1)
        else:
            conc_f1 = f1_val
            conc_f2 = f2_val

        if my_strategy == 'mix':
            doses_f1s1 = 0.5*conc_f1*np.ones(n_seasons)
            doses_f1s2 = 0.5*conc_f1*np.ones(n_seasons)
            doses_f2s1 = 0.5*conc_f2*np.ones(n_seasons)
            doses_f2s2 = 0.5*conc_f2*np.ones(n_seasons)
                
        elif my_strategy == 'alt_12':
            doses_f1s1 = conc_f1*np.ones(n_seasons)
            doses_f1s2 = np.zeros(n_seasons)
            doses_f2s1 = np.zeros(n_seasons)
            doses_f2s2 = conc_f2*np.ones(n_seasons)
        
        elif my_strategy == 'alt_21':
            doses_f1s1 = np.zeros(n_seasons)
            doses_f1s2 = conc_f1*np.ones(n_seasons)
            doses_f2s1 = conc_f2*np.ones(n_seasons)
            doses_f2s2 = np.zeros(n_seasons)
        
        else:
            raise Exception("incorrect strategy named")
        
        fung1_doses = dict(
                    spray_1 = doses_f1s1,
                    spray_2 = doses_f1s2
                    )
          
        fung2_doses = dict(
            spray_1 = doses_f2s1,
            spray_2 = doses_f2s2
            )
        return fung1_doses, fung2_doses



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
        
        inoc_vec = np.zeros(n_years+1)

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
            primary_inoculum = self.primary_calculator(Config, res_props['f1'], res_props['f2'])
        

        selection_vec_dict = dict(
            f1 = [],
            f2 = []
            )

        # post-sex from previous year
        start_of_season = {}
        # pre-sex
        end_of_season = {}
        for key in primary_inoculum.keys():
            start_of_season[key] = np.zeros(n_years+1)
            end_of_season[key] = np.zeros(n_years+1)
        
        inoc_vec[0]   = PARAMS.init_den

        strain_freqs = primary_inoculum

        for i in range(n_years):
            
            # stop the solver after we drop below threshold
            if not (i>0 and yield_vec[i-1] < self.yield_stopper): 

                inoc_in = inoc_vec[i]
                    
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
                    model_inoc_in[key] = inoc_in*strain_freqs[key]
                    start_of_season[key][i] = strain_freqs[key]
                
                output = self.simulator.solve_ode(fung1_doses, fung2_doses, model_inoc_in)

                yield_vec[i] = 100*(output['yield_integral']/self.dis_free_yield)
                
                freqs_out = output['props_out']
                for key in freqs_out.keys():
                    end_of_season[key][i] = freqs_out[key]
                
                # sex/asex after each season
                res_prop_1_end = output['props_out']['RR'] + output['props_out']['RS']
                res_prop_2_end = output['props_out']['RR'] + output['props_out']['SR']

                # get next year's primary inoc - after SR step
                strain_freqs = self.primary_calculator(
                                        Config,
                                        res_prop_1_end,
                                        res_prop_2_end,
                                        freqs_out)

                inoc_vec[i+1] = output['inoc']
                
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
                'start_of_season': start_of_season,
                'end_of_season': end_of_season,
                'yield_vec': yield_vec,
                'inoc_vec': inoc_vec, 
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

        start_freqs = {}
        for key in ['RR', 'RS', 'SR', 'SS']:
            start_freqs[key] = np.zeros((n_doses, n_doses, n_seasons+1))
        
        inoc_array  = np.zeros((n_doses,n_doses,n_seasons+1))
        
        for f1_ind in tqdm(range(n_doses)):
            for f2_ind in range(n_doses):

                fung1_doses, fung2_doses = self.get_doses(
                                ConfigG.strategy, f1_ind, f2_ind, n_doses, n_seasons)               

                ConfRun.fung1_doses = fung1_doses
                ConfRun.fung2_doses = fung2_doses

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                
                LTY[f1_ind,f2_ind] = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                TY[f1_ind,f2_ind] = self.total_yield(one_tact_output['yield_vec'])
                FY[f1_ind,f2_ind] = one_tact_output['failure_year']


                attr = {
                    'yield_vec': yield_array,
                    'res_vec_dict': res_arrays,
                    'start_of_season': start_freqs,
                    'selection_vec_dict': selection_arrays,
                    'inoc_vec': inoc_array
                    }

                # update these variables
                for key in attr.keys():
                    if isinstance(attr[key], dict):
                        for key_ in attr[key].keys():
                            attr[key][key_][f1_ind,f2_ind,:] = one_tact_output[key][key_]
                    else:
                        attr[key][f1_ind,f2_ind,:] = one_tact_output[key]

        t_vec = one_tact_output['t_vec']
        
        grid_output = {'LTY': LTY,
                    'TY': TY,
                    'FY': FY,
                    'yield_array': yield_array,
                    'res_arrays': res_arrays,
                    'start_freqs': start_freqs,
                    'selection_arrays': selection_arrays,
                    'inoc_array': inoc_array,
                    't_vec': t_vec}
        
        filename = ConfigG.config_string
        object_dump(filename, grid_output)

        return grid_output

    
    @staticmethod
    def constant_effect(x, cont_radial, cont_perp):
        print("nb this uses PARAMS curvatures and omega=1")
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

        start_freqs = {}
        for key in ['RR', 'RS', 'SR', 'SS']:
            start_freqs[key] = np.zeros((n_doses, n_doses, n_seasons))
        
        inoc_array  = np.zeros((n_doses,n_doses,n_seasons+1))

        # contour_perp = [2**(xx) for xx in np.linspace(-6, 0, n_doses)]
        contour_perp = np.linspace(0, 1, n_doses+2)[1:-1]
        contours_radial = [2**(n) for n in range(-floor(n_doses/2), floor(n_doses/2)+1, 1)]

        f1_vals = np.zeros((n_doses, n_doses))
        f2_vals = np.zeros((n_doses, n_doses))
        
        for i in tqdm(range(n_doses)):
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

                fung1_doses, fung2_doses = self.get_doses(
                                        ConfigDS.strategy, f1_val, f2_val, None, n_seasons)

                ConfRun.fung1_doses = fung1_doses
                ConfRun.fung2_doses = fung2_doses

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                LTY[i,j] = self.lifetime_yield(one_tact_output['yield_vec'],one_tact_output['failure_year'])
                TY[i,j] = self.total_yield(one_tact_output['yield_vec'])
                FY[i,j] = one_tact_output['failure_year']
                f1_vals[i,j] = f1_val
                f2_vals[i,j] = f2_val


                attr = {
                    'yield_vec': yield_array,
                    'res_vec_dict': res_arrays,
                    'start_of_season': start_freqs,
                    'selection_vec_dict': selection_arrays,
                    'inoc_vec': inoc_array
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
                    'start_freqs': start_freqs,
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
        n_angles = ConfigDS.n_angles
        n_radii = ConfigDS.n_radii

        ConfRun = copy.copy(ConfigDS)

        if self.dis_free_yield is None:
            self.dis_free_yield = self.simulator.find_disease_free_yield()
                
        angles = np.linspace(0, pi/2, n_angles)
        radii = np.linspace(0, 2**(0.5), n_radii+1)[1:]

        row_list = []
        
        for angle in tqdm(angles):
            for radius in radii:

                f1_val = radius*cos(angle)
                f2_val = radius*sin(angle)

                if radius==2**0.5 and angle==pi/4:
                    # floating point error meant this important point missed without this
                    f1_val=1
                    f2_val=1

                if f1_val>1 or f2_val>1:
                    continue
                
                fung1_doses, fung2_doses = self.get_doses(
                                        ConfigDS.strategy, f1_val, f2_val, None, n_seasons)

                ConfRun.fung1_doses = fung1_doses
                ConfRun.fung2_doses = fung2_doses

                one_tact_output =  self.master_loop_one_tactic(ConfRun)
                
                lty = self.lifetime_yield(one_tact_output['yield_vec'], one_tact_output['failure_year'])
                ty = self.total_yield(one_tact_output['yield_vec'])
                fy = one_tact_output['failure_year']

                row_list.append(dict(d1=f1_val,
                            d2=f2_val,
                            LTY=lty,
                            TY=ty,
                            FY=fy,
                            angle=angle,
                            radius=radius,
                            ))
                
        df_out = pd.DataFrame(row_list)

        
        filename = ConfigDS.config_string
        filename = filename.replace("grid", "radial")
        object_dump(filename, df_out)

        return df_out






# End of RunModel class






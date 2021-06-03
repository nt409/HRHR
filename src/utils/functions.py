import itertools
import numpy as np
from math import exp, ceil, floor, pi, sin, cos, log10
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
# Param Scan fns
# class Simulator
# class RunSingleTactic
# classes based on RunMultipleTactics


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


def logit10(x):
    if x>0 and x<1:
        return log10(x/(1-x))
    else:
        raise Exception(f"x={x} - invalid value")

def log10_difference(x1, x2):
    return log10(x1) - log10(x2)

def logit10_difference(x1, x2):
    return logit10(x1) - logit10(x2)


# End of utility functions


# * changing doses fns

def get_SR_by_doses(doses, freqs):
    outputs = {}
    for dose, rf in itertools.product(doses, freqs):
        ConfigSingleRun = SingleConfig(1, rf, rf, dose, dose, dose, dose)
        output = RunSingleTactic().run_single_tactic(ConfigSingleRun)
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
    """
    Convert dataframe into x,y,z output to plot
    """

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
                delta_1 = PARAMS.delta_1,
                delta_2 = PARAMS.delta_2,
            )
            output = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
            
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
                delta_1 = PARAMS.delta_1,
                delta_2 = PARAMS.delta_2,
            )
            output = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
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
                delta_1 = PARAMS.delta_1,
                delta_2 = PARAMS.delta_2,
            )
            output = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
            FY = output['failure_year']
            rows.append(dict(sex_p=sex_p, curve=curve, asymp=asymp, failure_year=FY))

        df = pd.DataFrame(rows)

        object_dump(filename, df)
    
    x, y, z = process_changing_fcide(df, 'sex_p', 'asymp', 'curve', wname="failure_year")
    
    return x, y, z, filename_img

# End of changing fcide fns

#----------------------------------------------------------------------------------------------
class Fungicide:
    def __init__(self, omega, theta, delta):
        self.omega = omega
        self.theta = theta
        self.delta = delta

    def effect(self, conc):
        effect = 1 - self.omega*(1 - exp(- self.theta * conc))
        return effect


# * Simulator Class

class Simulator:
    def __init__(self, fungicide_params):
        if fungicide_params is None:
            omega_1 = PARAMS.omega_1
            omega_2 = PARAMS.omega_2
            
            theta_1 = PARAMS.theta_1
            theta_2 = PARAMS.theta_2

            delta_1 = PARAMS.delta_1
            delta_2 = PARAMS.delta_2

        else:
            omega_1 = fungicide_params['omega_1']
            omega_2 = fungicide_params['omega_2']
            
            theta_1 = fungicide_params['theta_1']
            theta_2 = fungicide_params['theta_2']

            delta_1 = fungicide_params['delta_1']
            delta_2 = fungicide_params['delta_2']

        self.fcide1 = Fungicide(omega_1, theta_1, delta_1)
        self.fcide2 = Fungicide(omega_2, theta_2, delta_2)
    

    @staticmethod
    def _growth(A, t):
        if t>=PARAMS.T_emerge:
            grw = PARAMS.r*(PARAMS.k-A)
            return grw
        else:
            return 0


    @staticmethod
    def _senescence(t):
        if t>=PARAMS.T_GS61:
            out = 0.005*((t-PARAMS.T_GS61)/(PARAMS.T_GS87-PARAMS.T_GS61)) + 0.1*exp(-0.02*(PARAMS.T_GS87-t))
            return out
        else:
            return 0


    def _ode_system(self, t, y):

        S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,conc_1,conc_2 = y

        A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R

        dydt = [self._growth(A,t) - (self._senescence(t))*S -  S * (PARAMS.beta/A) * (  
                  (IR + PR)
                + (IRS + PRS) * (self.fcide2.effect(conc_2))
                + (ISR + PSR) * (self.fcide1.effect(conc_1))
                + (IS  +  PS) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2))),
            
            S*(PARAMS.beta/A) * (IR + PR) - (self._senescence(t)) * ER  - PARAMS.gamma * ER,
            S*(PARAMS.beta/A) * (IRS + PRS) * (self.fcide2.effect(conc_2)) - (self._senescence(t)) * ERS - PARAMS.gamma * (self.fcide2.effect(conc_2)) * ERS,
            S*(PARAMS.beta/A) * (ISR + PSR) * (self.fcide1.effect(conc_1)) - (self._senescence(t)) * ESR - PARAMS.gamma * (self.fcide1.effect(conc_1)) * ESR,
            S*(PARAMS.beta/A) * (IS  +  PS) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2)) - (self._senescence(t)) * ES  - PARAMS.gamma * (self.fcide1.effect(conc_1))*(self.fcide2.effect(conc_2)) * ES,
            
            PARAMS.gamma * ER   -  PARAMS.mu * IR,
            PARAMS.gamma * (self.fcide2.effect(conc_2)) * ERS  -  PARAMS.mu * IRS,
            PARAMS.gamma * (self.fcide1.effect(conc_1)) * ESR  -  PARAMS.mu * ISR,
            PARAMS.gamma * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2)) * ES   -  PARAMS.mu * IS,
            
            PARAMS.mu * (IR + IRS + ISR + IS)   +  (self._senescence(t)) * (S + ER + ERS + ESR + ES),
            
            - PARAMS.nu * PR,
            - PARAMS.nu * PRS,
            - PARAMS.nu * PSR,
            - PARAMS.nu * PS,
            
            -self.fcide1.delta * conc_1,
            -self.fcide2.delta * conc_2
            ]

        return dydt


    

    
    def _final_value_res_props(self):
        """
        Uses final value (end of season) to determine the res props. 
        These are used for next season (with a SR step in between if sr_prop=/=0)
        """

        disease = (self.solution[-1,PARAMS.IR_ind] + 
                        self.solution[-1,PARAMS.IRS_ind] +
                        self.solution[-1,PARAMS.ISR_ind] + 
                        self.solution[-1,PARAMS.IS_ind])
            
        Res_disease_1 = self.solution[-1,PARAMS.IR_ind] + self.solution[-1,PARAMS.IRS_ind]
        Res_disease_2 = self.solution[-1,PARAMS.IR_ind] + self.solution[-1,PARAMS.ISR_ind]

        res_props_out = dict(
            f1 = Res_disease_1/disease,
            f2 = Res_disease_2/disease,
            )
        
        strain_frequencies = dict(
            RR = self.solution[-1,PARAMS.IR_ind]/disease,
            RS = self.solution[-1,PARAMS.IRS_ind]/disease,
            SR = self.solution[-1,PARAMS.ISR_ind]/disease,
            SS = self.solution[-1,PARAMS.IS_ind]/disease
        )
        
        return res_props_out, strain_frequencies

    
    def _integrated_res_props(self):
        """
        Not used in most recent version 
        
        - idea was to integrate over end of the season,
        rather than just take the final value
        
        """

        disease = simps(self.solution[:,PARAMS.IR_ind] + self.solution[:,PARAMS.IRS_ind]
                         + self.solution[:,PARAMS.ISR_ind] + self.solution[:,PARAMS.IS_ind], self.solutiont)
            
        Res_disease_1 = simps(self.solution[:,PARAMS.IR_ind] + self.solution[:,PARAMS.IRS_ind], self.solutiont)
        Res_disease_2 = simps(self.solution[:,PARAMS.IR_ind] + self.solution[:,PARAMS.ISR_ind], self.solutiont)
        
        res_props_out = dict(
            f1 = Res_disease_1/disease,
            f2 = Res_disease_2/disease,
            )
        
        strain_frequencies = dict(
            RR = simps(self.solution[:,PARAMS.IR_ind], self.solutiont)/disease,
            RS = simps(self.solution[:,PARAMS.IRS_ind], self.solutiont)/disease,
            SR = simps(self.solution[:,PARAMS.ISR_ind], self.solutiont)/disease,
            SS = simps(self.solution[:,PARAMS.IS_ind], self.solutiont)/disease,
            )

        return res_props_out, strain_frequencies
    
    
    #----------------------------------------------------------------------------------------------
    def _calculate_res_props(self, method=None):

        if method is None or method=='final_value':
            return self._final_value_res_props()

        if method=='integrated':
            return self._integrated_res_props()


    def _integrated_inoc_val(self):
        disease = simps(self.solution[:,PARAMS.IR_ind]+self.solution[:,PARAMS.IRS_ind] +
                self.solution[:,PARAMS.ISR_ind]+self.solution[:,PARAMS.IS_ind],self.solutiont)
        return PARAMS.init_den*((1-PARAMS.last_year_prop) + 
                        PARAMS.last_year_prop * PARAMS.inoc_frac_integral * disease)


    def _final_inoc_value(self):
        return PARAMS.init_den*((1-PARAMS.last_year_prop) +
                    PARAMS.last_year_prop * PARAMS.inoc_frac*
                    (self.solution[-1,PARAMS.IR_ind]+self.solution[-1,PARAMS.IRS_ind]+
                        self.solution[-1,PARAMS.ISR_ind]+self.solution[-1,PARAMS.IS_ind])
                    ) 

    #----------------------------------------------------------------------------------------------
    def _get_inoculum_value(self):
        
        method = PARAMS.res_prop_calc_method
        
        if method is None or method=='final_value':
            return self._final_inoc_value()
        
        elif method=='integrated':
            return self._integrated_inoc_val()
        
        else:
            raise Exception("inoculum method incorrect")



    def _solve_ode(self,
            fung1_doses,
            fung2_doses,
            primary_inoc):

        initial_res_dict = dict(
            f1 = primary_inoc['RR'] + primary_inoc['RS'],
            f2 = primary_inoc['RR'] + primary_inoc['SR']
            ) 

        y0 = [PARAMS.S_0] + [0]*9 + [primary_inoc['RR'],  primary_inoc['RS'], primary_inoc['SR'], primary_inoc['SS']] + [0]*2
        
        sol = ode(self._ode_system,jac=None).set_integrator('dopri5',nsteps= PARAMS.nstepz)

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
        self.solutiont = np.concatenate(t_list)
        solutionTranspose  = np.concatenate(y_list,axis=1)
        self.solution  = np.transpose(solutionTranspose)
        

        # #----------------------------------------------------------------------------------------------  
        final_res_dict, props_out = self._calculate_res_props(method=PARAMS.res_prop_calc_method)
        
        inoc = self._get_inoculum_value()

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

        yield_integral = simps(y_yield[PARAMS.S_ind,:] + 
                            y_yield[PARAMS.ER_ind,:] 
                            + y_yield[PARAMS.ERS_ind,:] + y_yield[PARAMS.ESR_ind,:]
                            + y_yield[PARAMS.ES_ind,:],
                            t_yield)

        out = dict(selection=selection,
                final_res_dict=final_res_dict,
                inoc=inoc,
                props_out=props_out,
                yield_integral=yield_integral,
                solution=self.solution,
                solutiont=self.solutiont)
        
        return out


# * End of Simulator Class







# * Dose class

class FungicideStrategy:
    def __init__(self, my_strategy, n_seasons):
        self.n_seasons = n_seasons
        self.my_strategy = my_strategy


    def _get_mixed_doses(self):        
        # did half 0.5*
        # but Hobbelen paper just says it means twice as much
        self.fung1_doses = dict(
            spray_1 = self.conc_f1*np.ones(self.n_seasons),
            spray_2 = self.conc_f1*np.ones(self.n_seasons)
            )
        self.fung2_doses = dict(
            spray_1 = self.conc_f2*np.ones(self.n_seasons),
            spray_2 = self.conc_f2*np.ones(self.n_seasons)
            )


    def _get_alt_12_doses(self):
        self.fung1_doses = dict(
            spray_1 = self.conc_f1*np.ones(self.n_seasons),
            spray_2 = np.zeros(self.n_seasons)
            )
        self.fung2_doses = dict(
            spray_1 = np.zeros(self.n_seasons),
            spray_2 = self.conc_f2*np.ones(self.n_seasons)
            )
    

    def _get_alt_21_doses(self):
        self.fung1_doses = dict(
            spray_1 = np.zeros(self.n_seasons),
            spray_2 = self.conc_f1*np.ones(self.n_seasons)
            )
        self.fung2_doses = dict(
            spray_1 = self.conc_f2*np.ones(self.n_seasons),
            spray_2 = np.zeros(self.n_seasons)
            )



    def _set_grid_concs(self, f1_val, f2_val, n_doses):
        self.conc_f1 = f1_val/(n_doses-1)
        self.conc_f2 = f2_val/(n_doses-1)


    def _set_regular_concs(self, f1_val, f2_val):
        self.conc_f1 = f1_val
        self.conc_f2 = f2_val


    def _set_concs(self, f1_val, f2_val, n_doses):
        if n_doses is not None:
            self._set_grid_concs(f1_val, f2_val, n_doses)
        else:
            self._set_regular_concs(f1_val, f2_val)


    def _dose_for_this_strategy(self):
        if self.my_strategy=='mix':
            self._get_mixed_doses()
                
        elif self.my_strategy=='alt_12':
            self._get_alt_12_doses()
        
        elif self.my_strategy=='alt_21':
            self._get_alt_21_doses()
        
        else:
            raise Exception("incorrect strategy named")



    def get_doses(self, f1_val, f2_val, n_doses):

        self._set_concs(f1_val, f2_val, n_doses)

        self._dose_for_this_strategy()      

        return self.fung1_doses, self.fung2_doses



# * End of Dose class



# * RunSingleTactic class

class RunSingleTactic:
    def __init__(self, fcide_parms=None):

        self.sim = Simulator(fcide_parms)

        self._find_disease_free_yield()

        self.yield_stopper = 95

        self.PATHOGEN_STRAIN_NAMES = ['RR', 'RS', 'SR', 'SS']


    def _find_disease_free_yield(self):

        y0   = [PARAMS.S_0] + [0]*(PARAMS.no_variables-1)

        sol  = ode(self.sim._ode_system, jac=None).set_integrator('dopri5', nsteps=PARAMS.nstepz)
        
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


    def _initialise_freqs(self):
        out = {}
        for key in self.PATHOGEN_STRAIN_NAMES:
            out[key] = np.zeros(self.n_years+1)

        return out




    def _initialise_res_vec_dict(self, res_props):
        out = {}
        keys = ['f1', 'f2']
        
        for key in keys:
            out[key] = np.zeros(self.n_years+1)
            # set first year
            out[key][0] = res_props[key]

        return out



    def _initialise_sel_vec_dict(self):
        out = {}
        keys = ['f1', 'f2']
        
        for key in keys:
            out[key] = np.zeros(self.n_years+1)

        return out
    
    
    def _initialise_variables_single_run(self, Config):
        
        res_props = Config.res_props

        primary_inoculum = Config.primary_inoculum

        self.n_years = len(Config.fung1_doses['spray_1'])
        
        self.failure_year = 0
        
        self.yield_vec = np.zeros(self.n_years)
        
        self.inoc_vec = np.zeros(self.n_years+1)
        
        if primary_inoculum is None:
            primary_inoculum = self._primary_calculator(Config, res_props['f1'], res_props['f2'])

        self.res_vec_dict = self._initialise_res_vec_dict(res_props)


        # array that has solution for each state variable for each year.
        self.sol_array = np.zeros((PARAMS.t_points, 
                                PARAMS.no_variables,
                                self.n_years)) 
        
        self.t_vec = np.zeros(PARAMS.t_points)


        self.selection_vec_dict = self._initialise_sel_vec_dict()

        # post-sex from previous year
        self.start_of_season = self._initialise_freqs()
        # pre-sex
        self.end_of_season = self._initialise_freqs()
        
        self.inoc_vec[0] = PARAMS.init_den

        self.strain_freqs = primary_inoculum



    def update_start_of_season(self, yr):
        for key in self.PATHOGEN_STRAIN_NAMES:
            self.start_of_season[key][yr] = self.strain_freqs[key]
    


    def get_initial_freqs(self, yr):
        out = {}
        for key in self.PATHOGEN_STRAIN_NAMES:
            out[key] = self.inoc_vec[yr]*self.strain_freqs[key]
        return out



    @staticmethod
    def _load_single_tactic(filename):
        if os.path.isfile(filename) and "single" in filename:
            loaded_run = pickle.load(open(filename, 'rb'))
            return loaded_run



    @staticmethod
    def _primary_calculator(
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
    def _get_single_fung_dose(dose_vec, yr):
        return dict(
                    spray_1 = dose_vec['spray_1'][yr],
                    spray_2 = dose_vec['spray_2'][yr]
                    )


    def _update_failure_year(self, yr):
        """
        Set failure year if:
        - yield is below threshold
        - started above threshold
        - is first time it has dropped below threshold
        """
        if ((self.yield_vec[yr]<PARAMS.yield_threshold) and 
                (self.yield_vec[0]>PARAMS.yield_threshold) and 
                (self.failure_year==0)):
            self.failure_year = yr+1




    def _run_single_year(self, Config, yr):
        
        fung1_doses = self._get_single_fung_dose(Config.fung1_doses, yr)
        
        fung2_doses = self._get_single_fung_dose(Config.fung2_doses, yr)
        
        self.update_start_of_season(yr)

        model_inoc_in = self.get_initial_freqs(yr)
        
        output = self.sim._solve_ode(fung1_doses, fung2_doses, model_inoc_in)

        self._process_single_output(output, Config, yr)

        self._update_failure_year(yr)



    
    def _loop_over_years(self, Config):
        for yr in range(self.n_years):
            # stop the solver after we drop below threshold
            if not (yr>0 and self.yield_vec[yr-1]<self.yield_stopper):
                self._run_single_year(Config, yr)
        
        if min(self.yield_vec)>PARAMS.yield_threshold:
            self.failure_year = -1



    def _update_selection_vec_dict(self, output, yr):
        for key in ['f1', 'f2']:
            self.selection_vec_dict[key][yr] = output['selection'][key]


    def _update_res_vec_dict(self, output, yr):
        for key in ['f1', 'f2']:
            self.res_vec_dict[key][yr] = output['final_res_dict'][key]


    def _update_end_of_season(self, freqs_out, yr):
        for key in freqs_out.keys():
            self.end_of_season[key][yr] = freqs_out[key]
        




    def _process_single_output(self, output, Config, yr):
        """
        Update variables so that this year's contribution included
        """

        self.yield_vec[yr] = 100*(output['yield_integral']/self.dis_free_yield)

        freqs_out = output['props_out']
                
        self._update_end_of_season(freqs_out, yr)

        # sex/asex after each season
        res_prop_1_end = output['props_out']['RR'] + output['props_out']['RS']
        res_prop_2_end = output['props_out']['RR'] + output['props_out']['SR']

        # get next year's primary inoc - including SR step
        self.strain_freqs = self._primary_calculator(
                                Config,
                                res_prop_1_end,
                                res_prop_2_end,
                                freqs_out)

        self.inoc_vec[yr+1] = output['inoc']

        self._update_selection_vec_dict(output, yr+1)

        self._update_res_vec_dict(output, yr+1)
        
        self.sol_array[:,:,yr] = output['solution']

        self.t_vec = output['solutiont']






    def _save_single_run(self):
        model_output = {
                'res_vec_dict': self.res_vec_dict,
                'start_of_season': self.start_of_season,
                'end_of_season': self.end_of_season,
                'yield_vec': self.yield_vec,
                'inoc_vec': self.inoc_vec, 
                'selection_vec_dict': self.selection_vec_dict, 
                'failure_year': self.failure_year, 
                'sol_array': self.sol_array, 
                't_vec': self.t_vec
                }

        if "single" in self.filename:
            object_dump(self.filename, model_output)
        
        return model_output


    def run_single_tactic(self, Config):
        """
        Run HRHR model for one strategy
        """
        
        self.filename = Config.config_string
        
        if Config.load_saved:
            loaded_run = self._load_single_tactic(self.filename)
            if loaded_run is not None:
                return loaded_run

        self._initialise_variables_single_run(Config)

        self._loop_over_years(Config)
        
        model_output = self._save_single_run()
        
        return model_output


# End of RunSingleTactic





# * RunMultipleTactics

class RunMultipleTactics:
    def __init__(self, fcide_parms=None):
        self.sing_tact = RunSingleTactic(fcide_parms)


    @staticmethod
    def _lifetime_yield(Y_vec, F_y):
        return sum(Y_vec[:(F_y+1)])/100

    @staticmethod
    def _this_year_profit(yield_, dose1, dose2):
        
        tons_per_ha = 10
        
        # £/tonne
        price_per_ton = 117.14

        # £32.40/ha for full dose, and two applications
        c1 = 2*32.4
        c2 = 2*32.4

        # £/ha other machinery costs through the year?
        breakeven_yield = 0.95
        
        # if breakeven yield is 95%, then breakeven if
        # spray full dose and obtain 95% yield
        other_costs = tons_per_ha*price_per_ton*breakeven_yield - 20 - c1 - c2
        
        if dose1+dose2>0:
            # £20 tractor costs
            other_costs += 20
        
        tonnes = (yield_/100)*tons_per_ha
        revenue = price_per_ton*tonnes
        
        costs = dose1*c1 + dose2*c2 + other_costs

        profit = revenue - costs

        # need other_costs>0.79*tons_per_ha*price_per_ton (yield without spraying)
        
        return profit

    
    def _economic_life(self, Y_vec, dose1, dose2):
        total_profit = 0
        profit = self._this_year_profit(Y_vec[0], dose1, dose2)
        # if profit>0:
        #     total_profit += profit
        
        i = 1
        while profit>0 and i<len(Y_vec):
            total_profit += profit
            profit = self._this_year_profit(Y_vec[i], dose1, dose2)
            i += 1
        
        return total_profit



    @staticmethod
    def _total_yield(Y_vec):
        return sum(Y_vec)/100



    @staticmethod
    def _load_multi_tactic(filename):
        if os.path.isfile(filename):
            loaded_run = pickle.load(open(filename, 'rb'))
            return loaded_run
        else:
            return None



    @staticmethod
    def _get_dict_of_zero_arrays(keys, shape):
        out = {}
        for key in keys:
            out[key] = np.zeros(shape)
        return out



    def _initialise_multi_vars(self, n_doses, n_years):

        self.LTY = np.zeros((n_doses, n_doses))
        self.TY = np.zeros((n_doses, n_doses))
        self.FY = np.zeros((n_doses, n_doses))
        self.econ = np.zeros((n_doses, n_doses))
        
        self.yield_array = np.zeros((n_doses, n_doses, n_years))
        self.inoc_array  = np.zeros((n_doses, n_doses, n_years+1))
        
        fung_keys = ['f1', 'f2']
        self.res_arrays = self._get_dict_of_zero_arrays(fung_keys, (n_doses, n_doses, n_years+1))
        self.selection_arrays = self._get_dict_of_zero_arrays(fung_keys, (n_doses, n_doses, n_years+1))

        strain_keys = ['RR', 'RS', 'SR', 'SS']
        self.start_freqs = self._get_dict_of_zero_arrays(strain_keys, (n_doses, n_doses, n_years+1))
        

    
    def _get_fung_strat(self, Conf):
        self.fung_strat = FungicideStrategy(Conf.strategy, Conf.n_years)



    def _update_dict_array_this_dose(self, to_update, data, f1_ind, f2_ind, key1):


        for key_ in to_update.keys():
            to_update[key_][f1_ind,f2_ind,:] = data[key1][key_]
        
        # print(to_update[key_][f1_ind, f2_ind, :])
        
        return to_update
    
       


    @staticmethod
    def _get_total_doses_applied_this_year(Conf):
        total_dose_f1 = Conf.fung1_doses['spray_1'][0] + Conf.fung1_doses['spray_2'][0]
        total_dose_f2 = Conf.fung2_doses['spray_1'][0] + Conf.fung2_doses['spray_2'][0]
        return total_dose_f1, total_dose_f2




    
    def _post_process_multi(self, data_this_dose, f1_ind, f2_ind, Conf):

        self.LTY[f1_ind,f2_ind] = self._lifetime_yield(data_this_dose['yield_vec'],data_this_dose['failure_year'])

        self.TY[f1_ind,f2_ind] = self._total_yield(data_this_dose['yield_vec'])

        total_dose_f1, total_dose_f2 = self._get_total_doses_applied_this_year(Conf)

        self.econ[f1_ind, f2_ind] = self._economic_life(data_this_dose['yield_vec'], total_dose_f1, total_dose_f2)        
        
        self.FY[f1_ind,f2_ind] = data_this_dose['failure_year']

        self.inoc_array[f1_ind,f2_ind,:] = data_this_dose["inoc_vec"]
        self.yield_array[f1_ind,f2_ind,:] = data_this_dose["yield_vec"]

        self.res_arrays = self._update_dict_array_this_dose(copy.copy(self.res_arrays), data_this_dose, f1_ind, f2_ind, "res_vec_dict")
        self.start_freqs = self._update_dict_array_this_dose(copy.copy(self.start_freqs), data_this_dose, f1_ind, f2_ind, "start_of_season")
        self.selection_arrays = self._update_dict_array_this_dose(copy.copy(self.selection_arrays), data_this_dose, f1_ind, f2_ind, "selection_vec_dict")


    






class RunGrid(RunMultipleTactics):
    def __init__(self, fcide_parms=None):
        super().__init__(fcide_parms)

    def _run_the_grid(self, Conf):
        
        for f1_ind in tqdm(range(Conf.n_doses)):
            for f2_ind in range(Conf.n_doses):

                Conf.fung1_doses, Conf.fung2_doses = self.fung_strat.get_doses(
                                                    f1_ind, f2_ind, Conf.n_doses)

                one_tact_output =  self.sing_tact.run_single_tactic(Conf)
                
                self._post_process_multi(one_tact_output, f1_ind, f2_ind, Conf)

        self.t_vec = one_tact_output['t_vec']




    def _save_grid(self):
        grid_output = {'LTY': self.LTY,
                    'TY': self.TY,
                    'FY': self.FY,
                    'yield_array': self.yield_array,
                    'res_arrays': self.res_arrays,
                    'start_freqs': self.start_freqs,
                    'selection_arrays': self.selection_arrays,
                    'inoc_array': self.inoc_array,
                    't_vec': self.t_vec,
                    'econ': self.econ,
                    }
        
        object_dump(self.filename, grid_output)
        
        return grid_output





    def grid_of_tactics(self, ConfigG):
        """
        Run across grid
        """

        self.filename = ConfigG.config_string

        if ConfigG.load_saved:
            loaded_run = self._load_multi_tactic(self.filename)
            if loaded_run is not None:
                return loaded_run

        Conf = copy.copy(ConfigG)

        self._get_fung_strat(Conf)

        self._initialise_multi_vars(Conf.n_doses, Conf.n_years)

        self._run_the_grid(Conf)

        grid_output = self._save_grid()

        return grid_output

    







class RunDoseSpace(RunMultipleTactics):
  
    @staticmethod
    def constant_effect(x, cont_radial, cont_perp):
        print("nb this uses PARAMS curvatures and omega=1")
        out = (1- exp(-PARAMS.theta_1*x)) * (1- exp(-PARAMS.theta_2*cont_radial*x)) - cont_perp
        return out


    def _initialise_dose_space(self, n_doses):

        self.contour_perp = np.linspace(0, 1, n_doses+2)[1:-1]
        self.contours_radial = [2**(n) for n in range(-floor(n_doses/2), floor(n_doses/2)+1, 1)]

        self.f1_vals = np.zeros((n_doses, n_doses))
        self.f2_vals = np.zeros((n_doses, n_doses))


    def _run_dose_space(self, Conf):
        
        self._get_fung_strat(Conf)

        self._initialise_dose_space(Conf.n_doses)     
        
        for i in tqdm(range(Conf.n_doses)):
            for j in range(Conf.n_doses):

                f1_val = fsolve(self.constant_effect, 
                                args=(self.contours_radial[j],
                                      self.contour_perp[i]), x0=0.001)[0]

                f2_val = self.contours_radial[j]*f1_val

                if f1_val>1 or f2_val>1:
                    self.LTY[i,j] = None
                    self.TY[i,j] = None
                    self.FY[i,j] = None
                    continue

                Conf.fung1_doses, Conf.fung2_doses = self.fung_strat.get_doses(
                                                            f1_val, f2_val, None)

                one_tact_output =  self.sing_tact.run_single_tactic(Conf)
                
                self.f1_vals[i,j] = f1_val
                self.f2_vals[i,j] = f2_val
                
                self._post_process_multi(one_tact_output, i, j, Conf)
                
        self.t_vec = one_tact_output['t_vec']

    def _save_dose_space(self):
        ds_output = {'LTY': self.LTY,
                    'TY': self.TY,
                    'FY': self.FY,
                    'f1_vals': self.f1_vals,
                    'f2_vals': self.f2_vals,
                    'contour_perp': self.contour_perp,
                    'contours_radial': self.contours_radial,
                    'yield_array': self.yield_array,
                    'res_arrays': self.res_arrays,
                    'start_freqs': self.start_freqs,
                    'selection_arrays': self.selection_arrays,
                    'inoc_array': self.inoc_array,
                    't_vec': self.t_vec,
                    'econ': self.econ,
                    }

        object_dump(self.filename, ds_output)

        return ds_output



    def run_loop(self, ConfigDS):
        """
        Run across grid
        """
        
        filename = ConfigDS.config_string
        self.filename = filename.replace("grid", "dose_space")

        if ConfigDS.load_saved:
            loaded_run = self._load_multi_tactic(self.filename)
            if loaded_run is not None:
                return loaded_run

        Conf = copy.copy(ConfigDS)

        self._initialise_multi_vars(Conf.n_doses, Conf.n_years)

        self._run_dose_space(Conf)
        
        ds_output = self._save_dose_space()

        return ds_output



class RunRadial(RunMultipleTactics):

    def _run_angle_radius_combos(self, ConfigDS):

        Conf = copy.copy(ConfigDS)

        self._get_fung_strat(Conf)
        
        angles = np.linspace(0, pi/2, Conf.n_angles)
        radii = np.linspace(0, 2**(0.5), Conf.n_radii+1)[1:]

        self.row_list = []

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
                
                Conf.fung1_doses, Conf.fung2_doses = self.fung_strat.get_doses(
                                                        f1_val, f2_val, None)

                one_tact_output =  self.sing_tact.run_single_tactic(Conf)
                
                lty = self._lifetime_yield(one_tact_output['yield_vec'], one_tact_output['failure_year'])
                ty = self._total_yield(one_tact_output['yield_vec'])
                fy = one_tact_output['failure_year']
                
                total_dose_f1 = Conf.fung1_doses['spray_1'][0] + Conf.fung1_doses['spray_2'][0]
                total_dose_f2 = Conf.fung2_doses['spray_1'][0] + Conf.fung2_doses['spray_2'][0]

                econ = self._economic_life(one_tact_output['yield_vec'], total_dose_f1, total_dose_f2)

                self.row_list.append(dict(d1=f1_val,
                            d2=f2_val,
                            LTY=lty,
                            TY=ty,
                            FY=fy,
                            Econ=econ,
                            angle=angle,
                            radius=radius,
                            ))



    def master_loop_radial(self, Conf):
        """
        Run radially in dose space
        """

        filename = Conf.config_string
        filename = filename.replace("grid", "radial")

        if Conf.load_saved:
            loaded_run = self._load_multi_tactic(filename)
            if loaded_run is not None:
                return loaded_run

        self._run_angle_radius_combos(Conf)
                
        df_out = pd.DataFrame(self.row_list)

        object_dump(filename, df_out)

        return df_out



class RunRfRatio(RunSingleTactic):
    def __init__(self, grid_output):
        # super.__init__(None)

        self.FY_array = grid_output['FY'] 
        self.res_arrays = grid_output['res_arrays']
        self.start_freqs = grid_output['start_freqs']

        self.N = self.FY_array.shape[0]

        self._find_opt_FY()


    def _find_opt_FY(self):
        """
        Get failure year for res freq @ breakdown strat
        """
        self.opt_FY_rfb = int(np.amax(self.FY_array))


    def _get_rfb_df(self, f2):
        delt_rf_list = []
        f1_list = []
        f2_list = []
        
        for f1 in range(self.N):
            fy = int(self.FY_array[f1, f2])
            
            rf1 = self.res_arrays['f1'][f1,f2,fy]
            rf2 = self.res_arrays['f2'][f1,f2,fy]
            
            if fy>0:
                delt_rf_list.append(logit10_difference(rf1, rf2))
                f1_list.append(f1)
                f2_list.append(f2)
        
        df = pd.DataFrame(dict(metric=delt_rf_list, f1=f1_list, f2=f2_list))
        return df
    

    def _get_eq_sel_df(self, f2):
        eq_sel_list = []
        f1_list = []
        f2_list = []
        
        for f1 in range(self.N):
            fy = int(self.FY_array[f1, f2])
            
            sr0 = self.start_freqs['SR'][f1,f2,0]
            sr1 = self.start_freqs['SR'][f1,f2,1]
            
            rs0 = self.start_freqs['RS'][f1,f2,0]
            rs1 = self.start_freqs['RS'][f1,f2,1]
            
            s1 = rs1/rs0
            s2 = sr1/sr0

            if fy>0:
                eq_sel_list.append(-1 + s1/s2)
                f1_list.append(f1)
                f2_list.append(f2)
        
        df = pd.DataFrame(dict(metric=eq_sel_list, f1=f1_list, f2=f2_list))
        return df


    @staticmethod
    def _interpolate_dfs(lessthan0, morethan0):
        lower_f1 = list(lessthan0.f1)[-1]
        upper_f1 = list(morethan0.f1)[0]
        
        lower_metric = list(lessthan0.metric)[-1]
        upper_metric = list(morethan0.metric)[0]

        f1_out = lower_f1 + (upper_f1-lower_f1) * (-lower_metric) / (upper_metric - lower_metric)
        return f1_out


    def get_f1_from_df(self, df):
        lessthan0 = df[df['metric']<0]
        equals0 = df[df['metric']==0]
        morethan0 = df[df['metric']>0]

        if len(equals0):
            print("is 0")
            if len(equals0.f2)>1:
                raise Exception("shouldn't be more than 1 dose with equal one selection")
            else:
                return equals0.f1

        else:
            if len(lessthan0) and len(morethan0):
                f1_out = self._interpolate_dfs(lessthan0, morethan0)
                return f1_out
            else:
                return None
    @staticmethod
    def _is_decreasing(L):
        return all(x>y for x, y in zip(L, L[1:]))
    
    @staticmethod
    def _is_increasing(L):
        return not all(x>y for x, y in zip(L, L[1:]))

    def get_f1_forK_from_df(self, df, K):
        lessthanK = df[df['metric']<K]
        equalsK = df[df['metric']==K]
        morethanK = df[df['metric']>K]
        
        # check monotonicity
        if not self._is_decreasing(list(df['metric'])):
            # raise Exception(L, [x<y for x, y in zip(L, L[1:])], "non monotone")
            print("not decreasing")
            # return None
        
        if not self._is_increasing(list(df['metric'])):
            # raise Exception(L, [x<y for x, y in zip(L, L[1:])], "non monotone")
            print("not increasing")
            # return None
        
        if self._is_increasing(list(df['metric'])):
            # raise Exception(L, [x<y for x, y in zip(L, L[1:])], "non monotone")
            print("IS increasing")
            # return None


        if len(equalsK):
            print("is K")
            if len(equalsK.f2)>1:
                raise Exception("shouldn't be more than 1 dose")
            else:
                return equalsK.f1

        else:
            if len(lessthanK) and len(morethanK):
                f1_out = self._interpolate_dfs(lessthanK, morethanK)
                if f1_out>30:
                    print(df, K, "wow")
                return f1_out
            else:
                return None


    def _find_RFB_contour(self):
        """
        Find contour along which rfs are equal in breakdown year
        """
        self.df_RFB = pd.DataFrame()
        for f2 in range(self.N):
            df = self._get_rfb_df(f2)
            f1 = self.get_f1_from_df(df)

            if f1 is not None:
                self.df_RFB = self.df_RFB.append(dict(f1_inds=f1, f2_inds=f2), ignore_index=True)
        
        self.df_RFB['f1'] = [i/(self.N-1) for i in self.df_RFB.f1_inds]
        self.df_RFB['f2'] = [i/(self.N-1) for i in self.df_RFB.f2_inds]

    def _find_RFB_contour_at_K(self, K):
        """
        Find contour along which logit rfs differ by K in breakdown year
        """
        df_RFB = pd.DataFrame()
        for f2 in range(self.N):
            df = self._get_rfb_df(f2)
            f1 = self.get_f1_forK_from_df(df, K)

            if f1 is not None:
                df_RFB = df_RFB.append(dict(f1_inds=f1, f2_inds=f2), ignore_index=True)
        
        if "f1_inds" in df_RFB.columns:
            df_RFB['f1'] = [i/(self.N-1) for i in df_RFB.f1_inds]
            df_RFB['f2'] = [i/(self.N-1) for i in df_RFB.f2_inds]
            return df_RFB
        else:
            return None


    def _find_EqS_contour(self):
        """
        Find contour along which rfs are equal in breakdown year
        """
        self.df_EqS = pd.DataFrame()
        for f2 in range(self.N):
            df = self._get_eq_sel_df(f2)
            f1 = self.get_f1_from_df(df)

            if f1 is not None:
                self.df_EqS = self.df_EqS.append(dict(f1_inds=f1, f2_inds=f2), ignore_index=True)
        
        self.df_EqS['f1'] = [i/(self.N-1) for i in self.df_EqS.f1_inds]
        self.df_EqS['f2'] = [i/(self.N-1) for i in self.df_EqS.f2_inds]


    def get_multi_RFB_contours(self, cont_list):
        out = []
        for k in cont_list:
            cont_df = self._find_RFB_contour_at_K(k)
            
            if cont_df is not None:
                x_list = list(cont_df.f1)
                y_list = list(cont_df.f2)
                out.append([x_list, y_list])
        return out


    def get_contours(self):
        self._find_RFB_contour()
        self._find_EqS_contour()
        return [self.df_RFB, self.df_EqS]
    

    def _run_RFB_contour(self, res_props):
        self._find_RFB_contour()

        d1s = self.df_RFB.f1
        d2s = self.df_RFB.f2

        rf1 = res_props['f1']
        rf2 = res_props['f2']

        FY_list = []

        for d1, d2 in zip(d1s, d2s):
            
            d1 = float(0.5*d1)
            d2 = float(0.5*d2)

            ConfigSingleRun = SingleConfig(30, rf1, rf2, d1, d1, d2, d2)
            output = RunSingleTactic().run_single_tactic(ConfigSingleRun)
            FY_list.append(output['failure_year'])
        
        self.opt_RFB = max(FY_list)

    def _run_EqS_contour(self, res_props):
        self._find_EqS_contour()

        d1s = self.df_EqS.f1
        d2s = self.df_EqS.f2

        rf1 = res_props['f1']
        rf2 = res_props['f2']

        FY_list = []

        for d1, d2 in zip(d1s, d2s):
            
            d1 = float(0.5*d1)
            d2 = float(0.5*d2)

            ConfigSingleRun = SingleConfig(30, rf1, rf2, d1, d1, d2, d2)
            output = RunSingleTactic().run_single_tactic(ConfigSingleRun)
            FY_list.append(output['failure_year'])
        
        self.opt_EqS = max(FY_list)



    def run_contours(self, res_props):
        self._run_RFB_contour(res_props)
        self._run_EqS_contour(res_props)

        return [self.opt_RFB, self.opt_EqS]
        

        



def compare_me_vs_hobb_by_ratio(ConfigGridRun, rf1, ratios):
    EqS_list = []
    RFB_list = []

    for ratio in ratios:
        rf2 = ratio*rf1 
        ConfigGridRun.res_props = dict(f1=rf1, f2=rf2)
        ConfigGridRun.add_string()
        grid = RunGrid().grid_of_tactics(ConfigGridRun)
        output = RunRfRatio(grid).run_contours(ConfigGridRun.res_props)
        RFB_list.append(output[0])
        EqS_list.append(output[1])

    df = pd.DataFrame(dict(ratio=ratios, EqS=EqS_list, RFB=RFB_list))
    return df

# End of RunMultipleTactics classes



class EqualResFreqBreakdownArray:
    def __init__(self, grid_output) -> None:
        self.FYs = grid_output['FY']
        self.res_arrays = grid_output['res_arrays']
        self.array = self._generate_RFB_array()
        self.is_valid = self._check_valid()
        
    def _generate_RFB_array(self):
        FYs = self.FYs

        out = np.ones(FYs.shape)
        
        for i, j in itertools.product(range(FYs.shape[0]), range(FYs.shape[1])):
            
            fy = FYs[i,j]

            if not fy>0:
                out[i,j] = None
            
            else:
                
                rf1 = self.res_arrays['f1'][i,j,int(fy)]
                rf2 = self.res_arrays['f2'][i,j,int(fy)]

                try:
                    out[i,j] = logit10_difference(rf1, rf2)
                except:
                    out[i,j] = None

        return out            


    def _check_valid(self):
        """
        Are there some doses where RFB favours fung A and 
        some where RFB favours fung B?
        """
        return (np.nanmax(self.array)>0
                            and np.nanmin(self.array)<0)




class EqualSelectionArray:
    def __init__(self, grid_output) -> None:
        
        self.FYs = grid_output['FY']

        self.start_freqs = grid_output['start_freqs']

        self.array = self._generate_EqSel_array()

        self.is_valid = self._check_valid()
        




    def _generate_EqSel_array(self):

        start_freqs = self.start_freqs
        
        out = np.ones(start_freqs['SR'][:,:,0].shape)
        
        for i, j in itertools.product(range(out.shape[0]), 
                                        range(out.shape[1])):
            
            fy = self.FYs[i,j]

            if not fy>0:
                out[i,j] = None
            
            else:
                
                sr1 = start_freqs['RS'][i,j,1]/start_freqs['RS'][i,j,0]
                sr2 = start_freqs['SR'][i,j,1]/start_freqs['SR'][i,j,0]

                try:
                    out[i,j] = sr1/(sr1+sr2)
                except:
                    out[i,j] = None

        return out


    def _check_valid(self):
        """
        Check if Equal Selection is a possible tactic:
        
        - are there dose pairs for which can select
        more strongly for either fcide?
        """
        return (np.nanmax(self.array)>0.5
                            and np.nanmin(self.array)<0.5)


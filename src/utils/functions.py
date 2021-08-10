import itertools
import numpy as np
from math import exp, ceil, floor, log10
from scipy.integrate import simps, ode
import pickle
import os
import copy
from tqdm import tqdm

from .params import PARAMS
from runHRHR.config_classes import SingleConfig

# * TOC
# Utility functions
# Changing doses fns
# cls Simulator
# cls RunSingleTactic
# cls RunGrid
# and other auxiliary classes


#----------------------------------------------------------------------------------------------
# Utility functions

def object_dump(file_name, object_to_dump):
    
    # check if file path exists - if not create
    outdir =  os.path.dirname(file_name)
    if not os.path.exists(outdir):
        os.makedirs(outdir,exist_ok=True) 
        
    with open(file_name, 'wb') as handle:
        pickle.dump(object_to_dump, handle, protocol=pickle.HIGHEST_PROTOCOL) # protocol?





def logit10(x):
    if x>0 and x<1:
        return log10(x/(1-x))
    else:
        raise Exception(f"x={x} - invalid value")


def logit10_difference(x1, x2):
    return logit10(x1) - logit10(x2)

def log10_difference(x1, x2):
    return log10(x1) - log10(x2)

# End of utility functions


# * changing doses fns

def get_SR_by_doses(doses, freqs):
    outputs = {}
    for dose, rf in itertools.product(doses, freqs):
        ConfigSingleRun = SingleConfig(1, rf, rf, dose, dose, dose, dose)
        output = RunSingleTactic().run(ConfigSingleRun)
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



#----------------------------------------------------------------------------------------------

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
    



    def run(self, fung1_doses, fung2_doses, primary_inoc):
        self.primary_inoc = primary_inoc
        self.fung1_doses = fung1_doses
        self.fung2_doses = fung2_doses

        self._solve_ode()
        
        out = self._post_process()

        return out


    def run_disease_free(self):

        y0 = [PARAMS.S_0] + [0]*(PARAMS.no_variables-1)

        sol = ode(self.ode_system).set_integrator('dopri5', nsteps=PARAMS.nstepz)
        
        tim_1, t_yield = self._get_dis_free_t_vecs()

        _ = self._get_y_array_this_segment(y0, sol, tim_1)
        y_yield = self._get_y_array_this_segment(sol.y, sol, t_yield)
        y_yield[:, -1] = sol.y

        dis_free_yield = simps(y_yield[0,:], t_yield)

        return dis_free_yield

    @staticmethod
    def _get_dis_free_t_vecs():
        t0 = PARAMS.T_emerge
        t1 = PARAMS.T_GS61
        t2 = PARAMS.T_GS87
        
        n1 = 1 + (t1-t0)/PARAMS.dt
        n2 = 1 + (t2-t1)/PARAMS.dt
        
        c1 = ceil(n1-0.5)
        c2 = ceil(n2-0.5)
        
        tim_1 = np.linspace(t0,t1,c1)
        t_yield = np.linspace(t1,t2,c2)
        
        return tim_1, t_yield



    def _solve_ode(self):

        sol = ode(self.ode_system).set_integrator('dopri5', nsteps=PARAMS.nstepz)

        time_list = [PARAMS.T_emerge,
            PARAMS.T_GS32,
            PARAMS.T_GS39,
            PARAMS.T_GS61,
            PARAMS.T_GS87]
        
        y_list, t_list = self._get_y_and_t_lists(time_list, sol)
        
        self.soln_t = np.concatenate(t_list)
        
        solutionTranspose  = np.concatenate(y_list, axis=1)

        self.soln  = np.transpose(solutionTranspose)

        self.yield_integral = YieldFinder(y_list[-1], t_list[-1]).yield_




    def _get_y_and_t_lists(self, time_list, sol):
        """
        Solves the ODE in 4 stages:
        
        - before first spray
        - spray 1
        - spray 2
        - yield period (end of season)

        """
        y0_new = None
        
        y_list = []
        t_list = []

        segments = ["start", "spray_1", "spray_2", "end"]

        list_of_tvs = self._get_list_of_time_vecs(time_list, segments)

        for time_vec, segment in zip(list_of_tvs, segments):

            y0_new = self._get_y0_this_segment(segment, sol)

            y_array = self._get_y_array_this_segment(y0_new, sol, time_vec)

            if segment=="end":
                y_out, t_out = self._update_lists_final_seg(y_array, time_vec, sol)
            else:
                y_out, t_out = self._update_lists_other_segs(y_array, time_vec)

            y_list.append(y_out)
            t_list.append(t_out)

        return y_list, t_list




    @staticmethod
    def _get_list_of_time_vecs(time_list, segments):
        sum_ns = 0

        list_of_tvs = []

        for ii, segment in enumerate(segments):

            if segment=="end":
                # makes sure total number of points is PARAMS.t_points
                n = 3 + PARAMS.t_points - sum_ns
            else:
                n = 1 + (time_list[ii+1]-time_list[ii])/PARAMS.dt
                n = floor(n)
                sum_ns += n

            time_vec = np.linspace(time_list[ii], time_list[ii+1], n)

            list_of_tvs.append(time_vec)

        return list_of_tvs


    def _get_y0_this_segment(self, segment, sol):
        
        if segment=="start":
            PI = self.primary_inoc

            y0_new = [PARAMS.S_0] + [0]*9 + [PI['RR'], 
                        PI['RS'],
                        PI['SR'],
                        PI['SS']] + [0]*2

        else:
            y0_new = sol.y

            if segment in ["spray_1", "spray_2"]:
                y0_new[PARAMS.Fung1_ind] = y0_new[PARAMS.Fung1_ind] + self.fung1_doses[segment]
                y0_new[PARAMS.Fung2_ind] = y0_new[PARAMS.Fung2_ind] + self.fung2_doses[segment]
        
        return y0_new
    
    
    
    @staticmethod
    def _get_y_array_this_segment(y0_new, sol, time_vec):
        
        sol.set_initial_value(y0_new, time_vec[0])

        y_array  = np.zeros((PARAMS.no_variables, len(time_vec)))

        for index, t in enumerate(time_vec[1:]):
            if sol.successful():
                y_array[:,index] = sol.y
                sol.integrate(t)
            else:
                raise RuntimeError('ode solver unsuccessful')
        
        return y_array



    @staticmethod
    def _update_lists_other_segs(y_array, time_vec):
        y_out = y_array[:,:-1]
        t_out = time_vec[:-1]

        return y_out, t_out



    @staticmethod
    def _update_lists_final_seg(y_array, time_vec, sol):
        # final one of loop - need to add final time
        # rather than leave it for start condition 
        # of next run through loop
        y_array[:,-1] = sol.y

        y_out = y_array
        t_out = time_vec
        
        return y_out, t_out




    def _post_process(self):

        final_res_dict, props_out = ResPropFinder().calculate(self.soln)
        
        inoc = PARAMS.init_den

        selection = SelectionFinder(self.primary_inoc, final_res_dict).sel

        out = dict(final_res_dict=final_res_dict,
                props_out=props_out,
                inoc=inoc,
                selection=selection,
                yield_integral=self.yield_integral,
                solution=self.soln,
                solutiont=self.soln_t)
        
        return out






    def ode_system(self, t, y):

        S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,conc_1,conc_2 = y

        A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R

        dydt = [self._growth(A,t)
             - (self._senescence(t))*S
             -  S * (PARAMS.beta/A) * (
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
            
            - self.fcide1.delta * conc_1,
            - self.fcide2.delta * conc_2
            ]

        return dydt

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

# * End of Sim cls


class Fungicide:
    def __init__(self, omega, theta, delta):
        self.omega = omega
        self.theta = theta
        self.delta = delta

    def effect(self, conc):
        effect = 1 - self.omega*(1 - exp(- self.theta * conc))
        return effect

# * End of Fcide cls




class ResPropFinder:
    def __init__(self) -> None:
        pass
    
    def calculate(self, solution):
        """
        Uses final value (end of season) to determine the res props. 

        These are used for next season (with a SR step in between if sr_prop=/=0)
        """

        disease = (solution[-1,PARAMS.IR_ind] + 
                        solution[-1,PARAMS.IRS_ind] +
                        solution[-1,PARAMS.ISR_ind] + 
                        solution[-1,PARAMS.IS_ind])
            
        Res_disease_1 = solution[-1,PARAMS.IR_ind] + solution[-1,PARAMS.IRS_ind]
        Res_disease_2 = solution[-1,PARAMS.IR_ind] + solution[-1,PARAMS.ISR_ind]

        res_props_out = dict(
            f1 = Res_disease_1/disease,
            f2 = Res_disease_2/disease,
            )
        
        strain_frequencies = dict(
            RR = solution[-1,PARAMS.IR_ind]/disease,
            RS = solution[-1,PARAMS.IRS_ind]/disease,
            SR = solution[-1,PARAMS.ISR_ind]/disease,
            SS = solution[-1,PARAMS.IS_ind]/disease
            )
        
        return res_props_out, strain_frequencies

# * End of RPFinder cls
    

class SelectionFinder:    
    def __init__(self, primary_inoc, final_res_dict) -> None:
        self._get_init_res_dict(primary_inoc)

        self.final_res_dict = final_res_dict

        self._get_selection()


    def _get_init_res_dict(self, primary_inoc):
        self.initial_res_dict = dict(
            f1 = primary_inoc['RR'] + primary_inoc['RS'],
            f2 = primary_inoc['RR'] + primary_inoc['SR']
            ) 


    def _get_selection(self):
        self.sel = dict(f1=1, f2=1)

        for key in ['f1','f2']:
            self._get_sel_this_fung(key)
        


    def _get_sel_this_fung(self, key):
        ird = self.initial_res_dict
        frd = self.final_res_dict
        
        if ird[key] > 0:
            self.sel[key] = frd[key] / (ird[key]/PARAMS.init_den)

# * End of SelFinder cls





class YieldFinder:
    def __init__(self, y, t) -> None:
        """
        y and t are from final stage of growing season
        """
        self.y = y
        self.t = t
        self.yield_ = self._get_yield()
    
    def _get_yield(self):
        out = simps(self.y[PARAMS.S_ind,:] + 
                            self.y[PARAMS.ER_ind,:] +
                            self.y[PARAMS.ERS_ind,:] +
                            self.y[PARAMS.ESR_ind,:] +
                            self.y[PARAMS.ES_ind,:],
                            self.t)
        return out 

# * End of YldFinder cls



class FungicideStrategy:
    def __init__(self, my_strategy, n_seasons):
        self.my_strategy = my_strategy
        self.n_seasons = n_seasons



    def get_doses(self, f1_val, f2_val, n_doses):

        self._set_concs(f1_val, f2_val, n_doses)

        self._dose_for_this_strategy()      

        return self.fung1_doses, self.fung2_doses


    def _set_concs(self, f1_val, f2_val, n_doses):
        if n_doses is not None:
            self._set_grid_concs(f1_val, f2_val, n_doses)
        else:
            self._set_regular_concs(f1_val, f2_val)


    def _set_grid_concs(self, f1_val, f2_val, n_doses):
        self.conc_f1 = f1_val/(n_doses-1)
        self.conc_f2 = f2_val/(n_doses-1)


    def _set_regular_concs(self, f1_val, f2_val):
        self.conc_f1 = f1_val
        self.conc_f2 = f2_val



    def _dose_for_this_strategy(self):
        if self.my_strategy=='mix':
            self._get_mixed_doses()
                
        elif self.my_strategy=='alt_12':
            self._get_alt_12_doses()
        
        elif self.my_strategy=='alt_21':
            self._get_alt_21_doses()
        
        else:
            raise Exception("incorrect strategy named")




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

# * End of FcideStrt cls






class RunSingleTactic:
    def __init__(self, fcide_parms=None):

        self.sim = Simulator(fcide_parms)

        self.dis_free_yield = self.sim.run_disease_free()

        self.yield_stopper = 95

        self.PATHOGEN_STRAIN_NAMES = ['RR', 'RS', 'SR', 'SS']





    def run(self, Config):
        """
        Run HRHR model for one strategy
        """
        
        self.filename = Config.config_string
        
        if Config.load_saved:
            loaded_run = self._load_single_tactic()
            if loaded_run is not None:
                return loaded_run

        self._initialise_variables_single_run(Config)

        self._loop_over_years(Config)
        
        model_output = self._save_single_run()
        
        return model_output

    




    
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

        self.selection_vec_dict = self._initialise_dict_of_vecs(['f1', 'f2'], self.n_years+1)
        
        # post-sex from previous year
        self.start_of_season = self._initialise_dict_of_vecs(self.PATHOGEN_STRAIN_NAMES, self.n_years+1)

        # pre-sex
        self.end_of_season = self._initialise_dict_of_vecs(self.PATHOGEN_STRAIN_NAMES, self.n_years+1)
        
        self.inoc_vec[0] = PARAMS.init_den

        self.strain_freqs = primary_inoculum





    def _initialise_dict_of_vecs(self, keys, length):
        out = {}
        
        for key in keys:
            out[key] = np.zeros(length)

        return out
    

    def _initialise_res_vec_dict(self, res_props):
        out = {}
        keys = ['f1', 'f2']
        
        for key in keys:
            out[key] = np.zeros(self.n_years+1)
            # set first year
            out[key][0] = res_props[key]

        return out





    def _load_single_tactic(self):
        filename = self.filename
        
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
    


    


    def _run_single_year(self, Config, yr):
        
        fung1_doses = self._get_single_fung_dose(Config.fung1_doses, yr)
        
        fung2_doses = self._get_single_fung_dose(Config.fung2_doses, yr)
        
        self.update_start_of_season(yr)

        model_inoc_in = self.get_initial_freqs(yr)
        
        output = self.sim.run(fung1_doses, fung2_doses, model_inoc_in)

        self._process_single_output(output, Config, yr)

        self._update_failure_year(yr)



    @staticmethod
    def _get_single_fung_dose(dose_vec, yr):
        return dict(spray_1 = dose_vec['spray_1'][yr],
                    spray_2 = dose_vec['spray_2'][yr])



    def update_start_of_season(self, yr):
        for key in self.PATHOGEN_STRAIN_NAMES:
            self.start_of_season[key][yr] = self.strain_freqs[key]
    


    def get_initial_freqs(self, yr):
        out = {}
        for key in self.PATHOGEN_STRAIN_NAMES:
            out[key] = self.inoc_vec[yr]*self.strain_freqs[key]
        return out




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









    
    def _loop_over_years(self, Config):
        for yr in range(self.n_years):
            # stop the solver after we drop below threshold
            if not (yr>0 and self.yield_vec[yr-1]<self.yield_stopper):
                self._run_single_year(Config, yr)
        
        if min(self.yield_vec)>PARAMS.yield_threshold:
            self.failure_year = -1









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





    def _update_selection_vec_dict(self, output, yr):
        for key in ['f1', 'f2']:
            self.selection_vec_dict[key][yr] = output['selection'][key]


    def _update_res_vec_dict(self, output, yr):
        for key in ['f1', 'f2']:
            self.res_vec_dict[key][yr] = output['final_res_dict'][key]


    def _update_end_of_season(self, freqs_out, yr):
        for key in freqs_out.keys():
            self.end_of_season[key][yr] = freqs_out[key]
        








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



# * End of RunSingleTactic






class RunGrid:
    def __init__(self, fcide_parms=None):
        self.sing_tact = RunSingleTactic(fcide_parms)



    def run(self, ConfigG):
        """
        Run across grid
        """

        self.filename = ConfigG.config_string

        if ConfigG.load_saved:
            loaded_run = self._load_multi_tactic(self.filename)
            if loaded_run is not None:
                return loaded_run

        Conf = copy.copy(ConfigG)

        self._initialise_multi_vars(Conf.n_doses, Conf.n_years)

        self._run_the_grid(Conf)

        grid_output = self._save_grid()

        return grid_output



    def _run_the_grid(self, Conf):
        fs = FungicideStrategy(Conf.strategy, Conf.n_years)
        
        for f1_ind in tqdm(range(Conf.n_doses)):
            for f2_ind in range(Conf.n_doses):

                Conf.fung1_doses, Conf.fung2_doses = fs.get_doses(
                                                    f1_ind, f2_ind, Conf.n_doses)

                one_tact_output =  self.sing_tact.run(Conf)
                
                self._post_process_multi(one_tact_output, f1_ind, f2_ind, Conf)

        self.t_vec = one_tact_output['t_vec']




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
        self.end_freqs = self._get_dict_of_zero_arrays(strain_keys, (n_doses, n_doses, n_years+1))
        

    

    def _update_dict_array_this_dose(self, to_update, data, f1_ind, f2_ind, key1):

        for key_ in to_update.keys():
            to_update[key_][f1_ind,f2_ind,:] = data[key1][key_]
        
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
        self.end_freqs = self._update_dict_array_this_dose(copy.copy(self.end_freqs), data_this_dose, f1_ind, f2_ind, "end_of_season")
        self.selection_arrays = self._update_dict_array_this_dose(copy.copy(self.selection_arrays), data_this_dose, f1_ind, f2_ind, "selection_vec_dict")


    









    def _save_grid(self):
        grid_output = {'LTY': self.LTY,
                    'TY': self.TY,
                    'FY': self.FY,
                    'yield_array': self.yield_array,
                    'res_arrays': self.res_arrays,
                    'start_freqs': self.start_freqs,
                    'end_freqs': self.end_freqs,
                    'selection_arrays': self.selection_arrays,
                    'inoc_array': self.inoc_array,
                    't_vec': self.t_vec,
                    'econ': self.econ,
                    }
        
        object_dump(self.filename, grid_output)
        
        return grid_output


# End of RunGrid class






class EqualResFreqBreakdownArray:
    def __init__(self, grid_output) -> None:
        self.FYs = grid_output['FY']
        self.end_freqs = grid_output['end_freqs']
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
                
                r1 = self.end_freqs['RS'][i,j,int(fy)-1]
                r2 = self.end_freqs['SR'][i,j,int(fy)-1]

                try:
                    out[i,j] = logit10_difference(r1, r2)

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



# * End of ERFB cls



class EqualSelectionArray:
    def __init__(self, grid_output) -> None:

        """
        NB have changed so it compares end of season and start of season,
        not start of consecutive seasons. 
        
        This is because sexual reproduction can cause a change that wasn't due 
        to the tactic but was just down to ratios starting away from those expected
        in a perfectly mixed population.
        """


        
        self.FYs = grid_output['FY']

        self.start_freqs = grid_output['start_freqs']
        
        self.end_freqs = grid_output['end_freqs']

        self.array = self._generate_EqSel_array()

        self.is_valid = self._check_valid()
        




    def _generate_EqSel_array(self):

        start_freqs = self.start_freqs
        end_freqs = self.end_freqs
        
        out = np.ones(start_freqs['SR'][:,:,0].shape)
        
        for i, j in itertools.product(range(out.shape[0]), 
                                        range(out.shape[1])):
            
            fy = self.FYs[i,j]

            if not fy>0:
                out[i,j] = None
            
            else:
                
                sr1 = end_freqs['RS'][i,j,0]/start_freqs['RS'][i,j,0]
                sr2 = end_freqs['SR'][i,j,0]/start_freqs['SR'][i,j,0]

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



# * End of ES cls

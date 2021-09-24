import numpy as np
from math import exp, log10, floor
from scipy.integrate import simps
import os
import pickle


from model.params import PARAMS


# * Utility functions

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

# * End of utility functions




# * Simulatr functions

def yield_calculator(y, t):
    out = simps(y, t)
    return out




def res_prop_calculator(solution):
    """
    Uses final value of disease densities (end of season) to determine the res props.

    These are used for next season (with a SR step in between if sr_prop=/=0)
    """

    disease = (solution.IRR[-1] + 
                    solution.IRS[-1] +
                    solution.ISR[-1] + 
                    solution.ISS[-1])
        
    Res_disease_1 = solution.IRR[-1] + solution.IRS[-1]
    Res_disease_2 = solution.IRR[-1] + solution.ISR[-1]

    res_props_out = dict(
        f1 = Res_disease_1/disease,
        f2 = Res_disease_2/disease,
        )
    
    strain_frequencies = dict(
        RR = solution.IRR[-1]/disease,
        RS = solution.IRS[-1]/disease,
        SR = solution.ISR[-1]/disease,
        SS = solution.ISS[-1]/disease
        )
    
    return res_props_out, strain_frequencies




# * End of Simulatr functions



# * ODESystem fns

def growth(A, t):
    if t>=PARAMS.T_emerge:
        grw = PARAMS.r*(PARAMS.k-A)
        return grw
    else:
        return 0


def senescence(t):
    if t>=PARAMS.T_GS61:
        out = 0.005*((t-PARAMS.T_GS61)/(PARAMS.T_GS87-PARAMS.T_GS61)) + 0.1*exp(-0.02*(PARAMS.T_GS87-t))
        return out
    else:
        return 0


# * End of ODESystem fns


    
# * Classes

class Fungicide:
    def __init__(self, omega, theta, delta):
        self.omega = omega
        self.theta = theta
        self.delta = delta

    def effect(self, conc):
        effect = 1 - self.omega*(1 - exp(- self.theta * conc))
        return effect

# * End of Fcide cls










class FungicideStrategy:
    def __init__(self, my_strategy, n_seasons):
        self.my_strategy = my_strategy
        self.n_seasons = n_seasons



    def get_grid_doses(self, f1_val, f2_val, n_doses):

        self.conc_f1 = f1_val/(n_doses-1)
        self.conc_f2 = f2_val/(n_doses-1)

        self._get_doses_for_this_strategy()      

        return self.fung1_doses, self.fung2_doses



    def _get_doses_for_this_strategy(self):
        if self.my_strategy=='mix':
            self._get_mixed_doses()

        elif self.my_strategy=='alt_12':
            self._get_alt_12_doses()

        elif self.my_strategy=='alt_21':
            self._get_alt_21_doses()

        else:
            raise Exception(f"Invalid strategy named: {self.my_strategy}")


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
        in_res_dict = self.initial_res_dict
        fn_res_dict = self.final_res_dict
        
        if in_res_dict[key] > 0:
            self.sel[key] = fn_res_dict[key] / (in_res_dict[key]/PARAMS.init_den)

# * End of SelFinder cls









class ModelTimes:
    def __init__(self, params) -> None:
        self.params = params

        self.seg_times = [self.params.T_emerge,
                        self.params.T_GS32,
                        self.params.T_GS39,
                        self.params.T_GS61,
                        self.params.T_GS87]

        
        self.seg_names = ["start", "spray_1", "spray_2", "yield"]

        self.t_vecs = self._get_list_of_time_vecs()

        self.t = self._get_t()


    def _get_list_of_time_vecs(self):
        
        seg_ts = self.seg_times
        
        sum_ns = 0

        list_of_tvs = []

        for ii, segment in enumerate(self.seg_names):

            if segment=="yield":
                # makes sure total number of points is self.params.t_points
                n = 3 + self.params.t_points - sum_ns

            else:
                # make n so that values are approx self.params.dt apart
                n = 1 + (seg_ts[ii+1]-seg_ts[ii])/self.params.dt
                n = floor(n)
                sum_ns += n

            time_vec = np.linspace(seg_ts[ii], seg_ts[ii+1], n)

            if segment=="yield":
                self.t_yield = time_vec

            list_of_tvs.append(time_vec)

        return list_of_tvs


    def _get_t(self):
        tvs = self.t_vecs

        out = np.concatenate([tvs[ii][:-1] if ii!=3 else tvs[ii]
                                for ii in range(len(tvs))])
        return out



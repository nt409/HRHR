import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from math import log10, exp


# from tqdm import tqdm
# import itertools

# from runHRHR.config_classes import GridConfig, SingleConfig
# from utils.params import PARAMS
# from utils.plotting import dose_grid_heatmap, eq_RFB_contours

from utils.functions import RunSingleTactic





class ContourStratDFGetter:
    def __init__(self, RandPars, Strategy, ContQuantFinder,
                    grid_output, contour_level, n_points, strat_name) -> None:

        self.rand_pars = RandPars

        self.strat_obj = Strategy(grid_output)

        self.cont_quant_finder = ContQuantFinder
        
        self.max_grid_EL = np.amax(grid_output['FY'])

        self.level = contour_level
        
        self.n_points = n_points
        
        self.strat_name = strat_name


        self.df = self.get_cntr_outcome_df()
        
        self.summary = ThisStratSummaryDF(self.df, self.strat_name, self.level, 
                                self.max_grid_EL).df





    def get_cntr_outcome_df(self):

        cntr = self._find_cntr()
        
        return self._get_dict_this_cntr(cntr)




    
    def _find_cntr(self):
        
        if not self.strat_obj.is_valid:
            cntr_out = {}
            return cntr_out
        else:    
            DS_min_max = DoseSumExtremesGetter(self.strat_obj.array, self.level).min_max_DS
            
            cntr_out = ContourDoseFinder(self.rand_pars, self.cont_quant_finder,
                            DS_min_max, self.n_points, self.level).doses_out
            return cntr_out






    def _get_dict_this_cntr(self, contour):
        
        df_this_run = pd.DataFrame()

        if not contour:
            return df_this_run


        for dose1, dose2, level in zip(contour['x'], contour['y'], contour['contVal']): 

            outcome = self._get_outcome_this_dose_combo(dose1, dose2)
            
            data = {"contour_level": level,
                     "dose1": dose1,
                     "dose2": dose2,
                     "DS": dose1 + dose2,
                     **outcome}

            df_this_run = df_this_run.append(data, ignore_index=True)
        
        
        return df_this_run





    def _get_outcome_this_dose_combo(self, dose1, dose2):
        """
        Returns dict with
        - contQuant
        - EL
        """

        this_dose_conf = self.rand_pars._get_single_conf(dose1, dose2)
                
        sing_run = RunSingleTactic(self.rand_pars.fung_parms).run_single_tactic(this_dose_conf)
        
        out = {f"{self.strat_name}_contQuant": self.cont_quant_finder(sing_run).cont_quant,
                f"{self.strat_name}_EL": sing_run['failure_year']}
        # econ?
        
        return out














class DoseSumExtremesGetter:

    def __init__(self, z, level) -> None:

        cs = self._get_contour_set(z, level)
        
        self.min_max_DS = self._get_min_max_dose_sums(cs)


      


    def _get_contour_set(self, z, level):

        x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]

        cs = plt.contour(x, y, z, levels=[level])
        
        return cs





    def _get_min_max_dose_sums(self, cs):

        if not cs.allsegs:
            print("Warning: cs.allsegs was empty!")
            # why not?
            # suspect too many 'nan's to create a proper contour
            return {}

        else:
            cont = cs.allsegs[0][0]
            
            x_vals = np.asarray(cont[:,0])
            y_vals = np.asarray(cont[:,1])

            dose_sum = x_vals+y_vals

            return dict(min=min(dose_sum), max=max(dose_sum))











class ContourDoseFinder:
    def __init__(self, RandPars, ContQuantFinder, DS_min_max, n_points, level, tol=0.001) -> None:
        self.rand_pars = RandPars
        self.cont_quant_finder = ContQuantFinder
        self.DS_min_max = DS_min_max
        self.n_points = n_points
        self.level = level

        self.min_dist_from_contour = tol

        self.doses_out = self.get_doses_on_contour()
    

    def get_doses_on_contour(self):
        
        DS_bds = self.DS_min_max

        dose_sums = np.linspace(DS_bds['min'], DS_bds['max'], self.n_points)

        x_list = []
        y_list = []
        contVal = []

        for ds in dose_sums:
            self.dose_sum = ds

            k = self._get_doses_on_contour_single_DS(ds)

            # only add if have got close to the contour (but why doesn't it get close otherwise?)
            if abs(self.model-self.level) < self.min_dist_from_contour:
                x_list.append(0.5*ds - k)
                y_list.append(0.5*ds + k)
                contVal.append(self.model)
            else:
                print("\n")
                print("this run didn't get close to the contour?? ...")
                print("contour level:", self.model)
                print("dose sum:", self.dose_sum)
                print(self.DS_min_max)
                

        return dict(x=x_list, y=y_list, contVal=contVal)




    # def plot_obj_fn(self):
    #     y_list = []
        
    #     xx = np.linspace(-0.06, -0.04, 100)
    #     for k in xx:
    #         y = self.objective_fn([k])
    #         y_list.append(y)

    #     plt.scatter(x=xx, y=y_list)
    #     plt.show()
        

    
    def _get_doses_on_contour_single_DS(self, ds):
        
        lower, upper = self._get_bnds(ds)
        
        x0 = [0.5*(lower+upper)]
        
        bnds = ((lower, upper), )

        thisFit = minimize(self.objective_fn, x0, bounds=bnds)
        
        return thisFit.x[0]




    def _get_bnds(self, ds):
        d = 0.5*ds

        # pick lower/upper carefully so that each dose in [0,1]
        lower = max(-d, d - 1)
        upper = min( d,  1 - d)

        return lower, upper






    def objective_fn(self, param):

        k = param[0]

        dose1 = 0.5*self.dose_sum - k
        dose2 = 0.5*self.dose_sum + k

        this_dose_conf = self.rand_pars._get_single_conf(dose1, dose2)
        this_dose_conf.load_saved = False

        sing_run = RunSingleTactic(self.rand_pars.fung_parms).run_single_tactic(this_dose_conf)
        
        self.model = self.cont_quant_finder(sing_run).cont_quant

        dist = self._get_dist(k)

        return dist




    def _get_dist(self, k):
        """
        If quantity not defined (e.g. this dose combo is in unacceptable region)
        then return a high distance and one which pushes 
        the optimiser back towards the middle of the x+y=DS line
        """

        try:
            return (self.model - self.level)**2
        except:
            return abs(k) + 5




                








class ThisStratSummaryDF:

    def __init__(self, df, strat_name, level, max_grid_EL) -> None:
        self.df = df        
        self.strat_name = strat_name
        self.level = level
        self.max_grid_EL = max_grid_EL

        self.df = self.get_summary_df(df, strat_name)





    def get_summary_df(self, df, strat_name):
        df = self.df

        if not df.shape[0]:
            return pd.DataFrame()

        lowDS_val = self._get_low_dose_max_EL()
        medDS_val = self._get_med_dose_max_EL()
        highDS_val = self._get_high_dose_max_EL()

        
        min_opt_dist_from_cntr = self._get_min_opt_dist_from_contour()
        

        data= {
            f"{strat_name}_minContEL": min(df[f"{strat_name}_EL"]),
            f"{strat_name}_maxContEL": max(df[f"{strat_name}_EL"]),

            f"{strat_name}_min_opt_dist_from_contour": min_opt_dist_from_cntr,

            f"{strat_name}_minDS": min(df['DS']),
            f"{strat_name}_maxDS": max(df['DS']),
            
            f"{strat_name}_worked": self.max_grid_EL==max(df[f"{strat_name}_EL"]),
            
            f"{strat_name}_lowDoseMaxEL": lowDS_val,
            f"{strat_name}_medDoseMaxEL": medDS_val,
            f"{strat_name}_highDoseMaxEL": highDS_val,
            }

        return pd.DataFrame([data])


    


    
    
    def _get_low_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + 0.2*(max(df['DS']) - min(df['DS']))
        
        filt = df[df['DS'] < DS_thres_low]

        return self._get_max_if_df_non_empty(filt)
    





    def _get_med_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + 0.2*(max(df['DS']) - min(df['DS']))
        DS_thres_high = min(df['DS']) + 0.8*(max(df['DS']) - min(df['DS']))

        filt = df[((df['DS'] >= DS_thres_low) & (df['DS'] <= DS_thres_high))]

        return self._get_max_if_df_non_empty(filt)

    




    def _get_high_dose_max_EL(self):
        df = self.df

        DS_thres_high = min(df['DS']) + 0.8*(max(df['DS']) - min(df['DS']))

        filt = df[df['DS'] > DS_thres_high]
        
        return self._get_max_if_df_non_empty(filt)
        



    def _get_max_if_df_non_empty(self, df):
        if df.shape[0]:
            return max(df[f"{self.strat_name}_EL"])
        else:
            return "NA"




    def _get_min_opt_dist_from_contour(self):
        df = self.df

        opt_df = df[df[f"{self.strat_name}_EL"]==self.max_grid_EL]

        vec = abs(opt_df[f"{self.strat_name}_contQuant"] - self.level)

        if len(vec):
            out = min(vec)
            return out
        else:
            return "NA"
        



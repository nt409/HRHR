import itertools
import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
from math import log, exp
# import matplotlib.pyplot as plt
# from sklearn import linear_model


from utils.functions import RunGrid, logit10_difference
# from utils.plotting import dose_grid_heatmap, eq_RFB_contours

# TOC:

# Utility fns
# ParamScanGrid
# ParamScanRand
# post process fns
# - combine_PS_grid_outputs

from .functions import ParamScan, get_PS_rand_str, get_PS_grid_str



def get_PS_grid_str(config):
    string = ""
    
    for key in config.keys():
        attribute = config[key]

        if key=="SR":
            string += f"SR={attribute[0]}"
            continue


        if type(attribute)==float or type(attribute)==int:
            string += f"{str(key)[:2]}={attribute}"
        else:
            string += f"{str(key)[:2]}={len(attribute)},{min(attribute)},{max(attribute)}"

    string = string.replace(".", ",")

    return string






class ParamScanGrid(ParamScan):
    """
    almost certainly won't work now have updated ParamScan and moved on
    """

    def __init__(self, config):
        super().__init__()
        self.RFS1 = config["RFS1"]
        self.RFS2 = config["RFS2"]
        self.RFD = config["RFD"]
        self.asym1 = config["asym1"]
        self.asym2 = config["asym2"]
        self.dec_rate1 = config["dec_rate1"]
        self.dec_rate2 = config["dec_rate2"]
        self.SR = config["SR"]
        self.NDoses = config["NDoses"]
        self.load_saved = config["load_saved"]


    def _get_this_run_params_grid(self, sr_prop, run_index):

        sr_dict = {"sr_prop": sr_prop}

        run_params = {**self.fung_parms, **self.inoc, **sr_dict}

        run_params["run"] = f"S{sr_prop}_ind={run_index}"

        self.this_run_params = run_params

    @staticmethod
    def _get_max_ND_RFB(negativeRFB):
        try:
            return max(negativeRFB)
        except:
            return "NA"

    @staticmethod
    def _get_min_PD_RFB(positiveRFB):
        try:
            return min(positiveRFB)
        except:
            return "NA"

    @staticmethod
    def _get_min_AD_RFB(delta):
        try:
            return min([abs(e) for e in delta])
        except:
            return "NA"


    def _add_doses_to_PS_df(self):
        
        out = pd.DataFrame()

        for i in range(self.N):
            this_dose = {}
            for key in self.processed.keys():
                this_dose[key] = self.processed[key][i]

            new_row = {**self.calculated_cols, **self.this_run_params_dict, **this_dose}
            
            out = out.append(new_row, ignore_index=True)

        return out



    def _generate_PS_df(self):
        
        self.N = len(self.processed['EL'])

        
        # if not self.N:
        #     # print("self N:", self.N)
        #     self.df_this_run = None
        #     return None

        deltRFB_list = list(self.processed['delta_RFB'])
        
        positiveRFB = [e for e in deltRFB_list if e>=0]
        
        negativeRFB = [e for e in deltRFB_list if e<0]
        
        self.calculated_cols = dict(
                    
                    maxEL=max(self.processed['EL']),
                    minMS=min(self.processed['MS']),
                    maxMS=max(self.processed['MS']),
                    maxEcon=max(self.processed['econ']),
                    
                    minPosDeltaRFB = self._get_min_PD_RFB(positiveRFB),
                    maxNegDeltaRFB = self._get_max_ND_RFB(negativeRFB),
                    minAbsDeltaRFB = self._get_min_AD_RFB(deltRFB_list),
                    )

        self.df_this_run = self._add_doses_to_PS_df()





    def _get_df_this_run(self):

        self.grid_output = RunGrid(self.fung_parms).grid_of_tactics(self.GridConf)
        
        self._get_dict_of_output()

        self._generate_PS_df()




    def _get_dict_of_output(self):
        
        res_arrays_1 = self.grid_output['res_arrays']['f1']
        res_arrays_2 = self.grid_output['res_arrays']['f2']
        
        FY = self.grid_output['FY']

        econ_array = self.grid_output['econ']

        delt_rf_list = []
        MS_list = []
        f1_list = []
        f2_list = []
        f1_vals = []
        f2_vals = []
        EL_list = []
        econ_list = []
        
        for f1, f2 in itertools.product(range(FY.shape[0]), range(FY.shape[1])):
            fy = int(FY[f1, f2])

            f1 = int(f1)
            f2 = int(f2)
            
            rf1 = res_arrays_1[f1,f2,fy]
            rf2 = res_arrays_2[f1,f2,fy]

            if fy>0:
                try:
                    delt_rf_list.append(logit10_difference(rf1, rf2))
                except:
                    delt_rf_list.append("NA")
                
                f1_list.append(f1)
                f2_list.append(f2)
                
                f1_vals.append(f1/(FY.shape[0]-1))
                f2_vals.append(f2/(FY.shape[1]-1))
                
                MS_list.append(f1/int(FY.shape[0]-1) + f2/int(FY.shape[1]-1))
                EL_list.append(fy)

                econ_list.append(econ_array[f1,f2])
        

        self.processed = dict(delta_RFB=delt_rf_list,
                        f1=f1_list,
                        f2=f2_list,
                        f1_vals=f1_vals,
                        f2_vals=f2_vals,
                        MS=MS_list,
                        EL=EL_list,
                        econ=econ_list)




    def run_param_scan_grid(self):

        df = pd.DataFrame()
        
        run_index = 0

        for rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop in tqdm(itertools.product(
                                self.RFS1,
                                self.RFS2,
                                self.RFD,
                                self.asym1,
                                self.asym2,
                                self.dec_rate1,
                                self.dec_rate2,
                                self.SR)):
            
            self._get_inoc(rfs1, rfs2, rfD)
            
            self._get_fung_parms(om_1, om_2, delt_1, delt_2)

            self.GridConf = self._get_grid_conf(rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop, self.load_saved)

            run_index += 1

            self._get_this_run_params_grid(sr_prop, run_index)

            self._get_df_this_run()

            df = df.append(self.df_this_run, ignore_index=True)
        
        return df



    def run(self, index):
        config_use = copy.copy(self.config)

        config_use["SR"] = [self.config["SR"][index]]

        df = self.run_param_scan_grid()
        
        par_str = get_PS_grid_str(config_use)

        df.to_csv(f"./param_scan_cluster/outputs/param_grid_scan_{par_str}.csv", index=False)





















class ParamScanRandFYY(ParamScanRand):
    def __init__(self, config) -> None:
        super().__init__(config)




    def _find_contours_FYY(self):
        """
        Find doses along the contours of constant first year yield.
        """

        output = RunGrid(self.fung_parms).grid_of_tactics(self.FirstYearConf)

        fyy_array = output['yield_array'][:,:,0]

        levels = self._get_contour_levels_FYY(fyy_array)

        self.contours = _get_contours(fyy_array, levels, step=4)














    @staticmethod
    def _get_contour_levels_FYY(z, num=5):
        """
        takes an array z and returns contour levels from 
        95 to maximum value, logarithmically spaced
        """
        
        max_val = z[-1,-1] - 0.05

        # want min to be just above 95, but has to be less than np.amax(z)
        min_val = min(95.01, (95+np.amax(z))/2)

        if max_val<min_val:
            print("max was less than min:", min_val, max_val)
            max_val = (min_val+z[-1,-1])/2

        exp_levels = np.linspace(exp(min_val), exp(max_val), num=num)

        levels = [log(ii) for ii in exp_levels]

        return levels
    
    







    def run_param_scan_FYY(self, seed):

        df = pd.DataFrame()

        np.random.seed(seed)

        
        for run_index in tqdm(range(self.config["NIts"])):
            
            this_run_parms = self._get_random_pars()

            rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop = this_run_parms

            self._get_inoc(rfs1, rfs2, rfD)
            
            self._get_fung_parms(om_1, om_2, delt_1, delt_2)

            self.FirstYearConf = self._get_grid_conf(*this_run_parms, self.config["load_saved"],
                            n_years=1, n_doses=self.config["grid_number"])
            
            self._find_contours_FYY()

            self._get_this_run_params_rand(sr_prop, run_index)

            self._get_data_this_conf_and_contrs(*this_run_parms)
            
            df = df.append(self.df_this_run, ignore_index=True)
        

        return df











    def run(self, seed):
        """
        Run random scan over uniform dists
        """

        df = self.run_param_scan_FYY(seed)
        
        par_str = get_PS_rand_str(self.config)

        filename = f"./param_scan_cluster/outputs/rand/par_scan/seed={seed}_{par_str}.csv"
        
        print(f"Random Scan, saved as {filename}")
        
        df.to_csv(filename, index=False)












def combine_PS_grid_outputs(config):
    config_use = copy.copy(config)

    df = pd.DataFrame()

    for sr_prop in config["SR"]:
        config_use["SR"] = [sr_prop]
        
        par_str = get_PS_grid_str(config_use)

        temporary = pd.read_csv(f"./param_scan_cluster/outputs/param_grid_scan_{par_str}.csv")
        
        df = df.append(temporary, ignore_index=True)
    
    par_str = get_PS_grid_str(config)
    
    df.to_csv(f"./param_scan_cluster/outputs/PS_combined_grid_output_{par_str}.csv")

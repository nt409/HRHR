import itertools
from numpy.core.fromnumeric import sort
import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
from math import log10, ceil


from runHRHR.config_classes import GridConfig, SingleConfig
from utils.params import PARAMS
from utils.functions import non_decreasing, non_increasing, RunGrid, \
    RunSingleTactic, logit10_difference
from utils.plotting import dose_grid_heatmap, first_year_yield

# TOC:

# Utility fns
# ParamScan
# ParamScanGrid
# ParamScanRand
# post process fns
# - combine_PS_grid_outputs
# - combine_PS_rand_outputs
# PostProcess

#----------------------------------------------------------------------

# Utility fns

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




def get_PS_rand_str(config):
    string = ""
    
    for key in config.keys():
        attribute = config[key]

        if key=="load_saved":
            continue

        if type(attribute)==float or type(attribute)==int:
            string += f"{str(key)[:2]}={attribute}"
        else:
            string += f"{str(key)[:2]}=L{attribute[0]},U{attribute[1]}"

    string = string.replace(".", ",")

    return string



def _get_updated_conf_str(conf_str):
    saved_run_str = conf_str.replace(".png", ".pickle")
    return saved_run_str.replace("figures/", "saved_runs/")



def _get_param_scan_conf_str(rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, Conf):
    conf_str = Conf.config_string_img
    conf_str = conf_str.replace("grid", "param_scan")
    par_str = f"_fung_pars={om_1},{om_2},{delt_1},{delt_2},rf1={rfs1},rf2={rfs2},rfd={rfD}"
    par_str = par_str.replace(".", ",")
    conf_str = conf_str.replace(".png", par_str + ".png")
    
    return conf_str

# End of Utility fns

class ParamScan:
    def __init__(self):
        self.n_years = 20

    def _get_inoc(self, rfs1, rfs2, rfD):
        self.inoc = dict(
                    RR = rfD,
                    RS = rfs1,
                    SR = rfs2,
                    SS = 1 - rfD - rfs1 - rfs2
                    )

    

    def _get_fung_parms(self, omega_1, omega_2, delta_1, delta_2):
        self.fung_parms = dict(
            omega_1 = omega_1,
            omega_2 = omega_2,
            theta_1 = PARAMS.theta_1,
            theta_2 = PARAMS.theta_2,
            delta_1 = delta_1,
            delta_2 = delta_2,
            )       


    def _get_conf(self, rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop, load_saved):

        Conf = GridConfig(self.n_years, None, None, self.NDoses, primary_inoculum=self.inoc)
        
        Conf.sex_prop = sr_prop

        Conf.load_saved = load_saved

        Conf.zeroth_season_reproduction = False

        Conf.add_string()

        conf_str = _get_param_scan_conf_str(rfs1, rfs2, rfD, om_1, om_2,
                                                        delt_1, delt_2, Conf)
        
        Conf.config_string_img = conf_str

        Conf.config_string = _get_updated_conf_str(conf_str)

        self.GridConf = Conf





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

            new_row = {**self.calculated_cols, **self.this_run_params, **this_dose}
            
            out = out.append(new_row, ignore_index=True)

        return out



    def _generate_PS_df(self):
        
        self.N = len(self.processed['EL'])

        
        # if not self.N:
        #     # print("self N:", self.N)
        #     self.df_this_run = None
        #     return None

        delta_RFB_list = list(self.processed['delta_RFB'])
        
        positiveRFB = [e for e in delta_RFB_list if e>=0]
        
        negativeRFB = [e for e in delta_RFB_list if e<0]
        
        self.calculated_cols = dict(
                    
                    maxEL=max(self.processed['EL']),
                    minMS=min(self.processed['MS']),
                    maxMS=max(self.processed['MS']),
                    maxEcon=max(self.processed['econ']),
                    
                    minPosDeltaRFB = self._get_min_PD_RFB(positiveRFB),
                    maxNegDeltaRFB = self._get_max_ND_RFB(negativeRFB),
                    minAbsDeltaRFB = self._get_min_AD_RFB(delta_RFB_list),
                    )

        self.df_this_run = self._add_doses_to_PS_df()








class ParamScanGrid(ParamScan):
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



    def run_param_scan(self):

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

            self._get_conf(rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop, self.load_saved)

            run_index += 1

            self._get_this_run_params_grid(sr_prop, run_index)

            self._get_df_this_run()

            df = df.append(self.df_this_run, ignore_index=True)
        
        return df



    def run(self, index):
        config_use = copy.copy(self.config)

        config_use["SR"] = [self.config["SR"][index]]

        df = self.run_param_scan()
        
        par_str = get_PS_grid_str(config_use)

        df.to_csv(f"param_scan_cluster/outputs/param_grid_scan_{par_str}.csv", index=False)




class ParamScanRand(ParamScan):
    def __init__(self, config) -> None:
        super().__init__()
        self.NDoses = config["NDoses"]
        self.config = config
    


    def _check_validity(self, fungicide_params, pathogen_pars):

        self._get_inoc(pathogen_pars['rfs1'],
                    pathogen_pars['rfs2'],
                    pathogen_pars['rfD'])
        
        ConfigSingleRun = SingleConfig(1, None, None, 1, 1, 1, 1, primary_inoculum=self.inoc)
        
        ConfigSingleRun.sex_prop = pathogen_pars['sr_prop']
        ConfigSingleRun.load_saved = False

        full_dose_run = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
        
        full_dose_yield = full_dose_run['yield_vec'][0]
        
        if full_dose_yield>95:
            print("valid params")
            return True
        else:
            print("invalid params")
            return False



    def _get_random_pars(self):
        conf = self.config

        valid_pars = False

        while valid_pars is False:
        
            rfs1_power = np.random.uniform(low=conf["RFS1"][0], high=conf["RFS1"][1])
            rfs2_power = np.random.uniform(low=conf["RFS2"][0], high=conf["RFS2"][1]) 
            rfD_power = np.random.uniform(low=conf["RFD"][0], high=conf["RFD"][1])
            
            rfs1 = 10**(rfs1_power)
            rfs2 = 10**(rfs2_power)
            rfD = 10**(rfD_power)

            om_1 = np.random.uniform(low=conf["asym1"][0], high=conf["asym1"][1])
            om_2 = np.random.uniform(low=conf["asym2"][0], high=conf["asym2"][1])
            
            delt_1 = np.random.uniform(low=conf["dec_rate1"][0], high=conf["dec_rate1"][1])
            delt_2 = np.random.uniform(low=conf["dec_rate2"][0], high=conf["dec_rate2"][1])
            
            sr_prop = np.random.uniform(low=conf["SR"][0], high=conf["SR"][1])

            fungicide_params = dict(
                omega_1 = om_1,
                omega_2 = om_2,
                theta_1 = PARAMS.theta_1,
                theta_2 = PARAMS.theta_2,
                delta_1 = delt_1,
                delta_2 = delt_2,
            )

            pathogen_pars = dict(rfs1=rfs1,
                rfs2=rfs2,
                rfD=rfD, 
                sr_prop=sr_prop)

            valid_pars = self._check_validity(fungicide_params, pathogen_pars)
            


        return rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop


    def _get_this_run_params_rand(self, sr_prop, run_index):

        sr_dict = {"sr_prop": sr_prop}

        run_params = {**self.fung_parms, **self.inoc, **sr_dict}

        run_params["run"] = f"{run_index}"

        self.this_run_params = run_params




    def run_param_scan(self, seed):

        df = pd.DataFrame()

        np.random.seed(seed)

        
        for run_index in tqdm(range(self.config["NIts"])):
            
            rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop = self._get_random_pars()

            self._get_inoc(rfs1, rfs2, rfD)
            
            self._get_fung_parms(om_1, om_2, delt_1, delt_2)

            self._get_conf(rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop, self.config["load_saved"])
            
            self._get_this_run_params_rand(sr_prop, run_index)

            self._get_df_this_run()
            
            df = df.append(self.df_this_run, ignore_index=True)
        
        return df




    def run(self, seed):
        """
        Run random scan over uniform dists
        """

        df = self.run_param_scan(seed)
        
        par_str = get_PS_rand_str(self.config)
        
        print(f"param_scan_cluster/outputs/param_scan_rand_seed={seed}_{par_str}.csv")
        
        df.to_csv(f"param_scan_cluster/outputs/param_scan_rand_seed={seed}_{par_str}.csv", index=False)











# post process fns

def combine_PS_grid_outputs(config):
    config_use = copy.copy(config)

    df = pd.DataFrame()

    for sr_prop in config["SR"]:
        config_use["SR"] = [sr_prop]
        
        par_str = get_PS_grid_str(config_use)

        temporary = pd.read_csv(f"param_scan_cluster/outputs/param_grid_scan_{par_str}.csv")
        
        df = df.append(temporary, ignore_index=True)
    
    par_str = get_PS_grid_str(config)
    
    df.to_csv(f"param_scan_cluster/outputs/PS_combined_grid_output_{par_str}.csv")


def combine_PS_rand_outputs(config, seeds):

    df = pd.DataFrame()
    
    par_str = get_PS_rand_str(config)

    for seed in seeds:

        temporary = pd.read_csv(f"param_scan_cluster/outputs/param_scan_rand_seed={seed}_{par_str}.csv")
        
        df = df.append(temporary, ignore_index=True)
    
    df.to_csv(f"param_scan_cluster/outputs/PS_combined_rand_output_{par_str}.csv")





# End of post process fns






class PostProcess:

    def __init__(self, df, par_str):
        self.df = df
        self.par_str = par_str



    def _get_best_doses(self):
        self.best_doses = self.df.loc[self.df['EL']==self.df['maxEL']]



    def process_best_doses(self):

        self._get_best_doses()

        dfTest = self.best_doses

        self.delt_worked = dfTest[(dfTest["minAbsDeltaRFB"]==dfTest["delta_RFB"]) | (dfTest["minAbsDeltaRFB"]==-dfTest["delta_RFB"]) ]
        self.all_runs = dfTest["run"].unique()
        self.runs_worked = self.delt_worked["run"].unique()
        self.runs_failed = [e for e in self.all_runs if not e in self.runs_worked]

        self._print_runs_best_dose_not_minRFB()
        



    def _print_runs_best_dose_not_minRFB(self):
        """
        If best EL is not for minimum absolute value for delta RFB (but this is a bit arbitrary)
        """
        self.delt_failed = self.df[self.df["run"].isin(self.runs_failed)]

        print("Min RFB does not give the optimal EL:", self.delt_failed.run.unique())



    @staticmethod
    def _max_MS(data):
        return data["EL"].max()



    def check_max_EL_by_MS(self, type):
        grouped_df = self.df.groupby(["run", "MS"])

        grouped = grouped_df.pipe(self._max_MS)
        grouped.to_csv(f"param_scan_cluster/outputs/{type}/maxEL_thisMS/{self.par_str}.csv")


    @staticmethod
    def _monotone_RFB(data):
        sorted_data = data.sort_values(["delta_RFB"])

        dfPos = sorted_data[sorted_data["delta_RFB"]>=0]
        dfNeg = sorted_data[sorted_data["delta_RFB"]<0]

        ELPos = dfPos["EL"]
        ELNeg = dfNeg["EL"]
        
        out = non_increasing(ELPos) and non_decreasing(ELNeg)

        return out

    @staticmethod
    def _EL_by_RFB(data):
        sorted_data = data.sort_values(["delta_RFB"])

        effective_lives = [str(int(e)) for e in sorted_data["EL"]]

        string = ","
        out = string.join(effective_lives)

        return out

    @staticmethod
    def _rounded_RFB(data):
        sorted_data = data.sort_values(["delta_RFB"])

        delt_RFB = [str(round(e, int(3 - ceil(log10(abs(e)))))) for e in sorted_data["delta_RFB"]]

        string = "/"
        out = string.join(delt_RFB)

        return out




    def check_monotone_RFB(self, type):
        grouped_df = self.df.groupby(["run", "MS"])

        grouped_RFB_is_monotone = grouped_df.apply(self._monotone_RFB)
        
        grouped_effective_lives = grouped_df.apply(self._EL_by_RFB)
        
        grouped_RFBs = grouped_df.apply(self._rounded_RFB)
        
        combined = pd.concat([grouped_effective_lives, grouped_RFB_is_monotone, grouped_RFBs], axis=1)

        combined.columns = ["ELs", "monotoneRFB", "deltRFBs"]
        
        combined = combined.reset_index()

        self.monotoneRFB_df = combined
        
        combined.to_csv(f"param_scan_cluster/outputs/{type}/monotoneRFB/{self.par_str}.csv")
        print(f"Monotone RFB works in all cases: {all(grouped_RFB_is_monotone)}")


    def get_params_for_specific_runs(self, which_runs):

        par_df = self.df

        par_df = par_df.drop(["Unnamed: 0", "EL", "MS", "econ",
            "delta_RFB", "f1", "f2", "f1_vals", "f2_vals"], axis=1)

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = par_df[par_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out


    def which_runs_worked(self, print_=False):

        failed = self.monotoneRFB_df[~self.monotoneRFB_df["monotoneRFB"]]
        succeeded = self.monotoneRFB_df[self.monotoneRFB_df["monotoneRFB"]]

        runs_that_failed = failed["run"].unique()
        runs_that_succeeded = succeeded["run"].unique()

        runs_that_succeeded = [e for e in runs_that_succeeded if e not in runs_that_failed]

        self.failed_pars = self.get_params_for_specific_runs(runs_that_failed)
        self.success_pars = self.get_params_for_specific_runs(runs_that_succeeded)

        n_fail = self.failed_pars.shape[0]
        n_success = self.success_pars.shape[0]
        self.failed_pars.to_csv(f"param_scan_cluster/outputs/failed_pars/failed_{n_fail}_out_of_{n_success+n_fail}.csv")



        if print_:
            print("\n")
            print("These runs didn't have monotone EL vs RFB for every MS :(")
            print("\n")
            print(self.failed_pars)

            print("\n")
            print("These runs had monotone EL vs RFB for every MS :)")
            print("\n")
            print(self.success_pars)


    def re_run_failures(self, NDoses, failed_run_indices=[0]):

        for ii in failed_run_indices:

            pars = self.failed_pars.iloc[int(ii),:]
        
            PS = ParamScan()
            
            PS.NDoses = NDoses

            PS._get_inoc(pars["RS"], pars["SR"], pars["RR"])
            
            PS._get_fung_parms(pars["omega_1"], pars["omega_2"], 
                pars["delta_1"], pars["delta_2"])

            PS._get_conf(pars["RS"], pars["SR"], pars["RR"],
                pars["omega_1"], pars["omega_2"], 
                pars["delta_1"], pars["delta_2"],
                pars["sr_prop"], True
                )
            
            output = RunGrid(PS.fung_parms).grid_of_tactics(PS.GridConf)

            conf_str = PS.GridConf.config_string_img
            dose_grid_heatmap(output, PS.GridConf, "FY", conf_str)
            first_year_yield(output, PS.GridConf)

        



# End of Param Scan fns
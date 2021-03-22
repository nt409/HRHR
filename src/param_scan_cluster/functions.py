import itertools
import pandas as pd
from tqdm import tqdm
import copy

from runHRHR.config_classes import GridConfig
from utils.params import PARAMS
from utils.functions import non_decreasing, non_increasing, RunGrid, \
    logit10_difference

# TOC:

# Utility fns
# ParamScan
# combine_PS_outputs
# PostProcess

#----------------------------------------------------------------------

# Utility fns

def get_par_str(config):
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
    def __init__(self, config):
        self.RFS1 = config["RFS1"]
        self.RFS2 = config["RFS2"]
        self.RFD = config["RFD"]
        self.asym1 = config["asym1"]
        self.asym2 = config["asym2"]
        self.dec_rate1 = config["dec_rate1"]
        self.dec_rate2 = config["dec_rate2"]
        self.SR = config["SR"]
        self.NDoses = config["NDoses"]



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
                            delt_1, delt_2, sr_prop)

            run_index += 1

            out = self._get_df_this_run(sr_prop, run_index)
            
            df = df.append(out, ignore_index=True)
        
        return df



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
                            delt_1, delt_2, sr_prop):

        Conf = GridConfig(20, None, None, self.NDoses, primary_inoculum=self.inoc)
        
        Conf.sex_prop = sr_prop

        Conf.add_string()

        conf_str = _get_param_scan_conf_str(rfs1, rfs2, rfD, om_1, om_2,
                                                        delt_1, delt_2, Conf)
        Conf.config_string_img = conf_str

        Conf.config_string = _get_updated_conf_str(conf_str)

        self.GridConf = Conf
        



    def _get_df_this_run(self, sr_prop, run_index):
        

        self.grid_output = RunGrid(self.fung_parms).grid_of_tactics(self.GridConf)
        
        processed = self._get_dict_of_output()

        run_params = {**self.fung_parms, **self.inoc}
        run_params["run"] = f"S{sr_prop}_ind={run_index}"

        return self._generate_PS_df(processed, sr_prop, run_params)





    def _get_dict_of_output(self):
        
        res_arrays_1 = self.grid_output['res_arrays']['f1']
        res_arrays_2 = self.grid_output['res_arrays']['f2']
        
        FY = self.grid_output['FY']

        delt_rf_list = []
        MS_list = []
        f1_list = []
        f2_list = []
        f1_vals = []
        f2_vals = []
        EL_list = []
        
        for f1, f2 in itertools.product(range(FY.shape[0]), range(FY.shape[1])):
            fy = int(FY[f1, f2])
            
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
                
                MS_list.append(f1/(FY.shape[0]-1) + f2/(FY.shape[1]-1))
                EL_list.append(fy)
        

        return dict(delta_RFB=delt_rf_list,
                        f1=f1_list,
                        f2=f2_list,
                        f1_vals=f1_vals,
                        f2_vals=f2_vals,
                        MS=MS_list,
                        EL=EL_list)


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


    def _get_PS_df(self, processed, other_cols, run_params):
        out = pd.DataFrame()

        for i in range(self.N):
            this_dose = {}
            for key in processed.keys():
                this_dose[key] = processed[key][i]

            new_row = {**other_cols, **run_params, **this_dose}
            out = out.append(new_row, ignore_index=True)
        return out



    def _generate_PS_df(self, processed, sr_prop, run_params):
        
        self.N = len(processed['EL'])
        
        if not self.N:
            return None

        delta_RFB_list = list(processed['delta_RFB'])
        
        positiveRFB = [e for e in delta_RFB_list if e>=0]
        
        negativeRFB = [e for e in delta_RFB_list if e<0]
        

        other_cols = dict(sr_prop=sr_prop,            
                    
                    maxEL=max(processed['EL']),
                    minMS=min(processed['MS']),
                    maxMS=max(processed['MS']),
                    
                    minPosDeltaRFB = self._get_min_PD_RFB(positiveRFB),
                    maxNegDeltaRFB = self._get_max_ND_RFB(negativeRFB),
                    minAbsDeltaRFB = self._get_min_AD_RFB(delta_RFB_list),
                    )

        out = self._get_PS_df(processed, other_cols, run_params)

        return out






# used in analyse

def combine_PS_outputs(config):
    config_use = copy.copy(config)

    df = pd.DataFrame()

    for sr_prop in config["SR"]:
        config_use["SR"] = [sr_prop]
        
        par_str = get_par_str(config_use)

        temporary = pd.read_csv(f"../outputs/csvs/param_scan_{par_str}.csv")
        
        df = df.append(temporary, ignore_index=True)
    
    par_str = get_par_str(config)
    
    df.to_csv(f"../outputs/csvs/PS_combined_output_{par_str}.csv")




class PostProcess:

    def __init__(self, df, par_str):
        self.df = df
        self.par_str = par_str

    def process_param_scan_df(self):
        self.best_doses = self.df.loc[self.df['EL']==self.df['maxEL']]



    def process_best_doses(self):
        dfTest = self.best_doses
        self.delt_worked = dfTest[dfTest["minAbsDeltaRFB"]==dfTest["delta_RFB"]]
        self.all_runs = dfTest["run"].unique()
        self.runs_worked = self.delt_worked["run"].unique()
        self.runs_failed = [e for e in self.all_runs if not e in self.runs_worked]
        

    def show_failed_runs(self):
        self.delt_failed = self.df[self.df["run"].isin(self.runs_failed)]


    @staticmethod
    def _max_MS(data):
        return data["EL"].max()



    def check_max_EL_by_MS(self):
        grouped = (self.df.groupby(["run", "MS"]).pipe(self._max_MS))
        grouped.to_csv(f"../outputs/csvs/maxEL_thisMS_{self.par_str}.csv")


    @staticmethod
    def _monotone_RFB(data):
        sorted_data = data.sort_values(["delta_RFB"])

        dfPos = sorted_data[sorted_data["delta_RFB"]>=0]
        dfNeg = sorted_data[sorted_data["delta_RFB"]<0]

        ELPos = dfPos["EL"]
        ELNeg = dfNeg["EL"]
        
        out = non_increasing(ELPos) and non_decreasing(ELNeg)

        return out



    def check_monotone_RFB(self):
        grouped = (self.df.groupby(["run", "MS"]).apply(self._monotone_RFB))
        grouped.to_csv(f"../outputs/csvs/monotoneRFB_{self.par_str}.csv")
        print(f"Monotone RFB works in all cases: {all(grouped)}")
        
# End of Param Scan fns
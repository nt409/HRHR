import itertools
from numpy import random
import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
from math import log10, ceil
import matplotlib.pyplot as plt


from runHRHR.config_classes import GridConfig, SingleConfig
from utils.params import PARAMS
from utils.functions import RunGrid, RunSingleTactic, \
    logit10_difference, EqualResFreqBreakdownArray, \
    EqualSelectionArray
from utils.plotting import dose_grid_heatmap, eq_RFB_contours

# TOC:

# Utility fns
# ParamScanRand
# RandomPars
# ContourFinder
# ConfigsParScan

# post process fns:
# - combine_PS_rand_outputs
# PostProcess
# MaxAlongContourDF

#----------------------------------------------------------------------

# Utility fns




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










# End of Utility fns












class ParamScanRand:
    def __init__(self, config) -> None:
        self.config = config
    





    def run(self, seed):
        """
        Run random scan over uniform dists
        """

        df = self._run_param_scan_RFB(seed)

        self.save_output(df, seed)
        

    
    
    
    def save_output(self, df, seed):
        
        par_str = get_PS_rand_str(self.config)

        filename = f"./param_scan_cluster/outputs/rand/par_scan/seed={seed}_{par_str}.csv"
        
        print(f"Random Scan, saved as:\n {filename}")
        
        df.to_csv(filename, index=False)







    def _run_param_scan_RFB(self, seed):

        df = pd.DataFrame()

        np.random.seed(seed)

        for run_index in tqdm(range(self.config["NIts"])):

            new_df = self._get_this_single_PS_run(run_index)

            df = pd.concat([df, new_df], axis=0)

        return df
    






    def _get_this_single_PS_run(self, run_index):

        # initialise
        self.df_this_run = pd.DataFrame()

        self._get_params_and_confs(run_index)
        
        self._get_grid_output_this_run()

        RFB_cntrs = self._find_contours_RFB()
    
        if self.rfb_obj.is_valid:
            self._get_data_this_conf_and_contrs(RFB_cntrs)
        
        grid_df = self._get_grid_output_df()

        df_RFB_and_EqSel = self._get_EqSel_columns()

        out = self._combine_all_dfs(df_RFB_and_EqSel, grid_df)

        return out






    def _get_params_and_confs(self, run_index):
    
        RP = RandomPars(self.config, run_index)

        RP.find_pars()

        self.rand_pars = RP
        
        self.ps_confs = ConfigsParScan(*self.rand_pars.pars[:-1],
                            sr_prop=self.rand_pars.pars[-1],
                            inoc=self.rand_pars.inoc,
                            load_saved=self.config["load_saved"],
                            n_years=self.config["n_years"])













    def _get_grid_output_this_run(self):
        
        grid_config = self.ps_confs._get_grid_conf(self.config["grid_number"])

        self.my_grid_output = RunGrid(self.rand_pars.fung_parms).grid_of_tactics(grid_config)










    def _find_contours_RFB(self):
        """
        Find doses along the contours of constant first year yield.
        """

        self.rfb_obj = EqualResFreqBreakdownArray(self.my_grid_output)
        
        if not self.rfb_obj.is_valid:            
            cntrs_out = []
            return cntrs_out
        else:
            cntrs_out = ContourFinder(self.rfb_obj.array, levels=[0]).cont_list
            return cntrs_out








    def _get_data_this_conf_and_contrs(self, cntrs):
        
        df_this_run = pd.DataFrame()

        if not cntrs:
            self.df_this_run = df_this_run
            return None


        for contour in cntrs:
            
            for dose1, dose2 in zip(contour['x'], contour['y']):            

                this_contour_dict = self._get_data_these_doses_on_contour(dose1, dose2)
               
                data = {"contour_level": contour['level'],
                            "max_dose_on_contour": contour['max_dose'],
                            "dose1": dose1,
                            "dose2": dose2,
                            **this_contour_dict}

                df_this_run = df_this_run.append(data, ignore_index=True)
        
        self.df_this_run = self._add_extra_RFB_EL_cols_for_this_cont(df_this_run)





    def _get_data_these_doses_on_contour(self, dose1, dose2):

        this_dose_conf = self.ps_confs._get_single_conf(dose1, dose2)
                
        single_run = RunSingleTactic(self.rand_pars.fung_parms).run_single_tactic(this_dose_conf)
        
        out = self._get_this_contour_dict(single_run)
        
        return out








    def _get_this_contour_dict(self, sing_run_output):

        fy_in = sing_run_output['failure_year']

        res_vecs = sing_run_output['res_vec_dict']
        
        fy = int(fy_in)

        rf1 = res_vecs['f1'][fy]
        
        rf2 = res_vecs['f2'][fy]

        try:
            delt_rf = logit10_difference(rf1, rf2)
        except:
            print(f"problem with calculating delta_RFB for: {rf1, rf2}")
            delt_rf = "NA"

        econ = "NA"
        
        out = dict(delta_RFB=delt_rf,
                                EL=fy,
                                econ=econ)
        
        return out







    @staticmethod
    def _add_extra_RFB_EL_cols_for_this_cont(df):
        
        n_empty_rows = df.shape[0]-1

        df['maxContEL'] = [max(df['EL'])] + [""]*(n_empty_rows)

        df['minDeltaRFB'] = [min(df['delta_RFB'])] + [""]*(n_empty_rows)

        df['maxDeltaRFB'] = [max(df['delta_RFB'])] + [""]*(n_empty_rows)

        return df







    def _get_grid_output_df(self):

        FYs = self.my_grid_output['FY']

        minEqDoseELVec = [float(FYs[i, i]) for i in range(FYs.shape[0])]

        data = dict(maxGridEL = np.amax(FYs),
                ERFB_Valid = self.rfb_obj.is_valid,
                GridMinRFB = np.nanmin(self.rfb_obj.array),
                GridMaxRFB = np.nanmax(self.rfb_obj.array),
                minEqDoseEL = self._get_first_non_zero_element(minEqDoseELVec),
                fullDoseEL = FYs[-1, -1],
                corner_00 = FYs[0, 0],
                corner_10 = FYs[-1, 0],
                corner_01 = FYs[0, -1],
                corner_11 = FYs[-1, -1])

        out = pd.DataFrame([data])

        return out








    def _get_EqSel_columns(self):
        """
        adds on two columns to existing dataframe:
        - maxEqSelEL
        - EqSelValid
        """
        
        df_ERFB_data = copy.copy(self.df_this_run)

        EqSel_cntrs = self._find_contours_EqSel()

        self._get_data_this_conf_and_contrs(EqSel_cntrs)

        maxEqSelEL = self._get_maxEL_EqSel()

        max_D = self._get_maxD(EqSel_cntrs)

        df_use = pd.DataFrame([dict(
                maxEqSelEL = maxEqSelEL,
                EqSelValid = self.eq_sel_obj.is_valid,
                max_dose_EL_cont = max_D
                )])
        
        out = pd.concat([df_ERFB_data, df_use], axis=1)

        return out
        






    def _find_contours_EqSel(self):
        
        self.eq_sel_obj = EqualSelectionArray(self.my_grid_output)
        
        if not self.eq_sel_obj.is_valid:
            cntrs_out = []
            return cntrs_out
        else:    
            cntrs_out = ContourFinder(self.eq_sel_obj.array, levels=[0.5]).cont_list
            return cntrs_out







    @staticmethod
    def _get_first_non_zero_element(vec):
        filtered = list(filter(lambda x: x>0, vec))

        if not filtered:
            return "NA"

        return filtered[0]






    def _get_maxEL_EqSel(self):
        
        data = self.df_this_run

        if data.shape[0]>0:
            return max(data['EL'])
        else:
            return "NA"





    @staticmethod
    def _get_maxD(cntrs):
        if not cntrs:
            return "NA"
        else:
            return cntrs[0]["max_dose"]











    def _combine_all_dfs(self, df_RFB_and_EqSel, grid_df):
        n_rows = df_RFB_and_EqSel.shape[0]
        
        par_df = self._get_params_and_run_index_df(n_rows)

        out = pd.concat([par_df,
                            df_RFB_and_EqSel,
                            grid_df],
                            axis=1)
        
        return out






    def _get_params_and_run_index_df(self, n_rows):

        parms_df = pd.DataFrame([{**self.rand_pars.this_run_params_dict}])

        run_index = int(parms_df.loc[0, 'run'])

        parms_df = parms_df.drop(["run"], axis=1)

        run_df = pd.DataFrame([{"run": run_index}]*n_rows)

        out = pd.concat([parms_df, run_df], axis=1)

        return out















class RandomPars:
    def __init__(self, config, run_index) -> None:

        self.config = config

        self.run_index = run_index



    def find_pars(self):
        """
        returns rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop
        """

        self.pars = self._get_random_pars()

        rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop = self.pars

        self._get_inoc_dict(rfs1, rfs2, rfD)
        
        self._get_fung_parms_dict(om_1, om_2, delt_1, delt_2)
      
        self._get_this_run_all_parms_dict(sr_prop)





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
                delta_2 = delt_2)

            pathogen_pars = dict(rfs1=rfs1,
                rfs2=rfs2,
                rfD=rfD, 
                sr_prop=sr_prop)

            valid_pars = self._check_validity(fungicide_params, pathogen_pars)
            
        return rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop





    def _check_validity(self, fungicide_params, pathogen_pars):

        self._get_inoc_dict(pathogen_pars['rfs1'],
                    pathogen_pars['rfs2'],
                    pathogen_pars['rfD'])
        
        for ii, jj in [[1,0], [0,1], [1,1]]:
            this_dose_yield = self._check_yield(fungicide_params, pathogen_pars, ii, jj)

            if not this_dose_yield>95:
                print("\n")
                print(f"invalid params; {round(this_dose_yield,2)}<=95, dose pair: {ii,jj}")
                return False
        
        print("\n")
        print(f"valid params; {round(this_dose_yield,2)}>95, dose pair: {ii,jj}")
        print("\n")
        return True





    def _check_yield(self, fungicide_params, pathogen_pars, d1, d2):

        ConfigSingleRun = SingleConfig(1, None, None, d1, d1, d2, d2, primary_inoculum=self.inoc)
        
        ConfigSingleRun.sex_prop = pathogen_pars['sr_prop']

        ConfigSingleRun.load_saved = False

        this_run = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
        
        yield_out = this_run['yield_vec'][0]

        return yield_out







    def _get_inoc_dict(self, rfs1, rfs2, rfD):
        self.inoc = dict(
                    RR = rfD,
                    RS = rfs1,
                    SR = rfs2,
                    SS = 1 - rfD - rfs1 - rfs2
                    )

    

    def _get_fung_parms_dict(self, omega_1, omega_2, delta_1, delta_2):
        self.fung_parms = dict(
            omega_1 = omega_1,
            omega_2 = omega_2,
            theta_1 = PARAMS.theta_1,
            theta_2 = PARAMS.theta_2,
            delta_1 = delta_1,
            delta_2 = delta_2,
            )     





    def _get_this_run_all_parms_dict(self, sr_prop):
        self.this_run_params_dict = {**self.fung_parms,
                **self.inoc,
                "sr_prop": sr_prop,
                "run": self.run_index}








class ContourFinder:
    def __init__(self, z, levels, min_n_pts_alng_cntr=12) -> None:
        
        cs = self._get_contour_set(z, levels)
        
        cntrs_out = self._get_contours(cs, levels)

        if cntrs_out and len(cntrs_out[0]['x'])<min_n_pts_alng_cntr:
            cntrs_out = self._interpolate_contours(cntrs_out, min_n_pts_alng_cntr)

        self.cont_list = cntrs_out



    @staticmethod
    def _get_contour_set(z, levels):

        x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]

        cs = plt.contour(x, y, z, levels=levels)

        return cs



    def _get_contours(self, cs, levels):
        """
        Takes a contour set and returns a list of dictionaries containing:
        - x values
        - y values
        - the contour level
        """

        output = []

        for level, conts in zip(levels, cs.allsegs):

            if not conts:
                # why not?
                # suspect too many 'nan's to create a proper contour
                continue
            
            this_cont_dict = self._get_contour_dict(level, conts)

            output.append(this_cont_dict)
        
        return output


    @staticmethod
    def _get_contour_dict(level, conts):
        
        cont = conts[0]
        
        # was x_list, y_list
        x_vals = cont[:,0]
        y_vals = cont[:,1]

        # x_vals = [x_list[0]] + list(x_list[1:-2]) + [x_list[-1]]
        # y_vals = [y_list[0]] + list(y_list[1:-2]) + [y_list[-1]]

        x_max = max(x_vals)
        y_max = max(y_vals)

        print(f"max dose along contour: (x,y) = {round(x_max,2), round(y_max,2)}")

        out = dict(
                x = x_vals,
                y = y_vals,
                max_dose = max(x_max, y_max),
                level = level
                )

        return out

    


    def _interpolate_contours(self, cntrs, num=15):
        """
        Takes old contours and adds in some interpolated values
        
        Result is more densely populated list of values (approximately) 
        along the contour
        """

        xx = copy.copy(cntrs[0]['x'])
        yy = copy.copy(cntrs[0]['y'])

        cntrs[0]['x'] = self._interp_vector(xx, num)
        cntrs[0]['y'] = self._interp_vector(yy, num)
        
        return cntrs
    



    @staticmethod
    def _interp_vector(old, num):
        
        nn = len(old)
        
        fp = list(range(nn))
        
        start = float(fp[0])
        stop = float(fp[-1])

        to_find = np.linspace(start, stop, num)
        
        interp = np.interp(to_find, fp, old)

        # include the old ones as well as interpolated ones
        out = list(old) + list(interp)[1:-1]

        return out









class ConfigsParScan:
    def __init__(self, *params, sr_prop, inoc, load_saved, n_years):

        # params = rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2
        self.pars = params

        self.sr_prop = sr_prop

        self.inoc = inoc

        self.load_saved = load_saved

        self.n_years = n_years
  


    def _get_grid_conf(self, n_doses):

        conf = GridConfig(self.n_years, None, None, n_doses, 
                                    primary_inoculum=self.inoc)

        config_out = self._process_conf(conf)

        return config_out



    def _get_single_conf(self, dose1, dose2):

        conf = SingleConfig(self.n_years, None, None, 
                                dose1, dose1, dose2, dose2,
                                primary_inoculum=self.inoc)
        
        config_out = self._process_conf(conf)

        return config_out





    def _process_conf(self, Conf):

        Conf.sex_prop = self.sr_prop

        Conf.load_saved = self.load_saved

        Conf.add_string()

        config_out = self._update_par_scan_conf_str(Conf)
        
        return config_out



    
    def _update_par_scan_conf_str(self, Conf):
        
        rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2 = self.pars

        conf_str = Conf.config_string_img
        conf_str = conf_str.replace("grid", "param_scan")

        par_str = f"_fung_pars={om_1},{om_2},{delt_1},{delt_2},rf1={rfs1},rf2={rfs2},rfd={rfD}"
        par_str = par_str.replace(".", ",")
        conf_str = conf_str.replace(".png", par_str + ".png")
        
        Conf.config_string_img = conf_str
        
        saved_run_str = conf_str.replace(".png", ".pickle")
        
        Conf.config_string = saved_run_str.replace("figures/", "saved_runs/")

        return Conf




















# post process fns


def combine_PS_rand_outputs(config, seeds):

    df = pd.DataFrame()
    
    par_str = get_PS_rand_str(config)

    for seed in seeds:

        temporary = pd.read_csv(f"./param_scan_cluster/outputs/rand/par_scan/seed={seed}_{par_str}.csv")

        temporary["run"] = [seed*config["NIts"] + e for e in temporary["run"]]
        
        df = df.append(temporary, ignore_index=True)

    df.to_csv(f"./param_scan_cluster/outputs/rand/combined/output_{par_str}.csv")


# End of post process fns












class PostProcess:

    def __init__(self, par_str):
        self.df = pd.read_csv(f"param_scan_cluster/outputs/rand/combined/output_{par_str}.csv")
        self.par_str = par_str







    def get_maximum_along_contour_df(self):
        
        MDF = MaxAlongContourDF(self.df)

        df_out = MDF.df

        filename = f"./param_scan_cluster/outputs/rand/analysis/max_along_contour/df_{len(df_out)}.csv"

        print(f"Saving maximum along contour csv to: \n{filename}")

        df_out.to_csv(filename)

        self.max_along_contour_df = df_out


    











    def analyse_max_contour_df(self):

        df = copy.copy(self.max_along_contour_df)

        df['min_corner'] = df[["corner_01", "corner_10"]].min(axis=1)

        df = df[df['min_corner']>0]

        self.filtered_dataframe_outcome(df, "ERFB_Valid", "maxCont%")
        self.filtered_dataframe_outcome(df, "EqSelValid", "maxEqSel%")
        self.filtered_dataframe_outcome(df, None, "fullDose%")
        self.filtered_dataframe_outcome(df, None, "minEqDose%")




    @staticmethod
    def filtered_dataframe_outcome(df_in, strategy_valid, strategy):
    
        if strategy_valid is None:
            df = df_in
        else:
            df = df_in[df_in[strategy_valid]]

        mean = df[strategy].mean()

        conditional_mean = df[df[strategy]<100][strategy].mean()

        worked = df[df[strategy]>=100].shape[0]
        greater = df[df[strategy]>100].shape[0]
        lesser = df[df[strategy]<100].shape[0]
        equal = df[df[strategy]==100].shape[0]

        total = df.shape[0]

        sum_total = sum([greater, lesser, equal])

        out = dict(worked=worked, 
                    greater=greater,
                    lesser=lesser,
                    equal=equal,
                    total=total,
                    sum_total=sum_total,
                    work_pc=round(100*worked/sum_total,1),
                    mean=round(mean,1),
                    conditional_mean=round(conditional_mean,1),
                    strategy=strategy
                    )
        
        print(out)

        return out









    def analyse_failed(self):

        df = copy.copy(self.max_along_contour_df)

        fail = self._get_failed_runs(df)

        print("\n")
        print("These runs failed:\n")

        print(fail[['diff_from_opt',
                    'run',
                    'count',
                    # "corner_01",
                    # "corner_10",
                    "min_corner",
                    'maxGridEL',
                    'maxContEL',
                    'delta_RFB',
                    # 'max_cont_d1',
                    # 'max_cont_d2',
                    'max_cont_dose_either_fung',
                    'max_sing_res']].to_string())

        
    @staticmethod
    def _get_failed_runs(df):

        fail = df[df['maxCont%']<100]

        fail = fail.assign(diff_from_opt=lambda d: d['maxGridEL'] - d['maxContEL'])
        
        fail['max_sing_res'] = fail[["SR", "RS"]].max(axis=1)
        
        fail['min_corner'] = fail[["corner_01", "corner_10"]].min(axis=1)

        # fail = fail.sort_values(by=['diff_from_opt', 'max_cont_dose_either_fung'])

        # fail = fail[(fail['min_corner']>0) & (fail['max_cont_dose_either_fung']<1)]

        return fail

    



    def which_runs_worked_max_cont(self):
        df = copy.copy(self.max_along_contour_df)

        failed = df[df['maxCont%']<100]

        failed = failed[failed['min_corner']>0]

        runs_that_failed = failed["run"].unique()

        self.failed_pars = self.get_params_for_specific_runs(runs_that_failed)

        n_fail = self.failed_pars.shape[0]

        self.failed_pars.to_csv(f"./param_scan_cluster/outputs/failed_pars/failed_maxCont_{n_fail}.csv")

    





    def get_params_for_specific_runs(self, which_runs):

        par_df = self.df

        par_df = par_df.drop(["Unnamed: 0", "EL",
            "contour_level", "delta_RFB"], axis=1)

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = par_df[par_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out








    @staticmethod
    def _analyse_hi_lo_df(df_in):
        df = df_in[df_in['ERFB_Valid']]

        




    def check_high_or_low_dose(self):

        my_df = copy.copy(self.df)

        my_df['FD_BetterThanMin'] = my_df['fullDoseEL'] >= my_df['minEqDoseEL']

        strats = ["minEqDose", "fullDose"]
        
        for string in strats:
            my_df[string + "%"] = 100*my_df[string + "EL"]/my_df["maxGridEL"]
        
        my_df = my_df.drop(['Unnamed: 0', 'EL', 'econ'], axis=1)
                
        grouped = my_df.groupby(["run"]).first()

        df_out = pd.DataFrame(grouped)

        df_out = df_out.reset_index()
        
        df_out = df_out.sort_values(['FD_BetterThanMin', 'sr_prop'])
                
        filename = f"./param_scan_cluster/outputs/rand/analysis/high_or_low_dose/df_{len(df_out)}.csv"

        print(f"Saving high or low dose csv to: \n{filename}")
        
        df_out.to_csv(filename)

        self._analyse_hi_lo_df(df_out)




  
    
    
    def re_run_failures(self, NDoses, failed_run_indices=None):

        if failed_run_indices is None:
            failed_run_indices = list(range(self.failed_pars.shape[0]))


        for ii in failed_run_indices:

            pars = self.failed_pars.iloc[int(ii),:]
           
            print("\nRe-running run:", self.failed_pars.iloc[int(ii),:].run, "\n")

            grid_config, fung_params = self._get_grid_config_and_fung_pars(pars, NDoses)

            output = RunGrid(fung_params).grid_of_tactics(grid_config)

            conf_str = grid_config.config_string_img


            # plot output
            dose_grid_heatmap(output, grid_config, "FY", conf_str)
            
            # first_year_yield(output, grid_config)

            eq_RFB_contours(output, grid_config, title=f"Run={str(pars.run)}")






    def _get_grid_config_and_fung_pars(self, pars, NDoses):

        RP = RandomPars(None, None)

        RP._get_inoc_dict(pars["RS"], pars["SR"], pars["RR"])
        
        RP._get_fung_parms_dict(pars["omega_1"], pars["omega_2"], 
            pars["delta_1"], pars["delta_2"])
        
        fung_params = RP.fung_parms
    
        CPS = ConfigsParScan(pars["RS"], pars["SR"], pars["RR"],
            pars["omega_1"], pars["omega_2"], 
            pars["delta_1"], pars["delta_2"],
            sr_prop=pars["sr_prop"],
            inoc=RP.inoc,
            load_saved=True, 
            n_years=35)

        grid_config = CPS._get_grid_conf(NDoses)

        return grid_config, fung_params







class MaxAlongContourDF:

    def __init__(self, df_input):
        self.df = self._generate_max_along_contour_df(df_input)
        
    

    def _generate_max_along_contour_df(self, df_input):
        
        df_inter = self._get_intermediate_df(df_input)

        calc_cols_df = self._calculated_cols_df(df_inter)
        
        pars_count_df = self._params_and_count_df(df_inter)

        untidy = pd.concat([pars_count_df, calc_cols_df], axis=1)
        
        out = self._tidy_df(untidy)

        return out




    def _get_intermediate_df(self, data):

        data.fillna(0)

        data['maxAlongContour'] = data['maxContEL'] >= data['maxGridEL']

        strats = ["maxCont", "minEqDose", "fullDose", "maxEqSel"]
        
        for string in strats:
            data[string + "%"] = 100*data[string + "EL"]/data["maxGridEL"]
        
        out = data.drop(['Unnamed: 0', 'econ'], axis=1)

        return out



    def _calculated_cols_df(self, df_in):
        DS_df = pd.DataFrame()

        DS_df['min_dose_sums'] = df_in.groupby(["run"]).apply(lambda data: min(data['dose1'] + data['dose2']))
        DS_df['max_dose_sums'] = df_in.groupby(["run"]).apply(lambda data: max(data['dose1'] + data['dose2']))
        
        DS_df['min_opt_DS'] = df_in.groupby(["run"]).apply(self._min_opt_DS)
        DS_df['max_opt_DS'] = df_in.groupby(["run"]).apply(self._max_opt_DS)
        
        DS_df['max_cont_d1'] = df_in.groupby(["run"]).apply(lambda data: max(data['dose1']))
        DS_df['max_cont_d2'] = df_in.groupby(["run"]).apply(lambda data: max(data['dose2']))
        DS_df['max_cont_dose_either_fung'] = DS_df[['max_cont_d1', 'max_cont_d2']].max(axis=1)

        DS_df['min_DS_best'] = np.where(DS_df['min_dose_sums']==DS_df['min_opt_DS'], True, DS_df['min_dose_sums']-DS_df['min_opt_DS'])
        DS_df['max_DS_best'] = np.where(DS_df['max_dose_sums']==DS_df['max_opt_DS'], True, DS_df['max_dose_sums']-DS_df['max_opt_DS'])
        DS_df['max_or_min_DS_best'] = np.where((DS_df['min_dose_sums']==DS_df['min_opt_DS'])
                                        | (DS_df['max_dose_sums']==DS_df['max_opt_DS']),
                                            True, False)
        return DS_df

    
   
    def _min_opt_DS(self, data):
        
        df = self._get_opt_df(data)

        if df is None or not df.shape[0]:
            return "NA"

        return min(df['dose1']+df['dose2'])




    def _max_opt_DS(self, data):
        
        df = self._get_opt_df(data)

        if df is None or not df.shape[0]:
            return "NA"
        
        return max(df['dose1']+df['dose2'])




    @staticmethod
    def _get_opt_df(data):
        if not data['maxContEL'].shape[0]:
            return None
        
        maxEL_this_run = float(list(data['maxContEL'])[0])

        df = data[data['EL']==maxEL_this_run]
        
        return df





    def _params_and_count_df(self, df_in):
        counting = df_in.groupby(["run"]).size()

        grouped = df_in.groupby(["run"]).first()

        grouped['count'] = counting
        
        out = grouped.reset_index()

        return out





    def _tidy_df(self, df):

        df['min_corner'] = df[["corner_01", "corner_10"]].min(axis=1)
        
        df = df[df["fullDoseEL"]>0]

        out = df.sort_values(['ERFB_Valid', 'min_corner', 'EqSelValid', 'maxCont%', 'maxEqSel%'])

        return out
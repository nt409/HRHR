import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
# import matplotlib.pyplot as plt
import itertools


from runHRHR.config_classes import GridConfig, SingleConfig
from utils.params import PARAMS
from utils.functions import RunGrid, RunSingleTactic, \
    logit10_difference, EqualResFreqBreakdownArray, \
    EqualSelectionArray
from utils.plotting import dose_grid_heatmap, eq_RFB_contours, \
    DoseSpaceScenariosPlot

from .cont_opt import ContourStratDFGetter

# TOC:

# Utility fns
# - get_PS_rand_str



# CLASSES:
# ParamScanRand

# SinglePSRun

# in cont_opt:
# ContourStratDFGetter
# ContourList
# SummaryDF

# FYSelFinder
# RFBFinder

# OtherStratDFGetter

# CheckStrategyUsingIVT
# ContourPassesThroughChecker

# RandomPars


# post process fns:
# - combine_PS_rand_outputs

# CLASSES:
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
        
        self.key_list = ["summary",
                            "RFB_df",
                            "EqSel_df"]
    





    def run(self, seed):
        """
        Run random scan over uniform dists
        """

        dict_of_dfs = self._get_scan_output(seed)

        self.save_output(dict_of_dfs, seed)
        

    
    
    
    def save_output(self, df_dict, seed):
        
        par_str = get_PS_rand_str(self.config)

        for key in self.key_list:

            filename = f"./param_scan_cluster/outputs/rand/par_scan/{key}_seed={seed}_{par_str}.csv"
            
            print("\n")
            print(f"Random Scan {key}, saved as:\n {filename}")

            df_dict[key].to_csv(filename, index=False)







    def _get_scan_output(self, seed):

        np.random.seed(seed)
        
        out = self.initialise_output_dict()

        N_ITS = self.config["NIts"]

        for run_index in tqdm(range(N_ITS)):

            new_run_ind = seed*N_ITS + run_index

            dict_of_dfs = SinglePSRun(self.config, new_run_ind).output

            out = self.update_dfs(out, dict_of_dfs)

        return out
    




    def initialise_output_dict(self):
        
        out = {}
        
        for key in self.key_list:
            out[key] = pd.DataFrame()
        
        return out





    @staticmethod
    def update_dfs(out, dict_of_dfs):

        for key in dict_of_dfs:
            new_df = dict_of_dfs[key]
            out[key] = pd.concat([out[key], new_df], axis=0)

        return out

















class SinglePSRun:
    def __init__(self, config, run_index) -> None:
        
        self.config = config
        
        self.run_index = run_index

        rand_pars = self._get_rand_pars_obj()
        
        grid_output = self._get_grid_output_this_run(rand_pars)
        
        without_RI = self._get_single_PS_run_dfs(rand_pars, grid_output)

        self.output = self._add_run_index_to_dfs(without_RI)






    def _get_single_PS_run_dfs(self, rand_pars, grid_output):

        par_df = self._get_par_df(rand_pars)
        
        RFB_dfs = ContourStratDFGetter(rand_pars, EqualResFreqBreakdownArray,
                                RFBFinder, grid_output, 0, self.config['n_cont_points'], "RFB")

        EqSel_dfs = ContourStratDFGetter(rand_pars, EqualSelectionArray, 
                                FYSelFinder, grid_output, 0.5, self.config['n_cont_points'], "EqSel")

        other_strats_df = OtherStratDFGetter(grid_output).df

        IVT_RFB_df = CheckStrategyUsingIVT(grid_output, 
                    EqualResFreqBreakdownArray, 0, "RFB").df

        IVT_EqSel_df = CheckStrategyUsingIVT(grid_output,
                    EqualSelectionArray, 0.5, "EqSel").df

        summary = pd.concat([par_df, RFB_dfs.summary, 
                    EqSel_dfs.summary, other_strats_df,
                    IVT_RFB_df, IVT_EqSel_df],
                    axis=1)

        return dict(
            summary = summary,
            RFB_df = RFB_dfs.df,
            EqSel_df = EqSel_dfs.df,
            )
        




    def _get_rand_pars_obj(self):
    
        RP = RandomPars(self.config, self.run_index)

        RP.find_pars()

        return RP




    def _get_grid_output_this_run(self, rand_pars):
        
        gridconfig = rand_pars._get_grid_conf(self.config["grid_number"])

        out= RunGrid(rand_pars.fung_parms).grid_of_tactics(gridconfig)

        return out 




    @staticmethod
    def _get_par_df(rand_pars):
        data = rand_pars._get_this_run_all_parms_dict()
        
        return pd.DataFrame([data])




    def _add_run_index_to_dfs(self, dict_dfs):

        for key in dict_dfs:
            dict_dfs[key]['run'] = [self.run_index]*dict_dfs[key].shape[0]
        
        return dict_dfs





























class FYSelFinder:
    def __init__(self, single_run) -> None:
        self.cont_quant = self._get_FY_selection(single_run)

    
    
    def _get_FY_selection(self, single_run):
        
        start_freqs = single_run['start_of_season']
        
        sr1 = start_freqs['RS'][1]/start_freqs['RS'][0]
        sr2 = start_freqs['SR'][1]/start_freqs['SR'][0]
        
        try:
            return sr1/(sr1+sr2)
        except:
            return None







class RFBFinder:

    def __init__(self, single_run) -> None:
        self.cont_quant = self._get_delta_RFB(single_run)


    def _get_delta_RFB(self, sing_run):
        
        fy_in = sing_run['failure_year']

        fy = int(fy_in)

        res_vecs = sing_run['res_vec_dict']
        
        rf1 = res_vecs['f1'][fy]
        
        rf2 = res_vecs['f2'][fy]

        try:
            return logit10_difference(rf1, rf2)
        except:
            print(f"problem with calculating delta_RFB for: {rf1, rf2}")
            return "NA"











class OtherStratDFGetter:
    def __init__(self, grid_output) -> None:
        FYs = grid_output['FY']
        
        self.df = self._get_strategy_outcomes(FYs)
    


    def _get_strategy_outcomes(self, FYs):
        
        minEqDoseELVec = [float(FYs[i, i]) for i in range(FYs.shape[0])]
        
        data = dict(
                max_grid_EL = np.amax(FYs),
                minEqDoseEL = self._get_first_non_zero_element(minEqDoseELVec),
                corner_00 = FYs[0, 0],
                corner_10 = FYs[-1, 0],
                corner_01 = FYs[0, -1],
                fullDoseEL = FYs[-1, -1])

        out = pd.DataFrame([data])

        return out
    
    
    
    @staticmethod
    def _get_first_non_zero_element(vec):
        filtered = list(filter(lambda x: x>0, vec))

        if not filtered:
            return "NA"

        return filtered[0]

































class CheckStrategyUsingIVT:
    """
    Use "intermediate value theorem" method:
    
    - find optimal region and see if contains points above and below contour
    - if so, then contour passes through (provided connected)

    """
    def __init__(self, grid_output, strategy_class, level, name) -> None:

        self.FYs = grid_output['FY']
        
        self.strategy_obj = strategy_class(grid_output)

        self.level = level

        self.name = name
        
        self.check_if_gives_optimum()

        self.find_best_value_this_strat()

        self.df = self.get_df(name)


    def check_if_gives_optimum(self):
        FYs = self.FYs
        
        opt_region = FYs!=np.amax(FYs)
        
        if not self.strategy_obj.is_valid:
            self.strat_works = f"Strategy {self.name} is invalid"
            return None
        
        strat_array = self.strategy_obj.array

        opt_strat = np.ma.masked_array(strat_array, mask=opt_region)
        
        self.strat_works = ContourPassesThroughChecker(opt_strat, self.level).passes_through



    def find_best_value_this_strat(self):
        FYs = self.FYs

        if not self.strategy_obj.is_valid:
            self.best_value = f"Strategy {self.name} is invalid"
            return None

        strat_array = self.strategy_obj.array

        for k in range(int(np.amax(FYs))):
            EL = np.amax(FYs)-k
            
            opt_region = FYs<EL

            opt_strat = np.ma.masked_array(strat_array, mask=opt_region)

            worked = ContourPassesThroughChecker(opt_strat, self.level).passes_through
        
            if worked:
                self.best_value = EL
                return None

        self.best_value = EL


    def get_df(self, name):
        data = {f"best_region_{name}": self.strat_works,
                f"best_value_{name}": self.best_value}

        return pd.DataFrame([data])












class ContourPassesThroughChecker:
    def __init__(self, array, level) -> None:
        
        self.passes_through = self.check_if_passes_through_region(array, level)




    def check_if_passes_through_region(self, array, level) -> bool:
        """
        Check if contour passes through optimal region.
        
        If False, might be because:
        - optimal region is not connected, or
        - all contour values are positive/negative in the region
        - grid sufficiently dense

        If True, contour passes through optimal region :)

        """
        
        includes_level = (np.amin(array)<level and np.amax(array)>level)

        if not includes_level:
            return False
        
        return self.check_is_connected(array, level)



    
    def check_is_connected(self, array, level) -> bool:

        for i, j in itertools.product(*(range(array.shape[ii]) for ii in [0,1])):

            if not array[i,j]:
                continue

            elif array[i,j]==level:
                return True

            else:
                valid = self.check_neighbours(array, i, j, level)
            
            if valid:
                return True

        return False
    


    def check_neighbours(self, array, i, j, level) -> bool:
        """
        For array only including optimal region:

        Check if a cell has any neighbours which have a differing sign for contour.

        If True, there is some intermediate value for doses between the two cells
        which is optimal and along the contour.

        If False, t
        """
        
        self_is_positive = array[i,j]>level

        cell_iterator = self.neighbours((i,j), array.shape[0])

        for ii, jj in cell_iterator:
                
            if not array[ii,jj]:
                continue

            neighbour_is_positive = array[ii,jj]>level

            if neighbour_is_positive!=self_is_positive:
                return True

        return False

    
    @staticmethod
    def neighbours(cell, shape):
        """
        Find neighbouring cells, excluding self and excluding any cells with:
        -: index < 0, or
        -: index > shape.

        Returns a generator

        """
        for new_cell in itertools.product(*(range(n-1, n+2) for n in cell)):
            if new_cell != cell and all(0 <= n < shape for n in new_cell):
                yield new_cell







        
        
























class RandomPars:
    """
    Generates random parameters for model run
    
    Can return
    - pars
    - grid config
    - single config

    """

    def __init__(self, config, run_index) -> None:

        self._config = config

        self._run_index = run_index



    def find_pars(self):
        """
        returns rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop
        """

        conf = self._config

        pars_are_valid = False

        while pars_are_valid is False:
        
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

            fungicide_params = self.get_fung_parms_dict(om_1, om_2, delt_1, delt_2)

            pathogen_pars = dict(rfs1=rfs1,
                rfs2=rfs2,
                rfD=rfD, 
                sr_prop=sr_prop)

            pars_are_valid = self._check_validity(fungicide_params, pathogen_pars)

        
        self.get_inoc_dict(rfs1, rfs2, rfD)
        
        self.fung_parms = self.get_fung_parms_dict(om_1, om_2, delt_1, delt_2)

        self.path_and_fung_pars = rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2

        self.sr_prop = sr_prop








    def _check_validity(self, fungicide_params, pathogen_pars):

        self.get_inoc_dict(pathogen_pars['rfs1'],
                    pathogen_pars['rfs2'],
                    pathogen_pars['rfD'])
        
        for ii, jj in [[1,0], [0,1], [1,1]]:
            this_dose_yield = self._check_yield(fungicide_params, pathogen_pars, ii, jj)

            if not this_dose_yield>95:
                return False
        
        return True





    def _check_yield(self, fungicide_params, pathogen_pars, d1, d2):

        ConfigSingleRun = SingleConfig(1, None, None, d1, d1, d2, d2, primary_inoculum=self.inoc)
        
        ConfigSingleRun.sex_prop = pathogen_pars['sr_prop']

        ConfigSingleRun.load_saved = False

        this_run = RunSingleTactic(fungicide_params).run_single_tactic(ConfigSingleRun)
        
        yield_out = this_run['yield_vec'][0]

        return yield_out





    def get_inoc_dict(self, rfs1, rfs2, rfD):
        self.inoc = dict(RR = rfD,
                         RS = rfs1,
                         SR = rfs2,
                         SS = 1 - rfD - rfs1 - rfs2)

    

    def get_fung_parms_dict(self, omega_1, omega_2, delta_1, delta_2):
        return dict(
                            omega_1 = omega_1,
                            omega_2 = omega_2,
                            theta_1 = PARAMS.theta_1,
                            theta_2 = PARAMS.theta_2,
                            delta_1 = delta_1,
                            delta_2 = delta_2)


    
    def _get_this_run_all_parms_dict(self):
        out = {**self.fung_parms,
                **self.inoc,
                "sr_prop": self.sr_prop,
                "run": self._run_index}
        return out





    # single and grid config finders:

    def _get_grid_conf(self, n_doses):
        
        conf = GridConfig(self._config['n_years'], None, None, n_doses,
                                    primary_inoculum=self.inoc)

        config_out = self._process_conf(conf)

        return config_out





    def _get_single_conf(self, dose1, dose2):

        conf = SingleConfig(self._config['n_years'], None, None,
                                dose1, dose1, dose2, dose2,
                                primary_inoculum=self.inoc)
        
        config_out = self._process_conf(conf)

        return config_out




    def _process_conf(self, Conf):

        Conf.sex_prop = self.sr_prop

        Conf.load_saved = self._config['load_saved']

        Conf.add_string()

        config_out = self._update_par_scan_conf_str(Conf)
        
        return config_out



    
    def _update_par_scan_conf_str(self, Conf):
        
        rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2 = self.path_and_fung_pars

        conf_str = Conf.config_string_img
        conf_str = conf_str.replace("grid", "param_scan")

        par_str = (f"_fung_pars={round(om_1,6)},{round(om_2,6)}," + 
                    f"{round(delt_1,6)},{round(delt_2,6)}," +
                    f"rf1={round(rfs1,10)},rf2={round(rfs2,10)},rfd={round(rfD,15)}")
        
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

        temporary = pd.read_csv(f"./param_scan_cluster/outputs/rand/par_scan/summary_seed={seed}_{par_str}.csv")

        df = df.append(temporary, ignore_index=True)

    df.to_csv(f"./param_scan_cluster/outputs/rand/combined/output_summary_{par_str}.csv")


# End of post process fns












class PostProcess:

    def __init__(self, par_str):
        df_in = pd.read_csv(f"param_scan_cluster/outputs/rand/combined/output_summary_{par_str}.csv")
        self.df = df_in.drop(["Unnamed: 0"], axis=1)

        self.par_str = par_str
        







    def get_maximum_along_contour_df(self):
        self.max_along_contour_df = MaxAlongContourDF(self.df).df


    











    def analyse_max_contour_df(self):

        df = copy.copy(self.max_along_contour_df)

        df = df[df['min_corner']>0]

        eq_sel_df = self._get_non_null_df(df, "EqSel_worked")

        self.filtered_dataframe_outcome(df, "RFB_maxCont%")
        self.filtered_dataframe_outcome(df, "fullDose%")
        self.filtered_dataframe_outcome(df, "minEqDose%")
        self.filtered_dataframe_outcome(eq_sel_df, "EqSel_maxCont%")
        self.filtered_dataframe_outcome(eq_sel_df, "EqSel_lowDoseMax%")




    def _get_non_null_df(self, df_in, col_to_check_if_null):
        invalid_runs = df_in[col_to_check_if_null].isnull()
        return df_in[~invalid_runs]



    def filtered_dataframe_outcome(self, df, strategy):

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
        
        if strategy=="RFB_maxCont%":
            self.check_IVT_method(df, "RFB")
        elif strategy=="EqSel_maxCont%":
            self.check_IVT_method(df, "EqSel")

        print(out)

        return out



    def check_IVT_method(self, df, method):

        if method=="EqSel":
            IVT_true = df[f'best_region_{method}']=="True"
            IVT_false = df[f'best_region_{method}']=="False"
            IVT_NA = df[f'best_region_{method}'].isin(["True", "False"])
        else:
            IVT_true = df[f'best_region_{method}']==True
            IVT_false = df[f'best_region_{method}']==False
            IVT_NA = df[f'best_region_{method}'].isin([True, False])

        opt_true = (df[f'{method}_worked']==True)
        opt_false = (df[f'{method}_worked']==False)
        opt_NA = (df[f'{method}_worked'].isnull())

        out = dict(
            IVT_method_T = df[IVT_true].shape[0],
            IVT_method_F = df[IVT_false].shape[0],
            IVT_method_NA = df[~IVT_NA].shape[0],
            
            opt_method_T = df[opt_true].shape[0],
            opt_method_F = df[opt_false].shape[0],
            opt_method_NA = df[opt_NA].shape[0],


            both_succeeded = df[(opt_true & IVT_true)].shape[0],
            both_failed = df[(~opt_true & ~IVT_true)].shape[0],
            
            either_succeeded = df[(opt_true | IVT_true)].shape[0],
            opt_but_not_IVT = df[(opt_true & ~IVT_true)].shape[0],
            IVT_but_not_opt = df[(~opt_true & IVT_true)].shape[0],
            )

        failed = df[(~opt_true & ~IVT_true)]
        

        print("\n")

        print(out)

        print("\n")
        
        print("These runs failed on both methods:")
        
        print(failed[[f'max_grid_EL', f'{method}_maxContEL', f'best_value_{method}']])






    def analyse_failed(self):

        df = copy.copy(self.max_along_contour_df)

        fail = self._get_failed_runs(df)

        fail = fail.assign(RFB_diff_from_opt=lambda d: d['max_grid_EL'] - d['RFB_maxContEL'])


        print("\n")
        print("These runs failed:\n")

        print(fail[['RFB_diff_from_opt',
                    'run',
                    "min_corner",
                    'max_grid_EL',
                    'RFB_maxContEL',
                    ]].to_string())

        
    @staticmethod
    def _get_failed_runs(df):
        return df[df['RFB_maxCont%']<100]

    



    def which_runs_worked_max_cont(self):
        df = copy.copy(self.max_along_contour_df)

        failed = df[df['RFB_maxCont%']<100]

        runs_that_failed = failed["run"].unique()

        failed_pars = self.get_params_for_specific_runs(runs_that_failed)

        n_fail = failed_pars.shape[0]

        failed_pars.to_csv(f"./param_scan_cluster/outputs/failed_pars/failed_{n_fail}.csv")

    





    def get_params_for_specific_runs(self, which_runs):

        par_df = copy.copy(self.df)

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = par_df[par_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out








        




    def check_high_or_low_dose(self):

        my_df = copy.copy(self.df)

        my_df['FD_BetterThanMin'] = my_df['fullDoseEL'] >= my_df['minEqDoseEL']

        strats = ["minEqDose", "fullDose"]
        
        for string in strats:
            my_df[string + "%"] = 100*my_df[string + "EL"]/my_df["max_grid_EL"]
        
        grouped = my_df.groupby(["run"]).first()

        df_out = pd.DataFrame(grouped)

        df_out = df_out.reset_index()
        
        df_out = df_out.sort_values(['FD_BetterThanMin', 'sr_prop'])
                
        filename = f"./param_scan_cluster/outputs/rand/analysis/high_or_low_dose/df_{len(df_out)}.csv"

        print(f"Saving high or low dose csv to: \n{filename}")
        
        df_out.to_csv(filename)




  
    
    
    def re_run(self, NDoses, run_indices):

        df_test = self.get_params_for_specific_runs(run_indices)


        for ii in range(df_test.shape[0]):

            pars = df_test.iloc[int(ii),:]
           
            print("\nRe-running run:", df_test.iloc[int(ii),:].run, "\n")

            grid_config, fung_params = self._get_grid_config_and_fung_pars(pars, NDoses)

            # grid_config.load_saved = False

            grid_output = RunGrid(fung_params).grid_of_tactics(grid_config)

            conf_str = grid_config.config_string_img

            FY = grid_output['FY']
            opt_region = FY == np.amax(FY)
                
            n_opt_doses = opt_region.sum()

            print(f"Number of optimal dose combos: {n_opt_doses}")

            # plot output
            # dose_grid_heatmap(grid_output, grid_config, "FY", conf_str)
            
            # eq_RFB_contours(grid_output, grid_config, title=f"Run={str(pars.run)}")
            DoseSpaceScenariosPlot(grid_output, conf_str)






    @staticmethod
    def _get_grid_config_and_fung_pars(pars, NDoses):

        config = {'load_saved': True, 'n_years': 35}

        RP = RandomPars(config, None)

        RP.get_inoc_dict(pars["RS"], pars["SR"], pars["RR"])
        
        fung_params = RP.get_fung_parms_dict(pars["omega_1"], pars["omega_2"], 
                                pars["delta_1"], pars["delta_2"])
        
        RP.sr_prop = pars["sr_prop"]

        RP.path_and_fung_pars = (pars["RS"], pars["SR"], pars["RR"], pars["omega_1"],
                    pars["omega_2"], pars["delta_1"], pars["delta_2"])

        grid_config = RP._get_grid_conf(NDoses)

        return grid_config, fung_params














class MaxAlongContourDF:

    def __init__(self, df_input):
        self.get_and_save(df_input)
    


    def get_and_save(self, df_input):
        df_inter = self._get_intermediate_df(df_input)

        self.df = self._tidy_df(df_inter)
        
        self._save_df()





    def _get_intermediate_df(self, data):

        data.fillna(0)

        data = pd.DataFrame(data)
        
        data['maxAlongContour'] = data['RFB_maxContEL'] >= data['max_grid_EL']

        strats = ["RFB_maxCont", "EqSel_maxCont", "minEqDose", "fullDose", "EqSel_lowDoseMax"]
        
        for string in strats:
            data.loc[:, string + "%"] = 100*data[string + "EL"]/data["max_grid_EL"]
        
        return data

    



    def _tidy_df(self, df):

        df['min_corner'] = df[["corner_01", "corner_10"]].min(axis=1)

        # avoid the "-1" case where has never failed:
        if df[df["fullDoseEL"]<=0].shape[0]:
            n_never_fail = df[df["fullDoseEL"]<=0].shape[0]
            print(f"{n_never_fail} runs never failed - need longer n_years. Filtering them out for now.")
            df = df[df["fullDoseEL"]>0]
        
        out = df.sort_values(['min_corner', 'RFB_maxCont%', 'EqSel_maxCont%'])

        return out




    def _save_df(self):
        df_out = self.df
        
        filename = f"./param_scan_cluster/outputs/rand/analysis/max_along_contour/df_{len(df_out)}.csv"

        print(f"Saving maximum along contour csv to: \n{filename}")

        df_out.to_csv(filename)
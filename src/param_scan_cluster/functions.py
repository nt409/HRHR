import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
import matplotlib.pyplot as plt
import itertools


from runHRHR.config_classes import GridConfig, SingleConfig
from utils.params import PARAMS
from utils.functions import RunGrid, RunSingleTactic, \
    logit10_difference, EqualResFreqBreakdownArray, \
    EqualSelectionArray
from utils.plotting import dose_grid_heatmap, eq_RFB_contours

# TOC:

# Utility fns
# ParamScanRand
# SinglePSRun

# ContourStratDFGetter
# FYSelFinder
# RFBFinder

# OtherStratDFGetter

# ContourPassesThroughChecker

# CheckStrategyUsingIVT

# RandomPars
# ContourList

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
        
        RFB_dfs = ContourStratDFGetter(rand_pars, grid_output,
                            EqualResFreqBreakdownArray, 0,
                            RFBFinder, "RFB")

        

        EqSel_dfs = ContourStratDFGetter(rand_pars, grid_output,
                            EqualSelectionArray, 0.5,
                            FYSelFinder, "EqSel")

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




















class ContourStratDFGetter:
    def __init__(self, rand_pars, grid_output, Strategy,
                    contour_level, ContQuantFinder, strat_name) -> None:

        self.rand_pars = rand_pars

        self.strat_name = strat_name

        self.strat_obj = Strategy(grid_output)

        self.max_grid_EL = np.amax(grid_output['FY'])

        self.ContQuantFinder = ContQuantFinder

        self.level = contour_level

        self.df = self._get_data()
        
        self.summary = SummaryDF(self.df, self.strat_name, self.level, 
                                self.max_grid_EL).df





    def _get_data(self):
        cntr = self._find_contours()
        
        return self._get_dict_this_cntr(cntr)




    
    def _find_contours(self):
        
        if not self.strat_obj.is_valid:
            cntr_out = {}
            return cntr_out
        else:    
            cntr_out = ContourList(self.strat_obj.array, level=self.level).cont_dict
            return cntr_out






    def _get_dict_this_cntr(self, contour):
        
        df_this_run = pd.DataFrame()

        if not contour:
            return df_this_run


        for dose1, dose2 in zip(contour['x'], contour['y']):            

            this_contour_dict = self._get_data_these_doses_on_contour(dose1, dose2)
            
            data = {"contour_level": contour['level'],
                     "dose1": dose1,
                     "dose2": dose2,
                     "DS": dose1 + dose2,
                     **this_contour_dict}

            df_this_run = df_this_run.append(data, ignore_index=True)
        
        
        
        return df_this_run





    def _get_data_these_doses_on_contour(self, dose1, dose2):
        """
        Returns dict with
        - contQuant
        - EL
        """

        this_dose_conf = self.rand_pars._get_single_conf(dose1, dose2)
                
        self.single_run = RunSingleTactic(self.rand_pars.fung_parms).run_single_tactic(this_dose_conf)
        
        out = self._get_these_doses_dict()
        
        return out








    def _get_these_doses_dict(self):
        """
        Returns dict with
        - contQuant
        - EL
        # - econ
        """
        sing_run = self.single_run

        strat_name = self.strat_name

        # econ="NA",

        out = {
                f"{strat_name}_contQuant": self.ContQuantFinder(sing_run).cont_quant,
                f"{strat_name}_EL": sing_run['failure_year'],
                }
        
        return out







class SummaryDF:

    def __init__(self, df, strat_name, level, max_grid_EL) -> None:
        self.df = df        
        self.strat_name = strat_name
        self.level = level
        self.max_grid_EL = max_grid_EL

        self.df = self._get_summary_df(df, strat_name)





    def _get_summary_df(self, df, strat_name):
        df = self.df

        lowDS_val = self._get_low_dose_max_EL()
        medDS_val = self._get_med_dose_max_EL()
        highDS_val = self._get_high_dose_max_EL()

        
        min_opt_dist_from_cntr = self._get_min_opt_dist_from_contour()
        

        data= {
            f"{strat_name}_minContEL": min(df[f"{strat_name}_EL"]),
            f"{strat_name}_maxContEL": max(df[f"{strat_name}_EL"]),

            f"{strat_name}_minContQuant": min(df[f"{strat_name}_contQuant"]),
            f"{strat_name}_maxContQuant": max(df[f"{strat_name}_contQuant"]),
            f"{strat_name}_min_opt_dist_from_contour": min_opt_dist_from_cntr,

            f"{strat_name}_minDS": min(df['DS']),
            f"{strat_name}_maxDS": max(df['DS']),
            
            f"{strat_name}_worked": self.max_grid_EL==max(df[f"{strat_name}_EL"]),
            
            f"{strat_name}_valid": ((min(df[f"{strat_name}_contQuant"])<self.level) and
                                        (max(df[f"{strat_name}_contQuant"])>self.level)),

            f"{strat_name}_lowDoseMaxEL": lowDS_val,
            f"{strat_name}_medDoseMaxEL": medDS_val,
            f"{strat_name}_highDoseMaxEL": highDS_val,
            }

        return pd.DataFrame([data])


    


    
    
    def _get_low_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + 0.2*(max(df['DS']) - min(df['DS']))
        
        filt = df[df['DS'] < DS_thres_low]

        return max(filt[f"{self.strat_name}_EL"])
    





    def _get_med_dose_max_EL(self):
        df = self.df

        DS_thres_low = min(df['DS']) + 0.2*(max(df['DS']) - min(df['DS']))
        DS_thres_high = min(df['DS']) + 0.8*(max(df['DS']) - min(df['DS']))

        filt = df[((df['DS'] >= DS_thres_low) & (df['DS'] <= DS_thres_high))]

        return max(filt[f"{self.strat_name}_EL"])

    




    def _get_high_dose_max_EL(self):
        df = self.df

        DS_thres_high = min(df['DS']) + 0.8*(max(df['DS']) - min(df['DS']))

        filt = df[df['DS'] > DS_thres_high]

        return max(filt[f"{self.strat_name}_EL"])





    def _get_min_opt_dist_from_contour(self):
        df = self.df

        opt_df = df[df[f"{self.strat_name}_EL"]==self.max_grid_EL]

        vec = abs(opt_df[f"{self.strat_name}_contQuant"] - self.level)
        if len(vec):
            out = min(vec)
            return out
        else:
            return "NA"
        










class FYSelFinder:
    def __init__(self, single_run) -> None:
        self.single_run = single_run
        self.cont_quant = self._get_FY_selection()

    
    
    def _get_FY_selection(self):
        
        start_freqs = self.single_run['start_of_season']
        
        sr1 = start_freqs['RS'][1]/start_freqs['RS'][0]
        sr2 = start_freqs['SR'][1]/start_freqs['SR'][0]
        
        try:
            return sr1/(sr1+sr2)
        except:
            return None







class RFBFinder:

    def __init__(self, single_run) -> None:
        self.single_run = single_run
        self.cont_quant = self._get_delta_RFB()


    def _get_delta_RFB(self):
        sing_run = self.single_run
        
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

        par_str = f"_fung_pars={om_1},{om_2},{delt_1},{delt_2},rf1={rfs1},rf2={rfs2},rfd={rfD}"
        par_str = par_str.replace(".", ",")
        conf_str = conf_str.replace(".png", par_str + ".png")
        
        Conf.config_string_img = conf_str
        
        saved_run_str = conf_str.replace(".png", ".pickle")
        
        Conf.config_string = saved_run_str.replace("figures/", "saved_runs/")

        return Conf




    










class ContourList:
    """
    Takes an array and a contour level (i.e. value to aim for) and generates 
    a dictionary with x, y values, the contour level and the max dose
    """

    def __init__(self, z, level, min_n_pts_alng_cntr=12, plt_conts=False) -> None:
        
        self.cont_dict = self.find_contours(z, level, min_n_pts_alng_cntr)

        if plt_conts:
            self.plot_it()

    
    
    
    def find_contours(self, z, level, min_n_pts_alng_cntr):

        cs = self._get_contour_set(z, level)
        
        cntr_out = self._get_contours(cs, level)

        if cntr_out and len(cntr_out['x'])<min_n_pts_alng_cntr:
            cntr_out = self._interpolate_contours(cntr_out, min_n_pts_alng_cntr)
        return cntr_out

      


    @staticmethod
    def _get_contour_set(z, level):

        x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]

        cs = plt.contour(x, y, z, levels=[level])
        
        return cs





    def _get_contours(self, cs, level):
        """
        Takes a contour set and returns a list of dictionaries containing:
        - x values
        - y values
        - the contour level
        - max dose along the contour
        """

        if not cs.allsegs:
            print("Warning: conts was False!")
            # why not?
            # suspect too many 'nan's to create a proper contour
            return {}

        else:
            cont = cs.allsegs[0][0]
            cont_dict = self._get_contour_dict(level, cont)
            return cont_dict




    @staticmethod
    def _get_contour_dict(level, cont):
        
        x_vals = cont[:,0]
        y_vals = cont[:,1]

        # x_max = max(x_vals)
        # y_max = max(y_vals)

        # print(f"max dose along contour: (x,y) = {round(x_max,2), round(y_max,2)}")

        out = dict(x = x_vals,
                   y = y_vals,
                #    max_dose = max(x_max, y_max),
                   level = level)

        return out

    


    def _interpolate_contours(self, cntr, num=15):
        """
        Takes old contours and adds in some interpolated values
        
        Result is more densely populated list of values (approximately) 
        along the contour
        """

        xx = copy.copy(cntr['x'])
        yy = copy.copy(cntr['y'])

        cntr['x'] = self._interp_vector(xx, num)
        cntr['y'] = self._interp_vector(yy, num)
        
        return cntr
    



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




    def plot_it(self):
        cont = self.cont_dict
        plt.scatter(cont['x'], cont['y'])
        plt.show()




















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

        if strategy=='maxEqSel%':
            print(df[strategy].isnull().sum())
            exit()

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

            grid_output = RunGrid(fung_params).grid_of_tactics(grid_config)

            conf_str = grid_config.config_string_img


            # plot output
            dose_grid_heatmap(grid_output, grid_config, "FY", conf_str)
            
            eq_RFB_contours(grid_output, grid_config, title=f"Run={str(pars.run)}")






    def _get_grid_config_and_fung_pars(self, pars, NDoses):

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

        data = pd.DataFrame(data)
        
        data['maxAlongContour'] = data['maxContEL'] >= data['maxGridEL']

        strats = ["maxCont", "minEqDose", "fullDose", "maxEqSel"]
        
        for string in strats:
            data.loc[:, string + "%"] = 100*data[string + "EL"]/data["maxGridEL"]
        
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

        # avoid the "-1" case where has never failed:
        if df[df["fullDoseEL"]<=0].shape[0]:
            n_never_fail = df[df["fullDoseEL"]<=0].shape[0]
            print(f"{n_never_fail} runs never failed - need longer n_years. Filtering them out for now.")
            df = df[df["fullDoseEL"]>0]
        
        out = df.sort_values(['ERFB_Valid', 'min_corner', 'EqSelValid', 'maxCont%', 'maxEqSel%'])

        return out
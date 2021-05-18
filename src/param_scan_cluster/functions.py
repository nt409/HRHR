import itertools
import pandas as pd
from tqdm import tqdm
import copy
import numpy as np
from math import log10, ceil, log, exp
import matplotlib.pyplot as plt
# from sklearn import linear_model


from runHRHR.config_classes import GridConfig, SingleConfig
from utils.params import PARAMS
from utils.functions import non_decreasing, non_increasing, RunGrid, \
    RunSingleTactic, logit10_difference
from utils.plotting import dose_grid_heatmap, eq_RFB_contours

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






def _get_contours(z, levels=None, step=1):
    """
    Takes an array z and returns a list of dictionaries containing:
    - x values
    - y values
    - the contour level
    """
    
    x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]

    cs = plt.contour(x, y, z, levels=levels)

    output = []

    for level, conts in zip(levels, cs.allsegs):
        
        if not conts:
            # why not?
            # suspect too many 'nan's to create a proper contour
            continue

        cont = conts[0]
        
        vals = dict(x=cont[:,0], y=cont[:,1])
            
        output.append(
            dict(
                x = [vals['x'][0]] + list(vals['x'][1:-2:step]) + [vals['x'][-1]],
                y = [vals['y'][0]] + list(vals['y'][1:-2:step]) + [vals['y'][-1]],
                # level = levels[ind]
                level = level
                )
            )
    
    if False:
        for ind in range(len(output)):
            plt.scatter(output[ind]['x'], output[ind]['y'])

        plt.show()
        
    return output









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


    def _get_grid_conf(self, rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop, load_saved,
                            n_years=None, n_doses=None):

        if n_years is None:
            n_years = self.n_years

        if n_doses is None:
            n_doses = self.NDoses


        Conf = GridConfig(n_years, None, None, n_doses, primary_inoculum=self.inoc)
        
        Conf.sex_prop = sr_prop

        Conf.load_saved = load_saved

        Conf.add_string()

        conf_str = _get_param_scan_conf_str(rfs1, rfs2, rfD, om_1, om_2,
                                                        delt_1, delt_2, Conf)
        
        Conf.config_string_img = conf_str

        Conf.config_string = _get_updated_conf_str(conf_str)

        return Conf







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






class ParamScanRand(ParamScan):
    def __init__(self, config) -> None:
        
        super().__init__()
        
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
            print("\n")
            print("valid params")
            print("\n")
            return True
        else:
            print("\n")
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
        self.this_run_params = {**self.fung_parms,
                **self.inoc,
                "sr_prop": sr_prop,
                "run": f"{run_index}"}





    def _get_data_this_conf_and_contrs(self, *params):
        
        df_this_run = pd.DataFrame()

        if not self.contours:
            self.df_this_run = df_this_run
            return None


        for contour in self.contours:
            
            for dose1, dose2 in zip(contour['x'], contour['y']):            

                self._get_data_this_contour(*params, dose1, dose2)
               
                data = {"contour_level": contour['level'],
                            "dose1": dose1,
                            "dose2": dose2,
                            **self.this_contour_df}

                df_this_run = df_this_run.append(data, ignore_index=True)
        



        df_this_run['maxContEL'] = max(df_this_run['EL'])

        df_this_run['minDeltaRFB'] = min(df_this_run['delta_RFB'])        

        df_this_run['maxDeltaRFB'] = max(df_this_run['delta_RFB'])

        self.df_this_run = df_this_run






    def _get_data_this_contour(self, *params):

        self.ThisDoseConf = self._get_single_conf(*params, n_years=None)
                
        self.this_run_output = RunSingleTactic(self.fung_parms).run_single_tactic(self.ThisDoseConf)

        self.fy = self.this_run_output['failure_year']
        
        self._get_this_contour_df()






    def _get_this_contour_df(self):

        res_vecs = self.this_run_output['res_vec_dict']
        
        fy = int(self.fy)

        rf1 = res_vecs['f1'][fy]
        
        rf2 = res_vecs['f2'][fy]

        try:
            delt_rf = logit10_difference(rf1, rf2)
        except:
            delt_rf = "NA"

        econ = "NA"
        
        self.this_contour_df = dict(delta_RFB=delt_rf,
                                    EL=fy,
                                    econ=econ)



    
    def _get_single_conf(self, rfs1, rfs2, rfD, om_1, om_2,
                            delt_1, delt_2, sr_prop,
                            dose1, dose2,
                            n_years=None):

        if n_years is None:
            n_years = self.n_years

        Conf = SingleConfig(n_years, None, None, 
                                dose1, dose1, dose2, dose2,
                                primary_inoculum=self.inoc)
        
        Conf.sex_prop = sr_prop

        Conf.load_saved = self.config['load_saved']

        Conf.add_string()

        conf_str = _get_param_scan_conf_str(rfs1, rfs2, rfD, om_1, om_2,
                                                        delt_1, delt_2, Conf)
        
        Conf.config_string_img = conf_str

        Conf.config_string = _get_updated_conf_str(conf_str)

        return Conf






    def _plot_df(self):
        my_df = copy.copy(self.df_this_run)

        for i in my_df.contour_level.unique():
            df_plot = my_df[my_df['contour_level']==i]

            df_plot.plot(x="delta_RFB", y="EL", title=i, kind="scatter")
        
        plt.show()








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














class ParamScanRandRFB(ParamScanRand):
    def __init__(self, config) -> None:
        super().__init__(config)


    @staticmethod
    def _get_first_non_zero_element(vec):
        filtered = list(filter(lambda x: x>0, vec))
        return filtered[0]







    def _get_grid_output_this_run(self):
        self.my_grid_output = RunGrid(self.fung_parms).grid_of_tactics(self.MyGridConf)


    def _interpolate_contours(self, num=15):
        """
        Takes old contours and adds in some interpolated values
        
        Result is more densely populated list of values (approximately) 
        along the contour
        """
        xx = copy.copy(self.contours[0]['x'])
        yy = copy.copy(self.contours[0]['y'])

        nn = len(xx)
        
        fp = list(range(nn))
        
        start = float(fp[0])
        stop = float(fp[-1])

        to_find = np.linspace(start, stop, num)
        
        int_x = np.interp(to_find, fp, xx)
        int_y = np.interp(to_find, fp, yy)

        # include the old ones as well as interpolated ones

        cont_x = list(xx) + list(int_x)[1:-1]
        cont_y = list(yy) + list(int_y)[1:-1]

        self.contours[0]['x'] = cont_x
        self.contours[0]['y'] = cont_y




    def _find_contours_EqSel(self):
        self._find_EqSel_array()

        self._check_EqSel_validity()

        if not self.EqSelValid:
            self.contours = []
        else:    
            self.contours = _get_contours(self.EqSel_array, levels=[0.5], step=1)

            if self.contours and len(self.contours[0]['x'])<8:
                self._interpolate_contours(8)



    def _find_EqSel_array(self):

        data = self.my_grid_output

        out = np.ones(data['start_freqs']['SR'][:,:,0].shape)
        
        for i, j in itertools.product(range(out.shape[0]), 
                                        range(out.shape[1])):
            
            self.fy = data['FY'][i,j]

            if not self.fy>0:
                out[i,j] = None
            
            else:
                
                sr1 = data['start_freqs']['RS'][i,j,1]/data['start_freqs']['RS'][i,j,0]
                sr2 = data['start_freqs']['SR'][i,j,1]/data['start_freqs']['SR'][i,j,0]

                try:
                    out[i,j] = sr1/(sr1+sr2)
                except:
                    out[i,j] = None

        self.EqSel_array = out
        


    
    
    def _check_RFB_validity(self):
        self.ERFB_valid = (np.nanmax(self.RFB_array)>0
                        and np.nanmin(self.RFB_array)<0)



    def _check_EqSel_validity(self):
        self.EqSelValid = (np.nanmax(self.EqSel_array)>0.5
                        and np.nanmin(self.EqSel_array)<0.5)



    def _find_contours_RFB(self):
        """
        Find doses along the contours of constant first year yield.
        """

        self._find_RFB_array()

        self._check_RFB_validity()

        if not self.ERFB_valid:            
            self.contours = []
        else:
            self.contours = _get_contours(self.RFB_array, levels=[0], step=1)
            if self.contours and len(self.contours[0]['x'])<8:
                self._interpolate_contours(8)




    def _find_RFB_array(self):
        
        data = self.my_grid_output
        
        FYs = data['FY'][:,:]

        out = np.ones(FYs.shape)
        
        for i, j in itertools.product(range(FYs.shape[0]), range(FYs.shape[1])):
            
            self.fy = data['FY'][i,j]

            if not self.fy>0:
                out[i,j] = None
            
            else:
                
                rf1 = data['res_arrays']['f1'][i,j,int(self.fy)]
                rf2 = data['res_arrays']['f2'][i,j,int(self.fy)]

                try:
                    out[i,j] = logit10_difference(rf1, rf2)
                except:
                    out[i,j] = None
            

        self.RFB_array = out

        # plt.imshow(out)
        # plt.show()






    def _get_grid_output_df(self):

        output = self.my_grid_output

        minEqDoseELVec = [output['FY'][i, i] for i in range(output['FY'].shape[0])]

        data = dict(maxGridEL = np.amax(output['FY']),
                ERFB_Valid = self.ERFB_valid,
                GridMinRFB = np.nanmin(self.RFB_array),
                GridMaxRFB = np.nanmax(self.RFB_array),
                minEqDoseEL = self._get_first_non_zero_element(minEqDoseELVec),
                fullDoseEL = output['FY'][-1, -1],
                corner_00 = output['FY'][0, 0],
                corner_10 = output['FY'][-1, 0],
                corner_01 = output['FY'][0, -1],
                corner_11 = output['FY'][-1, -1])

        self.grid_output_df = pd.DataFrame([data])










    def _get_maxEL_EqSel(self):
        
        data = self.df_this_run

        if data.shape[0]>0:
            self.maxEqSelEL = max(data['EL'])
        else:
            self.maxEqSelEL = "NA"





    def _get_EqSel_columns(self, *this_run_parms):
        """
        adds on two columns to existing dataframe:
        - maxEqSelEL
        - EqSelValid
        """
        
        df_ERFB_data = copy.copy(self.df_this_run)

        self._find_contours_EqSel()

        self._get_data_this_conf_and_contrs(*this_run_parms)

        self._get_maxEL_EqSel()

        self._check_EqSel_validity()

        df_use = pd.DataFrame([dict(
                maxEqSelEL = self.maxEqSelEL,
                EqSelValid = self.EqSelValid
            )])
        
        self.df_RFB_and_EqSel = pd.concat([df_ERFB_data, df_use], axis=1)
        





    def _get_params_df(self):
        parms_df = pd.DataFrame([{**self.this_run_params}])

        run_index = int(parms_df.loc[0, 'run'])

        parms_df = parms_df.drop(["run"], axis=1)

        n_rows = self.df_RFB_and_EqSel.shape[0]

        run_df = pd.DataFrame([{"run": run_index}]*n_rows)

        out = pd.concat([parms_df, run_df], axis=1)

        return out







    def _get_df_to_append(self, *this_run_parms):
        
        self._get_grid_output_df()

        self._get_EqSel_columns(*this_run_parms)

        par_df = self._get_params_df()        

        out = pd.concat([par_df,
                            self.df_RFB_and_EqSel,
                            self.grid_output_df,
                            ],
                            axis=1)

        return out








    def run_param_scan_RFB(self, seed):

        df = pd.DataFrame()

        np.random.seed(seed)

        for run_index in tqdm(range(self.config["NIts"])):
            
            # initialise as empty
            self.df_this_run = pd.DataFrame()

            this_run_parms = self._get_random_pars()

            rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop = this_run_parms

            self._get_inoc(rfs1, rfs2, rfD)
            
            self._get_fung_parms(om_1, om_2, delt_1, delt_2)

            self.MyGridConf = self._get_grid_conf(*this_run_parms,
                            self.config["load_saved"],
                            n_years=self.n_years,
                            n_doses=self.config["grid_number"])
            
            self._get_grid_output_this_run()
            
            self._find_contours_RFB()

            self._get_this_run_params_rand(sr_prop, run_index)

            if self.ERFB_valid:
                self._get_data_this_conf_and_contrs(*this_run_parms)
        
            new_df = self._get_df_to_append(*this_run_parms)

            df = pd.concat([df, new_df], axis=0)
            
            


        return df











    def run(self, seed):
        """
        Run random scan over uniform dists
        """

        df = self.run_param_scan_RFB(seed)
        
        par_str = get_PS_rand_str(self.config)

        filename = f"./param_scan_cluster/outputs/rand/par_scan/seed={seed}_{par_str}.csv"
        
        print(f"Random Scan, saved as:\n {filename}")
        
        df.to_csv(filename, index=False)



















# post process fns

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

    def __init__(self, df, par_str):
        self.df = df
        self.par_str = par_str



    def _get_best_doses(self):
        self.best_doses = self.df.loc[self.df['EL']==self.df['maxEL']]





    def plot_df(self, run=0):
        filtered_df = self.df[self.df.run==run]
        for i in filtered_df.contour_level.unique():
            df_plot = filtered_df[filtered_df['contour_level']==i]

            df_plot.plot(x="delta_RFB", y="EL", title=i, kind="scatter")
        
        plt.show()

    
    
    
    def testing_RFB(self, run=0):
        
        my_df = self.df[self.df['run']==run]

        grouped_df = my_df.groupby(["run", "contour_level"])

        _, ax = plt.subplots(1)

        for key, this_df in grouped_df:
            this_df.plot(x="delta_RFB", y="EL", kind="line", ax=ax, label=key)

        plt.show()
        



    def process_best_doses(self):

        self._get_best_doses()

        # dfTest = self.best_doses

        # grouped_df = dfTest.groupby(["run", "contour_level"])


        



    def _print_runs_best_dose_not_minRFB(self):
        """
        If best EL is not for minimum absolute value for delta RFB 
        (but this is a bit arbitrary)
        """
        self.delt_failed = self.df[self.df["run"].isin(self.runs_failed)]

        print("Min RFB does not give the optimal EL:", self.delt_failed.run.unique())



    @staticmethod
    def _max_EL(data):
        return data["EL"].max()



    def check_max_EL_by_contour(self, type):
        grouped_df = self.df.groupby(["run", "contour_level"])

        grouped = grouped_df.pipe(self._max_EL)
        grouped.to_csv(f"./param_scan_cluster/outputs/{type}/analysis/maxEL_thisContour/{self.par_str}.csv")





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
        grouped_df = self.df.groupby(["run", "contour_level"])

        grouped_RFB_is_monotone = grouped_df.apply(self._monotone_RFB)
        
        grouped_effective_lives = grouped_df.apply(self._EL_by_RFB)
        
        grouped_RFBs = grouped_df.apply(self._rounded_RFB)
        
        combined = pd.concat([grouped_effective_lives, grouped_RFB_is_monotone, grouped_RFBs], axis=1)

        combined.columns = ["ELs", "monotoneRFB", "deltRFBs"]
        
        combined = combined.reset_index()

        self.monotoneRFB_df = combined
        
        combined.to_csv(f"./param_scan_cluster/outputs/{type}/analysis/monotoneRFB/{self.par_str}.csv")
        print(f"Monotone RFB works in all cases: {all(grouped_RFB_is_monotone)}")






    def get_params_for_specific_runs(self, which_runs):

        par_df = self.df

        par_df = par_df.drop(["Unnamed: 0", "EL",
            "contour_level", "delta_RFB"], axis=1)

        out = pd.DataFrame()
        for rr in which_runs:
            this_run = par_df[par_df["run"]==rr].iloc[0,:]
            out = out.append(this_run, ignore_index=True)
        
        return out






    def which_runs_worked_monotone_RFB(self, print_=False):

        failed = self.monotoneRFB_df[~self.monotoneRFB_df["monotoneRFB"]]
        succeeded = self.monotoneRFB_df[self.monotoneRFB_df["monotoneRFB"]]

        runs_that_failed = failed["run"].unique()
        runs_that_succeeded = succeeded["run"].unique()

        runs_that_succeeded = [e for e in runs_that_succeeded if e not in runs_that_failed]

        self.failed_pars = self.get_params_for_specific_runs(runs_that_failed)
        self.success_pars = self.get_params_for_specific_runs(runs_that_succeeded)

        n_fail = self.failed_pars.shape[0]
        n_success = self.success_pars.shape[0]
        self.failed_pars.to_csv(f"./param_scan_cluster/outputs/failed_pars/failed_{n_fail}_out_of_{n_success+n_fail}.csv")



        if print_:
            print("\n")
            print("These runs didn't have monotone EL vs RFB for every contour :(")
            print("\n")
            print(self.failed_pars)

            print("\n")
            print("These runs had monotone EL vs RFB for every contour :)")
            print("\n")
            print(self.success_pars)

    
    



    def which_runs_worked_max_cont(self):
        df = copy.copy(self.max_along_contour_df)

        failed = df[df['maxCont%']<100]

        runs_that_failed = failed["run"].unique()

        self.failed_pars = self.get_params_for_specific_runs(runs_that_failed)

        n_fail = self.failed_pars.shape[0]

        self.failed_pars.to_csv(f"./param_scan_cluster/outputs/failed_pars/failed_maxCont_{n_fail}.csv")


    @staticmethod
    def get_opt_df(data):
        if not data['maxContEL'].shape[0]:
            return None
        
        maxEL_this_run = float(list(data['maxContEL'])[0])

        df = data[data['EL']==maxEL_this_run]
        
        return df




    def min_opt_DS(self, data):
        
        df = self.get_opt_df(data)

        if df is None or not df.shape[0]:
            return "NA"

        return min(df['dose1']+df['dose2'])



    def max_opt_DS(self, data):
        
        df = self.get_opt_df(data)

        if df is None or not df.shape[0]:
            return "NA"
        
        return max(df['dose1']+df['dose2'])

    
    def _generate_max_along_contour_df(self):
        my_df = copy.copy(self.df)

        my_df.fillna(0)

        my_df['maxAlongContour'] = my_df['maxContEL'] >= my_df['maxGridEL']

        strats = ["maxCont", "minEqDose", "fullDose", "maxEqSel"]
        
        for string in strats:
            my_df[string + "%"] = 100*my_df[string + "EL"]/my_df["maxGridEL"]
        
        my_df = my_df.drop(['Unnamed: 0', 'econ'], axis=1)

        counting = my_df.groupby(["run"]).size()

        DS_df = pd.DataFrame()

        DS_df['min_dose_sums'] = my_df.groupby(["run"]).apply(lambda data: min(data['dose1'] + data['dose2']))
        DS_df['max_dose_sums'] = my_df.groupby(["run"]).apply(lambda data: max(data['dose1'] + data['dose2']))
        DS_df['min_opt_DS'] = my_df.groupby(["run"]).apply(self.min_opt_DS)
        DS_df['max_opt_DS'] = my_df.groupby(["run"]).apply(self.max_opt_DS)

        DS_df['min_DS_best'] = np.where(DS_df['min_dose_sums']==DS_df['min_opt_DS'], True, DS_df['min_dose_sums']-DS_df['min_opt_DS'])
        DS_df['max_DS_best'] = np.where(DS_df['max_dose_sums']==DS_df['max_opt_DS'], True, DS_df['max_dose_sums']-DS_df['max_opt_DS'])
        DS_df['max_or_min_DS_best'] = np.where((DS_df['min_dose_sums']==DS_df['min_opt_DS'])
                                        | (DS_df['max_dose_sums']==DS_df['max_opt_DS']),
                                            True, False)

        grouped = my_df.groupby(["run"]).first()

        grouped['count'] = counting
        
        df_out = grouped.reset_index()

        out = pd.concat([df_out, DS_df], axis=1)

        out = out.sort_values(['ERFB_Valid', 'EqSelValid', 'maxCont%', 'maxEqSel%'])

        return out




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



    # @staticmethod
    def analyse_max_contour_df(self):
        df = copy.copy(self.max_along_contour_df)

        df['min_corner'] = df[["corner_01", "corner_10"]].min(axis=1)

        # filter so only consider those with min corner > 0

        df = df[df['min_corner']>0]

        n_rows = df.shape[0]
        n_erfb = self.filtered_dataframe_outcome(df, "ERFB_Valid", "maxCont%")
        n_eq_sel = self.filtered_dataframe_outcome(df, "EqSelValid", "maxEqSel%")
        n_fd = self.filtered_dataframe_outcome(df, None, "fullDose%")
        n_min_eq_dose = self.filtered_dataframe_outcome(df, None, "minEqDose%")




    def analyse_failed(self):
        df_in = copy.copy(self.max_along_contour_df)

        df = df_in[df_in['maxCont%']<100]
        
        # df = df.assign(delta_ratio=lambda d: d['delta_1']/d['delta_2'])
        
        # df = df.assign(omega_ratio=lambda d: d['omega_1']/d['omega_2'])

        # df = df.assign(sing_res_ratio=lambda d: d['SR']/d['RS'])
        
        df = df.assign(diff_from_opt=lambda d: d['maxGridEL'] - d['maxContEL'])
        
        df['max_sing_res'] = df[["SR", "RS"]].max(axis=1)
        
        df['min_corner'] = df[["corner_01", "corner_10"]].min(axis=1)

        df = df.sort_values(by=['min_corner', 'max_sing_res'])

        print("\n")
        print("These runs failed:\n")
        
        print(df[['diff_from_opt',
                    # 'sr_prop',
                    # 'SR',
                    # 'RS',
                    # 'RR',
                    'run',
                    'count',
                    # "corner_00",
                    "corner_01",
                    "corner_10",
                    "min_corner",
                    # "corner_11",
                    # 'omega_1',
                    # 'omega_2',
                    'maxGridEL',
                    # 'sing_res_ratio',
                    # 'delta_ratio',
                    # 'omega_ratio',
                    'max_sing_res']])

        


    def get_maximum_along_contour_df(self):
        
        df_out = self._generate_max_along_contour_df()

        filename = f"./param_scan_cluster/outputs/rand/analysis/max_along_contour/df_{len(df_out)}.csv"

        print(f"Saving maximum along contour csv to: \n{filename}")

        df_out.to_csv(filename)

        self.max_along_contour_df = df_out


        

    @staticmethod
    def _analyse_hi_lo_df(df_in):
        df = df_in[df_in['ERFB_Valid']]
        
        #  = df.loc[:,'minEqDoseEL']/(df.loc[:,'fullDoseEL']+df.loc[:,'minEqDoseEL'])
        # df['min_vs_full'] = df['minEqDoseEL']/(df['fullDoseEL']+df['minEqDoseEL'])
        
        # X = df[['sr_prop']]
        # Y = df['min_vs_full']
        
        # print(df)

        # regr = linear_model.LinearRegression()
        # regr.fit(X, Y)
        # print(vars(regr))

        




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



    def _check_sing_res_smaller_than(self, pars, level=0.001):
        if pars.SR > level or pars.RS > level:
            return True
        else: 
            return False
    
    
    
    def re_run_failures(self, NDoses, failed_run_indices=None):

        if failed_run_indices is None:
            failed_run_indices = list(range(self.failed_pars.shape[0]))


        for ii in failed_run_indices:

            pars = self.failed_pars.iloc[int(ii),:]

            big_single_res = self._check_sing_res_smaller_than(pars)

            if big_single_res:
                print("single res are too big")
                continue
            
            print("\nRe-running run:", self.failed_pars.iloc[int(ii),:].run, "\n")
        
            PS = ParamScan()
            
            PS.NDoses = NDoses

            PS._get_inoc(pars["RS"], pars["SR"], pars["RR"])
            
            PS._get_fung_parms(pars["omega_1"], pars["omega_2"], 
                pars["delta_1"], pars["delta_2"])

            PS.GridConf = PS._get_grid_conf(pars["RS"], pars["SR"], pars["RR"],
                pars["omega_1"], pars["omega_2"], 
                pars["delta_1"], pars["delta_2"],
                pars["sr_prop"], True
                )
            

            output = RunGrid(PS.fung_parms).grid_of_tactics(PS.GridConf)

            conf_str = PS.GridConf.config_string_img
            


            dose_grid_heatmap(output, PS.GridConf, "FY", conf_str)
            
            # first_year_yield(output, PS.GridConf)

            eq_RFB_contours(output, PS.GridConf, title=f"Run={str(pars.run)}")






    def maxEL_by_contour(self):
        grouped_df = self.df.groupby(['run', 'EL'])
        for key, _ in grouped_df:
            print(grouped_df.get_group(key), "\n\n")
        
        pass



# End of Param Scan fns
import numpy as np
import pandas as pd

from model.config_classes import GridConfig, SingleConfig
from model.simulator import RunSingleTactic





class RandomPars:
    """
    Generates random parameters for model run
    
    Can return
    - pars
    - grid config
    - single config
    """

    def __init__(self, config, run_index) -> None:

        self.config = config

        self._run_index = run_index



    def find_pars(self):
        """
        returns rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, sr_prop
        """

        conf = self.config

        pars_are_valid = False

        while pars_are_valid is False:
        
            rfs1_power = np.random.uniform(low=conf["RF_single"][0], high=conf["RF_single"][1])
            rfs2_power = np.random.uniform(low=conf["RF_single"][0], high=conf["RF_single"][1]) 
            rfD_power = np.random.uniform(low=conf["RF_double"][0], high=conf["RF_double"][1])
            
            rfs1 = 10**(rfs1_power)
            rfs2 = 10**(rfs2_power)
            rfD = 10**(rfD_power)

            om_1 = np.random.uniform(low=conf["omega"][0], high=conf["omega"][1])
            om_2 = np.random.uniform(low=conf["omega"][0], high=conf["omega"][1])
            
            delt_1 = np.random.uniform(low=conf["decay_rate"][0], high=conf["decay_rate"][1])
            delt_2 = np.random.uniform(low=conf["decay_rate"][0], high=conf["decay_rate"][1])
            
            theta_1 = np.random.uniform(low=conf["theta"][0], high=conf["theta"][1])
            theta_2 = np.random.uniform(low=conf["theta"][0], high=conf["theta"][1])
            
            sr_prop = np.random.uniform(low=conf["SR"][0], high=conf["SR"][1])

            fungicide_params = self.get_fung_parms_dict(om_1, om_2, delt_1, delt_2, theta_1, theta_2)

            pathogen_pars = dict(rfs1=rfs1,
                                rfs2=rfs2,
                                rfD=rfD, 
                                sr_prop=sr_prop)

            pars_are_valid = self._check_validity(fungicide_params, pathogen_pars)

        
        self.get_inoc_dict(rfs1, rfs2, rfD)
        
        self.fung_parms = fungicide_params

        self.path_and_fung_pars = rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, theta_1, theta_2

        self.sr_prop = sr_prop








    def _check_validity(self, fungicide_params, pathogen_pars):

        self.get_inoc_dict(pathogen_pars['rfs1'],
                    pathogen_pars['rfs2'],
                    pathogen_pars['rfD'])
        
        
        corners = [[1,0], [0,1], [1,1]]

        for ii, jj in corners:
            this_dose_yield = self._get_yield_these_doses(fungicide_params,
                                                pathogen_pars, ii, jj)

            if this_dose_yield<=95:
                return False
        
        return True





    def _get_yield_these_doses(self, fungicide_params, pathogen_pars, d1, d2):

        conf_single = SingleConfig(1, None, None, d1, d1, d2, d2, primary_inoculum=self.inoc)
        
        conf_single.bs_sex_prop = pathogen_pars['sr_prop']

        conf_single.load_saved = False

        conf_single.add_string()

        this_run = RunSingleTactic(fungicide_params).run(conf_single)
        
        yield_out = this_run.yield_vec[0]

        return yield_out





    def get_inoc_dict(self, rfs1, rfs2, rfD):
        self.inoc = dict(RR = rfD,
                         RS = rfs1,
                         SR = rfs2,
                         SS = 1 - rfD - rfs1 - rfs2)

    

    def get_fung_parms_dict(self, omega_1, omega_2, delta_1, delta_2, theta_1, theta_2):
        return dict(omega_1 = omega_1,
                    omega_2 = omega_2,
                    theta_1 = theta_1,
                    theta_2 = theta_2,
                    delta_1 = delta_1,
                    delta_2 = delta_2)


    
    def get_all_parms_df(self):
        out = {**self.fung_parms,
                **self.inoc,
                "sr_prop": self.sr_prop,
                "run": self._run_index}
        
        self.par_df = pd.DataFrame([out])





    # single and grid config finders:

    def get_grid_conf(self, n_doses):
        
        conf = GridConfig(self.config['n_years'], None, None, n_doses,
                                    primary_inoculum=self.inoc)

        config_out = self._process_conf(conf)

        self.grid_conf = config_out

        return config_out





    def get_single_conf(self, dose1, dose2):

        conf = SingleConfig(self.config['n_years'], None, None,
                                dose1, dose1, dose2, dose2,
                                primary_inoculum=self.inoc)
        
        config_out = self._process_conf(conf)

        self.sing_conf = config_out

        return config_out




    def _process_conf(self, conf):

        conf.bs_sex_prop = self.sr_prop

        conf.load_saved = self.config['load_saved']
        
        conf.save = self.config['save']

        conf.add_string()

        config_out = self._update_par_scan_conf_str(conf)

        return config_out



    
    def _update_par_scan_conf_str(self, conf):
        
        rfs1, rfs2, rfD, om_1, om_2, delt_1, delt_2, thet_1, thet_2 = self.path_and_fung_pars

        conf_str = conf.config_string_img
        conf_str = conf_str.replace("grid", "param_scan")

        par_str = (f"_fung_pars={round(om_1,3)},{round(om_2,3)},"
                f"{round(delt_1,3)},{round(delt_2,3)},"
                f"{round(thet_1,2)},{round(thet_2,2)},"
                f"rf1={round(rfs1,10)},rf2={round(rfs2,10)},rfd={round(rfD,15)}")
        
        par_str = par_str.replace(".", ",")
        conf_str = conf_str.replace(".png", par_str + ".png")
        
        conf.config_string_img = conf_str
        
        saved_run_str = conf_str.replace(".png", ".pickle")
        
        conf.config_string = saved_run_str.replace("figures/", "saved_runs/")

        return conf




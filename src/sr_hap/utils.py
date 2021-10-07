from math import floor, sqrt
import numpy as np
import pandas as pd


from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig, SingleConfig
from model.params import PARAMS


def get_haploid_outputs_res(n_variants, n_doses, double_freq_factor_lowest, index, bss):
    outputs = []
    
    dfp = get_sr_scan_params_res(n_variants, double_freq_factor_lowest, index)
    
    print(dfp)

    for bs in bss:

        output = get_this_run_output(n_doses, bs, dfp)

        outputs.append(output)
    return outputs



def get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, index, bss):
    outputs = []
    
    dfp = get_sr_scan_params_app(n_variants, double_freq_factor_lowest, index)
    
    print(dfp)

    for bs in bss:

        output = get_this_run_output(n_doses, bs, dfp)

        outputs.append(output)
    return outputs






def get_this_run_output(n_doses, bs, dfp):

    data = dfp.iloc[int(0),:]

    conf = GridConfig(30, None, None, n_doses)
    
    conf.primary_inoculum = dict(
        RR = data["RR"],
        RS = data["RS"],
        SR = data["SR"],
        SS = data["SS"],
        )

    fcide_pars = dict(
        omega_1 = data["omega_1"],
        omega_2 = data["omega_2"],
        theta_1 = data["theta_1"],
        theta_2 = data["theta_2"],
        delta_1 = data["delta_1"],
        delta_2 = data["delta_2"],
        )

    conf.load_saved = True

    conf.bs_sex_prop = bs
    conf.add_string()

    output = RunGrid(fcide_pars).run(conf)

    return output



def get_rfd(x, y):
    B = 1 - x - y
    C = x*y
    out = (B - sqrt(B**2 -4*C))/2
    return out







def get_sr_scan_params_app(n_variants, double_freq_factor_lowest, index):

    power_of_10 = floor(index/n_variants)

    double_freq_factor = double_freq_factor_lowest * 10**(power_of_10)

    index = index % n_variants

    np.random.seed(index)

    valid = False

    while not valid:

        rf1 = 10**(np.random.uniform(-8, -4))
        rf2 = 10**(np.random.uniform(-8, -4))
        
        rfd = get_rfd(rf1, rf2)
        rfd = double_freq_factor*rfd

        om1 = np.random.uniform(0.4,1)
        om2 = np.random.uniform(0.4,1)
        
        thet1 = np.random.uniform(4,12)
        thet2 = np.random.uniform(4,12)
        
        delta_factor1 = np.random.uniform(1/3,3)
        delta_factor2 = np.random.uniform(1/3,3)



        params = dict(
            RR=rfd,
            RS=rf1,
            SR=rf2,
            SS=1-rf1-rf2-rfd,
            omega_1 = om1,
            omega_2 = om2,
            theta_1 = thet1,
            theta_2 = thet2,
            delta_1 = PARAMS.delta_1*delta_factor1,
            delta_2 = PARAMS.delta_2*delta_factor2,
            double_freq_factor = double_freq_factor,
            )
        
        conf = SingleConfig(1, None, None,
                            1, 1, 1, 1,
                            primary_inoculum=None
                            )

        conf.primary_inoculum = dict(
            RR = params["RR"],
            RS = params["RS"],
            SR = params["SR"],
            SS = params["SS"],
            )

        fcide_pars = dict(
            omega_1 = params["omega_1"],
            omega_2 = params["omega_2"],
            theta_1 = params["theta_1"],
            theta_2 = params["theta_2"],
            delta_1 = params["delta_1"],
            delta_2 = params["delta_2"],
            )
        
        conf.load_saved = False
        conf.add_string()

        yld = RunSingleTactic(fcide_pars).run(conf).yield_vec[0]
        
        if yld>95:
            valid = True
            par_df = pd.DataFrame([params])

    return par_df







def get_sr_scan_params_res(n_its, double_freq_factor, index):

    np.random.seed(index)

    # double_freq_factor = double_freq_factors[int(floor(index/n_its))]

    valid = False

    while not valid:

        rf1 = 10**(np.random.uniform(-8, -4))
        rf2 = 10**(np.random.uniform(-8, -4))
        
        rfd = get_rfd(rf1, rf2)
        rfd = double_freq_factor*rfd

        om1 = np.random.uniform(0.4,1)
        om2 = np.random.uniform(0.4,1)
        
        thet1 = np.random.uniform(4,12)
        thet2 = np.random.uniform(4,12)
        
        delta_factor1 = np.random.uniform(1/3,3)
        delta_factor2 = np.random.uniform(1/3,3)



        params = dict(
            RR=rfd,
            RS=rf1,
            SR=rf2,
            SS=1-rf1-rf2-rfd,
            omega_1 = om1,
            omega_2 = om2,
            theta_1 = thet1,
            theta_2 = thet2,
            delta_1 = PARAMS.delta_1*delta_factor1,
            delta_2 = PARAMS.delta_2*delta_factor2,
            double_freq_factor = double_freq_factor,
            )
        
        conf = SingleConfig(1, None, None,
                            1, 1, 1, 1,
                            primary_inoculum=None
                            )

        conf.primary_inoculum = dict(
            RR = params["RR"],
            RS = params["RS"],
            SR = params["SR"],
            SS = params["SS"],
            )

        fcide_pars = dict(
            omega_1 = params["omega_1"],
            omega_2 = params["omega_2"],
            theta_1 = params["theta_1"],
            theta_2 = params["theta_2"],
            delta_1 = params["delta_1"],
            delta_2 = params["delta_2"],
            )
        
        conf.load_saved = False
        conf.add_string()

        yld = RunSingleTactic(fcide_pars).run(conf).yield_vec[0]
        
        if yld>95:
            valid = True
            par_df = pd.DataFrame([params])

    return par_df
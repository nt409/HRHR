from math import floor, sqrt
import numpy as np
import pandas as pd
import sys
# import itertools
# import os
# from tqdm import tqdm
import plotly.graph_objects as go

from model.simulator import RunGrid, RunSingleTactic
from model.config_classes import GridConfig, SingleConfig
from model.params import PARAMS


def get_sr_scan_df(n_its, n_sex_props, n_doses, double_freq_factors, index):
    
    dfp = get_sr_scan_params(n_its, double_freq_factors, index)

    sex_props = np.linspace(0, 1, n_sex_props)

    df = pd.DataFrame()



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

    conf.load_saved = False

    for bs in sex_props:

        conf.bs_sex_prop = bs
        conf.add_string()

        output = RunGrid(fcide_pars).run(conf)

        data["maxEL"] = np.amax(output.FY)
        data["bs_sex_prop"] = bs
        data["run"] = index
        
        df = df.append(data, ignore_index=True)


    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_hap/outputs/single/df_{n_its}_{n_sex_props}_{n_doses}_{dff_str}_{index}.csv"
    print(f"Saving df to: {filename}")
    df.to_csv(filename, index=False)
    
    return df







def get_sr_scan_params(n_its, double_freq_factors, index):

    np.random.seed(index)

    double_freq_factor = double_freq_factors[int(floor(index/n_its))]

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



def get_rfd(x, y):
    B = 1 - x - y
    C = x*y
    out = (B - sqrt(B**2 -4*C))/2
    return out





def plot_df(df):
    trcs = []

    for rr in range(n_its):
        data_use = df.loc[df["run"]==rr, :]
        xx = data_use["bs_sex_prop"]
        yy = data_use["maxEL"]
        trc = dict(x=xx, y=yy, name=f"run={rr}")

        trcs.append(trc)

    fig = go.Figure(data=trcs, layout=dict(template="simple_white"))
    fig.show()




if __name__=="__main__":
    n_its = 10
    n_sex_ps = 11
    n_doses = 21
    # n_doses = 6

    if len(sys.argv)!=2:
        raise Exception("Supply one argument: a run index")

    index = int(sys.argv[1])
    
    double_freq_factors = [1e-5, 1, 1e5]
    df = get_sr_scan_df(n_its, n_sex_ps, n_doses, double_freq_factors, index)

    print(df)

    # plot_df(df)

        




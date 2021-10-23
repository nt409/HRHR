import copy
import sys
import numpy as np
import pandas as pd


from model.simulator import RunGrid
from model.config_classes import GridConfig
from alt_scan.utils import get_alt_scan_params_rand


def main(n_doses, alt_strat, index):
   
    primary_inoc, fcide_parms = get_alt_scan_params_rand(index)

    model_run = RunGrid(fcide_parms)
    



    # MIXTURE
    # asex, mixture
    conf_grid_a_mix = GridConfig(40, None, None, n_doses)
    conf_grid_a_mix.load_saved = False
    
    conf_grid_a_mix.primary_inoculum = primary_inoc
    conf_grid_a_mix.add_string()

    asex_mix = model_run.run(conf_grid_a_mix)

    # sex, mixture
    conf_grid_s_mix = copy.copy(conf_grid_a_mix)
    conf_grid_s_mix.bs_sex_prop = 1
    conf_grid_s_mix.add_string()

    sex_mix = model_run.run(conf_grid_s_mix)




    # ALTERNATION
    # asex, alt
    conf_grid_a_alt = copy.copy(conf_grid_a_mix)
    conf_grid_a_alt.strategy = alt_strat
    conf_grid_a_alt.add_string()

    asex_alt = model_run.run(conf_grid_a_alt)

    # sex, alt
    conf_grid_s_alt = copy.copy(conf_grid_a_alt)
    conf_grid_s_alt.bs_sex_prop = 1
    conf_grid_s_alt.add_string()
    
    sex_alt = model_run.run(conf_grid_s_alt)




    # ADD TO DATAFRAME
    data = dict(
        RR =  primary_inoc["RR"],
        RS =  primary_inoc["RS"],
        SR =  primary_inoc["SR"],
        omega1 = fcide_parms["omega_1"],
        omega2 = fcide_parms["omega_2"],
        theta1 = fcide_parms["theta_1"],
        theta2 = fcide_parms["theta_2"],
        asex_alt = np.amax(asex_alt.FY),
        asex_mix = np.amax(asex_mix.FY),
        sex_alt = np.amax(sex_alt.FY),
        sex_mix = np.amax(sex_mix.FY),
        )
    
    df = pd.DataFrame([data])



    filename = f"./alt_scan/outputs/single/out_rand_{n_doses}_{alt_strat}_{index}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename)

    return df




if __name__=="__main__":

    if len(sys.argv)!=2:
        raise Exception("Supply one argument: the run index")

    index = int(sys.argv[1])

    df = main(n_doses=21, alt_strat="alt_21", index=index)
    df = main(n_doses=21, alt_strat="alt_12", index=index)



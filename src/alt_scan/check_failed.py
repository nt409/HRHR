import copy
import numpy as np
import pandas as pd

from model.params import PARAMS
from model.simulator import RunGrid
from model.config_classes import GridConfig




def main(rf1, rf2, om1, om2, thet1, thet2, alt_strat, n_doses):
    primary_inoc = dict(RR=rf1*rf2, RS=rf1, SR=rf2, SS=1-rf1-rf2-rf1*rf2)

    fcide_parms = dict(
        omega_1 = om1,
        omega_2 = om2,
        theta_1 = thet1,
        theta_2 = thet2,
        delta_1 = PARAMS.delta_1,
        delta_2 = PARAMS.delta_2,
        )

    model_run = RunGrid(fcide_parms)

    conf_grid_a_mix = GridConfig(40, None, None, n_doses)
    conf_grid_a_mix.load_saved = False
    conf_grid_a_mix.primary_inoculum = primary_inoc
    conf_grid_a_mix.add_string()
    asex_mix = model_run.run(conf_grid_a_mix)

    # print(np.amax(asex_mix.FY))

    conf_grid_a_alt = copy.copy(conf_grid_a_mix)
    conf_grid_a_alt.strategy = alt_strat
    conf_grid_a_alt.add_string()

    asex_alt = model_run.run(conf_grid_a_alt)

    return dict(
        asex_alt = np.amax(asex_alt.FY),
        asex_mix = np.amax(asex_mix.FY),
        )


if __name__=="__main__":
    
    alt_strat = "alt_12"
    rf1 = 1e-10
    rf2 = 1e-5

    om1 = 1
    om2 = 1

    thet1 = 10
    thet2 = 9.6

    n_doses = 21

    out = main(rf1, rf2, om1, om2, thet1, thet2, alt_strat, n_doses)

    print(out)

    df = pd.DataFrame([out])
    filename = f"./alt_scan/failed_run_{alt_strat}_{n_doses}_{rf1}_{rf2}_{om1}_{om2}_{thet1}_{thet2}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)
    
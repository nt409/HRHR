import copy
import numpy as np
import pandas as pd

# from model.params import PARAMS
from model.simulator import RunGrid
from model.config_classes import GridConfig
from alt_scan.utils import get_alt_scan_params




def main(n_doses, n_its, alt_strat, index):

    primary_inoc, fcide_parms, _, _, _ = get_alt_scan_params(n_its, index)

    model_run = RunGrid(fcide_parms)

    # run mix
    conf_grid_a_mix = GridConfig(40, None, None, n_doses)
    conf_grid_a_mix.load_saved = False
    conf_grid_a_mix.primary_inoculum = primary_inoc
    conf_grid_a_mix.add_string()
    asex_mix = model_run.run(conf_grid_a_mix)


    # run alt
    conf_grid_a_alt = copy.copy(conf_grid_a_mix)
    conf_grid_a_alt.strategy = alt_strat
    conf_grid_a_alt.add_string()

    asex_alt = model_run.run(conf_grid_a_alt)

    return dict(
        asex_alt = np.amax(asex_alt.FY),
        asex_mix = np.amax(asex_mix.FY),
        )


if __name__=="__main__":
    
    n_doses = 21
    n_its = 5
    alt_strat = "alt_12"
    index = 5

    # rf1 = 1e-10
    # rf2 = 1e-5

    # om1 = 1
    # om2 = 1

    # thet1 = 10
    # thet2 = 9.6

    out = main(n_doses, n_its, alt_strat, index)

    print(out)


    df = pd.DataFrame([out])
    filename = f"./alt_scan/outputs/failed_runs/{alt_strat}_{n_doses}_{n_its}_{index}.csv"
    print(f"saving df to: {filename}")
    df.to_csv(filename, index=False)
    
from model.simulator import RunGrid
from model.config_classes import GridConfig
from sr_scan.sf_ratio import get_sf_ratio_data
from model.params import PARAMS
from plotting.paper_figs import SRPlot
from sr_scan.get_df import get_sr_grid_df


n_doses = 21
n_sp = 41

# n_doses = 11
# n_sp = 11

# n_doses = 3
# n_sp = 3


filestr = f"Nd={n_doses}_Nsp={n_sp}"

def_df = get_sr_grid_df(n_doses, n_sp, load = False)

fcide_pars = dict(theta_1=8, theta_2=8,
                    omega_1=0.85, omega_2=0.85,
                    delta_1=PARAMS.delta_1, delta_2=PARAMS.delta_2)

low_df = get_sr_grid_df(n_doses, n_sp, fcide_pars, load = False)

    
sf_ratio_data = get_sf_ratio_data(n_doses, load = False)


rf1s, rf2s = 1e-5, 1e-5
primary_inoc_same = dict(RR=rf1s*rf2s, RS=rf1s, SR=rf2s, SS=1-rf1s-rf2s-rf1s*rf2s)

# same RFs, same fung
conf_grid = GridConfig(30, None, None, n_doses)
conf_grid.load_saved = False
conf_grid.primary_inoculum = primary_inoc_same
conf_grid.bs_sex_prop = 1
conf_grid.add_string()

grid_output = RunGrid().run(conf_grid)



SRPlot(def_df, low_df, sf_ratio_data, grid_output, filestr)

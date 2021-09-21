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
def_df = get_sr_grid_df(n_doses, n_sp)

fcide_pars = dict(theta_1=8, theta_2=8,
                    omega_1=0.85, omega_2=0.85,
                    delta_1=PARAMS.delta_1, delta_2=PARAMS.delta_2)

low_df = get_sr_grid_df(n_doses, n_sp, fcide_pars)

    
sf_ratio_data = get_sf_ratio_data(n_doses)

SRPlot(def_df, low_df, sf_ratio_data, filestr)

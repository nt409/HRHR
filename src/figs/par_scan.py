"""
Saves plots like this one:
../outputs/figures/paper_figs/par_scan_RF_bd=-10,-3_RF_bd=-15,-4_ome_bd=0,1_the
_bd=0,12_dec_bd=0037,0333_SR_bd=0,1_gri=51_cont=100_year=45_iter=5.png
"""

from plotting.paper_figs import ParamScanPlotMeVsHobb

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess


if __name__ == "__main__":
    par_str = config_rand['par_str']
    data = PostProcess(par_str).processed_df
    ParamScanPlotMeVsHobb(data, f"{par_str}.png")

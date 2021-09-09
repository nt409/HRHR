from plotting.paper_figs import ParamScanPlotMeVsHobb  

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess

def get_param_data(par_str):
    PP = PostProcess(par_str)
    # data = PP.processed_df
    # par_data = PP.par_df
    # return data #, par_data
    return PP.processed_df


par_str = config_rand['par_str']
# data = get_param_data(par_str)
data = PostProcess(par_str).processed_df
# data, par_data = get_param_data(par_str)
ParamScanPlotMeVsHobb(data, f"{par_str}.png")



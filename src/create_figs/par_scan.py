from plotting.paper_figs import ParamScanPlotMeVsHobb, ParamScanPlotHighLowDose  

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess

def get_param_data(par_str):
    PP = PostProcess(par_str)
    PP.get_maximum_along_contour_df()
    data = PP.max_along_contour_df
    return data


par_str = config_rand['par_str']
data = get_param_data(par_str)
ParamScanPlotMeVsHobb(data, f"{par_str}.png")
# ParamScanPlotHighLowDose(data, f"{par_str}.png")



from plotting.paper_figs import ParamScanPlotMeVsHobb  

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess



if __name__=="__main__":
    par_str = config_rand['par_str']
    data = PostProcess(par_str).processed_df
    ParamScanPlotMeVsHobb(data, f"{par_str}.png")



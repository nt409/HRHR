"""
Check why the single failed run actually failed
"""

from model.strategy_arrays import EqualResFreqBreakdownArray
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess




def get_data_for_appendix(grid_output):
    
    # one year before optimal failure year
    yy = np.amax(grid_output.FY) - 1

    yld = grid_output.yield_array[:, :, yy]
    best_yield = np.amax(yld)

    x = np.where(yld==best_yield)
    
    best_SR = grid_output.end_freqs_DA['SR'][x[0][0], x[1][0], yy]
    best_RS = grid_output.end_freqs_DA['RS'][x[0][0], x[1][0], yy]
    best_RR = grid_output.end_freqs_DA['RR'][x[0][0], x[1][0], yy]

    print(
          f"\n best dose RR: {best_RR}"
          f"\n best dose RS: {best_RS}"
          f"\n best dose SR: {best_SR}"
          f"\n best dose SS: {1 - best_SR - best_RS - best_RR}"
          f"\n best dose yield: {best_yield}"
        )



    RFB = EqualResFreqBreakdownArray(grid_output).array

    small_RFB = abs(RFB)<0.01

    ind_RFB = np.where(small_RFB)
    y = yld[small_RFB]

    which_ind = np.where(y==np.amax(y))[0][0]

    x1, x2 = ind_RFB[0][which_ind], ind_RFB[1][which_ind]

    RFB_SR = grid_output.end_freqs_DA['SR'][x1, x2, yy]
    RFB_RS = grid_output.end_freqs_DA['RS'][x1, x2, yy]
    RFB_RR = grid_output.end_freqs_DA['RR'][x1, x2, yy]
    RFB_yld = grid_output.yield_array[x1, x2, yy]
    
    
    print(
          f"\n ERFB RR: {RFB_RR}"
          f"\n ERFB SR: {RFB_SR}"
          f"\n ERFB RS: {RFB_RS}"
          f"\n ERFB SS: {1 - RFB_RR - RFB_RS - RFB_SR}"
          f"\n ERFB yield: {RFB_yld}"
        )


    print(f"year is {yy} or {yy+1}")



if __name__=="__main__":
    PP = PostProcess(config_rand['par_str'])

    grid_output = PP.re_run_grid(NDoses=201, run_indices=[183], plot=False)

    get_data_for_appendix(grid_output)

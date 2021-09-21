"""
Check why the single failed run actually failed
"""

from model.strategy_arrays import EqualResFreqBreakdownArray
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

from param_scan.fns.config import config_rand
from param_scan.fns.post_process import PostProcess



PP = PostProcess(config_rand['par_str'])

grid_output = PP.re_run_grid(NDoses=201, run_indices=[117], plot=False)

def heatmap(z, colorscale=None):
    z = np.transpose(z)
    trc = go.Heatmap(x=list(range(z.shape[0])), y=list(range(z.shape[1])), z=z, 
            colorscale=colorscale)
    fig = go.Figure(data=[trc])
    fig.show()


def yield_by_year():
    yld = grid_output.yield_array[80, 23, :13]

    xx = list(range(len(yld)))

    fig = go.Figure(
                data=[dict(x=xx, y=yld)]
                    )
    fig.show()
    return None


def get_yield_cs():
    cs = []
    pltly_clr_scale = list(px.colors.sequential.Inferno)

    my_grey = "rgb(0.4,0.4,0.4)"
    pal = [my_grey, my_grey] + pltly_clr_scale


    lim1 = 94/np.amax(yld)
    vals = [0, lim1] + list(np.linspace(lim1, 1, len(pal)-2))

    for val, col in zip(vals, pal):
        cs.append([val, col])
    return cs

def get_RFB_cs():
    cs = []

    my_grey = "rgb(0.4,0.4,0.4)"
    pal = [my_grey, my_grey, "red", "red"]


    lim1 = -np.nanmin(RFB)/(np.nanmax(RFB)-np.nanmin(RFB))
    vals = [0, lim1, lim1, 1]

    for val, col in zip(vals, pal):
        cs.append([val, col])
    return cs

def get_data_for_app():
    yld = grid_output.yield_array[:, :, 11]
    best_yield = np.amax(yld)

    x = np.where(yld==best_yield)
    
    best_SR = grid_output.end_freqs_DA['SR'][x[0][0], x[1][0], 11]
    best_RS = grid_output.end_freqs_DA['RS'][x[0][0], x[1][0], 11]
    best_RR = grid_output.end_freqs_DA['RR'][x[0][0], x[1][0], 11]

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

    RFB_SR = grid_output.end_freqs_DA['SR'][x1, x2, 11]
    RFB_RS = grid_output.end_freqs_DA['RS'][x1, x2, 11]
    RFB_RR = grid_output.end_freqs_DA['RR'][x1, x2, 11]
    RFB_yld = grid_output.yield_array[x1, x2, 11]
    
    
    print(
          f"\n ERFB RR: {RFB_RR}"
          f"\n ERFB SR: {RFB_SR}"
          f"\n ERFB RS: {RFB_RS}"
          f"\n ERFB SS: {1 - RFB_RR - RFB_RS - RFB_SR}"
          f"\n ERFB yield: {RFB_yld}"
        )




if __name__=="__main__":

    get_data_for_app()

        
    # optimal
    # sr = grid_output.end_freqs_DA['SR'][74:84, 20:26, 12]
    # rs = grid_output.end_freqs_DA['RS'][74:84, 20:26, 12]
    # rr = grid_output.end_freqs_DA['RR'][74:84, 20:26, 12]

    # optimal and 1 year off
    sr = grid_output.end_freqs_DA['SR'][64:152, 12:65, 11]
    rs = grid_output.end_freqs_DA['RS'][64:152, 12:65, 11]
    rr = grid_output.end_freqs_DA['RR'][64:152, 12:65, 11]
    yld = grid_output.yield_array[64:152, 12:65, 11]
    fy = grid_output.FY[64:152, 12:65]
    RFB = EqualResFreqBreakdownArray(grid_output).array[64:152, 12:65]


    # cs = get_RFB_cs()
    # heatmap(RFB, cs)
    # heatmap(sr)
    # heatmap(rs)
    # heatmap(rr)
    # cs = get_yield_cs()
    # heatmap(yld, cs)
    # heatmap(fy)


    # x = np.where(grid_output.FY == 12)
    # print(min(x[0]))
    # print(max(x[0]))
    # print(min(x[1]))
    # print(max(x[1]))

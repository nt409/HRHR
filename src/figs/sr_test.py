import plotly.graph_objects as go
import numpy as np

from model.simulator import RunGrid
from model.config_classes import GridConfig
from plotting.paper_figs import DoseSpaceScenarioSingle


n_doses = 6


filestr = f"Nd={n_doses}_.png"

rf1s, rf2s = 1e-5, 1e-5
rfds = np.logspace(-15,-5,11)

for rfd in rfds:

    primary_inoc_same = dict(RR=rfd, RS=rf1s, SR=rf2s, SS=1-rf1s-rf2s-rfd)

    # same RFs, same fung
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = False
    conf_grid.primary_inoculum = primary_inoc_same
    conf_grid.bs_sex_prop = 1
    conf_grid.add_string()

    grid_output = RunGrid().run(conf_grid)

    DoseSpaceScenarioSingle(grid_output, filestr)

    if np.amax(grid_output.FY)==grid_output.FY[-1,-1]:
        print(f"it worked, rfd={rfd}")
        DoseSpaceScenarioSingle(grid_output, filestr)






def linegraphs():
    # FY = grid_output.FY[-1,-1]
    FY = 10

    y_RR = grid_output.end_freqs_DA['RR'][-1, -1, :int(FY)]
    y_SR = grid_output.end_freqs_DA['SR'][-1, -1, :int(FY)]
    y_RS = grid_output.end_freqs_DA['RS'][-1, -1, :int(FY)]
    y_SS = grid_output.end_freqs_DA['SS'][-1, -1, :int(FY)]

    y_RR_s = grid_output.start_freqs_DA['RR'][-1, -1, :int(FY)]
    y_SR_s = grid_output.start_freqs_DA['SR'][-1, -1, :int(FY)]
    y_RS_s = grid_output.start_freqs_DA['RS'][-1, -1, :int(FY)]
    y_SS_s = grid_output.start_freqs_DA['SS'][-1, -1, :int(FY)]

    xx = list(range(len(y_RR)))

    trcs = [
        dict(x=xx, y=y_RR,  line=dict(dash="solid"), name="RR end"),
        dict(x=xx, y=y_SR,  line=dict(dash="solid"), name="SR end"),
        dict(x=xx, y=y_RS,  line=dict(dash="solid"), name="RS end"),
        dict(x=xx, y=y_SS,  line=dict(dash="solid"), name="SS end"),
        dict(x=xx, y=y_RR_s, line=dict(dash="dash"), name="RR start"),
        dict(x=xx, y=y_SR_s, line=dict(dash="dash"), name="SR start"),
        dict(x=xx, y=y_RS_s, line=dict(dash="dash"), name="RS start"),
        dict(x=xx, y=y_SS_s, line=dict(dash="dash"), name="SS start"),
        ]

    fig = go.Figure(data=trcs, layout=dict(template="simple_white", title=f"sex={conf_grid.bs_sex_prop}"))
    fig.show()




# DoseSpaceScenarioSingle(grid_output, filestr)

import plotly.graph_objects as go
import numpy as np

from .params import PARAMS

# TOC
# Utility functions
# Single Tactic
# Grid of tactics

attrs_dict = {
    str(PARAMS.S_ind): dict(name='Susceptible', colour='green'),
    str(PARAMS.ER_ind): dict(name='Latent (RR)', colour='rgb(0,0,150)'),
    str(PARAMS.ERS_ind): dict(name='Latent (RS)', colour='rgb(0,0,180)'),
    str(PARAMS.ESR_ind): dict(name='Latent (SR)', colour='rgb(0,0,210)'),
    str(PARAMS.ES_ind): dict(name='Latent (SS)', colour='rgb(0,0,240)'),
    str(PARAMS.IR_ind): dict(name='Infected (RR)', colour='rgb(150,0,0)'),
    str(PARAMS.IRS_ind): dict(name='Infected (RS)', colour='rgb(180,0,0)'),
    str(PARAMS.ISR_ind): dict(name='Infected (SR)', colour='rgb(240,0,0)'),
    str(PARAMS.IS_ind): dict(name='Infected (SS)', colour='rgb(240,0,0)'),
    str(PARAMS.R_ind): dict(name='Removed', colour='black'),
    str(PARAMS.PR_ind): dict(name='Primary inoc. (RR)', colour='rgb(150,0,150)'),
    str(PARAMS.PRS_ind): dict(name='Primary inoc. (RS)', colour='rgb(180,0,180)'),
    str(PARAMS.PSR_ind): dict(name='Primary inoc. (SR)', colour='rgb(210,0,210)'),
    str(PARAMS.PS_ind): dict(name='Primary inoc. (SS)', colour='rgb(240,0,240)'),
    str(PARAMS.Fung1_ind): dict(name='Fung. 1 conc.', colour='rgb(150,200,0)'),
    str(PARAMS.Fung2_ind): dict(name='Fung. 2 conc.', colour='rgb(200,150,0)'),
    }


#----------------------------------------------------------------------------------------------
# * Utility functions

def standard_layout(legend_on):
    return go.Layout(
            font = dict(size=22),
            template="plotly_white",
            width=1400,
            height=700,
            showlegend=legend_on,
            )


TITLE_MAP = dict(
    LTY = "Lifetime yield",
    TY = "Total yield",
    FY = "Failure year",
)
# End of utility functions

#----------------------------------------------------------------------------------------------
# * Single Tactic

def yield_by_year(data):
    traces = []
    
    y = data['yield_vec']
    x = list(range(1,1+len(y)))

    line = go.Scatter(
        x = x,
        y = y,
    )

    traces.append(line)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_xaxes(title="Year")
    fig.update_yaxes(title="Yield")

    return fig




def res_freqs_single_t_plot(data):
    traces = []

    titles = dict(
        f1 = "Fungicide 1",
        f2 = "Fungicide 2"
    )
    
    for key in ['f1', 'f2']:
        y = data['res_vec_dict'][key]
        x = list(range(len(y)))

        line = go.Scatter(
            x = x,
            y = y,
            name = titles[key]
        )

        traces.append(line)

    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_xaxes(title="Year")
    fig.update_yaxes(title="Resistant frequency")

    return fig



def single_year_plot(data):
    traces = []
    
    for ind in range(16):
        y = data['sol_array'][:, ind, 0]
        x = data['t_vec']

        line = go.Scatter(
            x = x,
            y = y,
            name = attrs_dict[str(ind)]['name'],
            line = dict(color=attrs_dict[str(ind)]['colour'])
        )

        traces.append(line)

    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_xaxes(title="Time (degree-days)")
    fig.update_yaxes(title="Amount")

    return fig



# End of single Tactic


#----------------------------------------------------------------------------------------------
# * Grid of tactics

def dose_grid_heatmap(data, Config, to_plot):
    traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])

    trace = go.Contour(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = TITLE_MAP[to_plot], # title here
            titleside = 'right',
            # titlefont=dict(
            #     size=14,
            #     family='Arial, sans-serif')
        )
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=660, height=600)

    fig.update_xaxes(title="Dose (fungicide. 1)")
    fig.update_yaxes(title="Dose (fungicide. 2)")

    return fig


# End of grid of tactics
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from math import log2, floor

from .params import PARAMS

# TOC
# Utility functions
# Single Tactic
# Changing dose
# Grid of tactics
# Dose space

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

LABEL_COLOR = "rgb(110,110,110)"

#----------------------------------------------------------------------------------------------
# * Utility functions

def standard_layout(legend_on):
    return go.Layout(
            font = dict(size=22),
            template="plotly_white",
            width=1400,
            height=700,
            showlegend=legend_on,
            xaxis=dict(showgrid=False),
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



def single_year_plot(data, indices):
    traces = []
    
    for ind in indices:
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
# * Changing dose

def SR_by_dose_plot(data):
    rows_list = []

    for key in data.keys():
        my_dict = dict(
            SR = data[key]['selection_vec_dict']['f1'][0],
            RF = float(key.split("rf=")[1]),
            dose = float(key.split("dose=")[1].split(",")[0]),
            yield_ = data[key]['yield_vec'][0],
            )
        rows_list.append(my_dict)
    
    df = pd.DataFrame(rows_list)
    

    traces = []

    all_RFs = list(df.RF.unique())
    all_RFs.reverse()

    for ind, rr in enumerate(all_RFs):
        filt_df = df[(df['RF']==rr) & (df['yield_']>95)]

        min_color = 0
        max_color = 255

        clr = min_color + (max_color - min_color)* ind/(len(all_RFs)-1)

        color = f"rgb({255-clr},{0},{clr})"

        x = np.asarray(filt_df.dose)
        y = np.asarray(filt_df.SR)
        
        line = go.Scatter(x=x, y=y, 
                name=f"Resistance frequency: {str(rr)}",
                line=dict(color=color))

        traces.append(line)
    
    annotz = []
    
    for text, show_arr, y_pos in zip(['Increasing<br>resistance<br>frequency', ''],
                                    [False, True],
                                    [0.94, 0.82],
                                    ):
        annotz.append(dict(x=1,
                y=y_pos,
                text=text,
                showarrow=show_arr,
                
                xref='paper',
                yref='paper',
                arrowcolor=LABEL_COLOR,
                arrowsize=2,
                arrowwidth=1,
                arrowhead=2,
                # arrowside="start",
                
                ax=0,
                ay=380,
                
                xanchor="center",
                yanchor="top",

                font=dict(
                        size=14,
                        color=LABEL_COLOR,
                    ),
                ))

    fig = go.Figure(data=traces, layout=standard_layout(True))
                    
    fig.update_layout(legend=dict(x=0.05, y=1, 
                        bgcolor='rgba(255,255,255,0.5)'),
                    annotations=annotz)

    fig.update_xaxes(title="Dose")
    fig.update_yaxes(title="Selection ratio")

    return fig

# End of Changing dose

#----------------------------------------------------------------------------------------------
# * Grid of tactics

def dose_grid_heatmap(data, Config, to_plot):
    traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])

    trace = go.Heatmap(
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


#----------------------------------------------------------------------------------------------
# * Dose space

def dose_space_contour(data, to_plot):
    traces = []
    
    x = [log2(xx) for xx in data['contours_radial']]
    # y = [log2(xx) for xx in data['contour_perp']]
    y = data['contour_perp']

    z = data[to_plot]
    # np.transpose?

    trace = go.Heatmap(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = TITLE_MAP[to_plot],
            titleside = 'right',
        )
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=660, height=600)

    spaces = max(floor(len(x)/5),1)
    xticks = x[0:len(x):spaces]
    xtick_text = list(map(str, data['contours_radial']))[0:len(x):spaces]

    fig.update_xaxes(title="Ratio (log scale)",
        tickvals=xticks,
        ticktext=xtick_text)

    
    # ytick_text = [str(round(xx,3)) for xx in data['contour_perp']][0:len(y):spaces]
    # yticks = y[0:len(y):spaces]

    fig.update_yaxes(title="'Strength'", # (log scale)",
        showgrid=False,
        # tickvals=yticks,
        # ticktext=ytick_text
        )

    return fig


# End of Dose space
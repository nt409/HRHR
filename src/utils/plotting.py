import plotly.graph_objects as go
from plotly.subplots import make_subplots
import itertools
import numpy as np
import pandas as pd
from math import log2, floor, log10, pi
import unicodedata

from .params import PARAMS

# TOC
# Utility functions
# Single Tactic
# Changing dose
# Grid of tactics
# Dose space

pi_str = unicodedata.lookup("GREEK SMALL LETTER PI")

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

def yield_res_freqs_plot(data):
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    
    rf_traces = []

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
            name = f"Resistance to {titles[key].lower()}"
        )

        rf_traces.append(line)
        fig.add_trace(line, row=2, col=1)

    fig.update_xaxes(title="Year", row=2, col=1, showgrid=False)
    fig.update_yaxes(title="Resistant<br>frequency", row=2, col=1)


    y = data['yield_vec']
    x = list(range(1,1+len(y)))

    line = go.Scatter(
        x = x,
        y = y,
        line=dict(color="green"),
        name="Yield",
        showlegend=False,
    )

    y_low = y[-1]-2
    
    annotz = [dict(
        xref="x1",
        yref="y1",
        x=5,
        y=0.5*(95+y_low),
        text="Unacceptable yield",
        showarrow=False,
        )]

    shape = go.Scatter(x=[0, 0, x[-1], x[-1]],
                        y=[y_low, 95, 95, y_low],
                        fill="toself",
                        mode="lines",
                        showlegend=False,
                        line=dict(width=0, color="rgb(150,150,150)"))
    
    fig.add_trace(shape, col=1, row=1)
    
    fig.add_trace(line, col=1, row=1)

    fig.update_yaxes(title="Yield<br>(% of disease free)", row=1, col=1)

    fig.update_layout(standard_layout(True))
    fig.update_layout(annotations=annotz)
    fig.update_layout(legend=dict(x=0.1,
                            y=0.35,
                            # yref="paper",
                            bgcolor="rgba(255,255,255,0.5)"))

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

    fig.update_xaxes(title="Dose (fungicide 1)")
    fig.update_yaxes(title="Dose (fungicide 2)")

    return fig





def dose_grid_heatmap_with_log_ratio(data, Config, to_plot):
    """
    Lines on heatmap, and then log ratio of RFs at break down
    """

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

    subplot1_traces = []
    subplot2_traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])

    heatmap = go.Heatmap(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = TITLE_MAP[to_plot], # title here
            titleside = 'right',
            # len=0.4,
            # y=0.58,
            # yanchor="bottom",
        )
    )

    subplot2_traces.append(heatmap)
    

    # add lines on plot
    up_to_2 = np.linspace(0,2,2*Config.n_doses-1)
    ind0 = 0
    mn = 0
    mx = 255

    colors = {}

    for ind in range(len(up_to_2)):

        dose_line_non_0 = False
        
        for i in range(ind):
            j = ind-i
            if i<z.shape[0] and j<z.shape[1] and z[i,j]>0:
                dose_line_non_0 = True

        if not dose_line_non_0:
            ind0 = ind
            continue

        clr = mn + (mx-mn)* (ind-ind0)/(len(up_to_2)-ind0-1)
        colors[ind] = f"rgba(0,{255-clr},{clr},0.8)"
        
        xx = up_to_2[:ind+1]
        xx = [x for x in xx if (x<=1 and x>=0)]
        yy = [up_to_2[ind] - xx[i] for i in range(len(xx))]

        if len(yy)==Config.n_doses:
            yy = [y for y in yy if (y<=1)]
            xx = xx[len(xx)-len(yy):]

        ds = round(up_to_2[ind], 2)
        subplot2_traces.append(
            go.Scatter(
                x = xx,
                y = yy,
                line=dict(color=colors[ind]),
                name=f"Dose sum: {ds}",
            )
            )

    # col 2    
    for trace in subplot2_traces:
        fig.add_trace(trace, row=1, col=2)
    

    FY = data["FY"]
    for ind in range(ind0,len(up_to_2)):
        x = []
        y = []

        for i in range(ind+1):
            
            j = ind-i

            if i<FY.shape[0] and j<FY.shape[1]:
                fy = int(FY[i,j])
                
                if fy>0:
                    s1 = data['res_arrays']['f1'][i,j,fy]
                    s2 = data['res_arrays']['f2'][i,j,fy]

                    x.append(log10(s1/s2))
                    y.append(fy)
                
        if not x:
            continue
        
        line = go.Scatter(x=x,
                        y=y,
                        showlegend=False,
                        # mode='markers',
                        line=dict(color=colors[ind]))
        subplot1_traces.append(line)
    
    for trace in subplot1_traces:
        fig.add_trace(trace, row=1, col=1)


    annotz = []

    for x_pos, y_pos, text, show_arr, arrow_length in zip(
            [0, 0.4, 0.34, 0.06],
            [-0.06, -0.06, -0.08, -0.08],
            ['More resistance<br>to f2', 'More resistance<br>to f1', '', ''],
            [False, False, True, True],
            [None, None, -200, 200]):
        annotz.append(dict(
            x=x_pos,
            y=y_pos,
            text=text,

            showarrow=show_arr,
            arrowcolor=LABEL_COLOR,
            arrowsize=2,
            arrowwidth=1,
            arrowhead=2,
            
            ax=arrow_length,
            ay=0,
                
            xref='paper',
            yref='paper',

            xanchor="center",
            yanchor="top",

            font=dict(
                    size=14,
                    color=LABEL_COLOR,
                ),
        ))


    fig.update_layout(standard_layout(True))
    fig.update_layout(width=1500, 
                        height=650,
                        annotations=annotz,
                        legend=dict(x=1.15,
                                    y=0.5,
                                    yanchor="middle",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20)
                                    )

    fig.update_xaxes(title="Log ratio of resistance<br>frequencies at breakdown", row=1, col=1)
    fig.update_yaxes(title="Failure year", row=1, col=1)
    
    fig.update_xaxes(title="Dose (fungicide 1)", row=1, col=2)
    fig.update_yaxes(title="Dose (fungicide 2)", row=1, col=2)

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


    fig.update_yaxes(title="'Strength'",
        showgrid=False,
        )

    return fig


def radial(radial_data, grid_data, Config):
    
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(grid_data['FY'])

    heatmap = go.Heatmap(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = TITLE_MAP['FY'],
            titleside = 'right',
        )
        )

    fig.add_trace(heatmap, row=1, col=2)

    angles = list(radial_data.angle.unique())
    
    for ind, angle in enumerate(angles):
        data_use = radial_data[radial_data["angle"]==angle]
        xx = list(data_use.d1)
        yy = list(data_use.d2)
        
        value = 255*ind/(len(angles)-1)

        clr = f"rgb(0,{value},{255-value})"

        heatmap_line = go.Scatter(x=xx,
                            y=yy,
                            line=dict(color=clr),
                            showlegend=False
                            )
        fig.add_trace(heatmap_line, row=1, col=2)


        x = list(data_use.radius)
        y = list(data_use.FY)
        
        ang = 2*90*angle/pi
        ang = round(ang, 0)

        name_string = f"Angle: {ang}" + u"\u00B0"
        
        line = go.Scatter(x=x,
                        y=y,
                        line=dict(color=clr),
                        name= name_string
                        )

        fig.add_trace(line, row=1, col=1)

    
    annotz = []

    for x_pos, y_pos, text, show_arr, arrow_length in zip(
            [0.4, 0.34],
            [-0.055, -0.08],
            ['Stronger fungicide<br>mixture', ''],
            [False, True],
            [None, -200]):
        annotz.append(dict(
            x=x_pos,
            y=y_pos,
            text=text,

            showarrow=show_arr,
            arrowcolor=LABEL_COLOR,
            arrowsize=2,
            arrowwidth=1,
            arrowhead=2,
            
            ax=arrow_length,
            ay=0,
                
            xref='paper',
            yref='paper',

            xanchor="center",
            yanchor="top",

            font=dict(
                    size=14,
                    color=LABEL_COLOR,
                ),
        ))

    fig.update_layout(standard_layout(True))
    fig.update_layout(width=1500,
                        height=650,
                        legend=dict(x=1.15,
                                    y=0.5,
                                    yanchor="middle",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20),
                        annotations=annotz
                        )


    fig.update_xaxes(title="Mixture strength<br>(radius measured from origin)", row=1, col=1)
    fig.update_yaxes(title="Failure year", row=1, col=1)

    fig.update_xaxes(title="Dose (fungicide 1)", row=1, col=2)
    fig.update_yaxes(title="Dose (fungicide 2)", row=1, col=2)

    return fig

# End of Dose space
from numpy.core.numerictypes import sctype2char
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from math import log2, floor, log10, pi, ceil

from .params import PARAMS

# TOC
# Utility functions
# Single Tactic
# Changing dose
# Changing fcide
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
NULL_HEATMAP_COLOUR = "rgb(100, 100, 100)"

strain_attrs = dict(SS=dict(color="rgb(34,140,34)", dash="dot", longname='Double sensitive'),
                RS=dict(color="rgb(20,20,200)", dash="dash", longname='Single resistant (RS)'),
                SR=dict(color="rgb(200,20,20)", dash="dash", longname='Single resistant (SR)'),
                RR=dict(color="rgb(50,50,50)", dash="solid", longname='Double resistant'),
                )

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

def yield_by_year(data, conf_str):
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

    fig.show()
    filename = conf_str.replace("/single/", "/single/yield_by_year/")
    fig.write_image(filename)




def res_freqs_single_t_plot(data, conf_str):
    traces = []

    titles = dict(
        f1 = "Fungicide A",
        f2 = "Fungicide B"
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

    fig.show()
    filename = conf_str.replace("/single/", "/single/res_freqs/")
    fig.write_image(filename)



def single_year_plot(data, indices, conf_str):
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

    fig.show()
    filename = conf_str.replace("/single/", "/single/within_season/")
    fig.write_image(filename)



def yield_res_freqs_plot(data, conf_str):
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    
    rf_traces = []

    titles = dict(
        f1 = "Fungicide A",
        f2 = "Fungicide B"
    )
    
    for key in ['f1', 'f2']:
        y = data['res_vec_dict'][key]
        x = list(range(len(y)))

        line = go.Scatter(
            x = x,
            y = y,
            mode="lines+markers",
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
        mode="lines+markers",
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

    fig.show()
    filename = conf_str.replace("/single/", "/single/yield_rf/")
    fig.write_image(filename)



def plot_frequencies(data, conf_str):
    names = list(data['end_of_season'].keys())
    names.remove("SS")
    traces = []

    year = 4

    legend_entries = [f"Start of season {year}",
                       f"End of season {year} (before sexual reproduction)", 
                       f"Start of season {year+1} (after SR step)"]
    
    for key, yr, legend_entry in zip(['start_of_season', 'end_of_season', 'start_of_season'],
                    [year, year, year+1],
                    legend_entries):
        
        y = []
        for ff in names:
            y.append(log10(data[key][ff][int(yr)]))

        bar = go.Bar(x=names,
                        y=y,
                        name=legend_entry
                    )
        
        traces.append(bar)

    fig = go.Figure(data=traces, layout=standard_layout(True))
    fig.update_layout(legend=dict(x=0.4,
                            y=0.1,
                            bgcolor="rgba(255,255,255,0.5)"))
    fig.update_layout(barmode='group')
    fig.update_xaxes(title="Pathogen strain")
    fig.update_yaxes(title="Frequency (log base 10)")
    
    fig.show()
    filename = conf_str.replace("/single/", "/single/strain_freqs/")
    fig.write_image(filename)



def plot_frequencies_over_time(data, conf_str):
    """
    Logit scale, all strain freqs vs year.

    Within season and between season points.
    """
    
    season_frac = 0.75
    
    traces = []
    shapes = []

    rf_s = data['start_of_season']
    rf_e = data['end_of_season']

    n_yr_run = len(rf_s['RR']) - 1
    
    keys = list(rf_s.keys())
    keys.reverse()

    # generate y upper/lower for shaded backgd
    min_ = 1
    max_ = 0
    for key in keys:
        min_ = min(np.amin(rf_s[key][:n_yr_run]), np.amin(rf_e[key][:n_yr_run]), min_)
        max_ = max(np.amax(rf_s[key][:n_yr_run]), np.amax(rf_e[key][:n_yr_run]), max_)
    
    shape_min = log10(min_/(1-min_)) - 0.2
    shape_max = log10(max_/(1-max_)) + 0.2
    

    for key in keys:
        y = []
        x = []

        x_scat = []
        y_scat = []
        
        for i in range(n_yr_run):
            x.append(i)
            y.append(log10(rf_s[key][i]/(1-rf_s[key][i])))

            y_scat.append(log10(rf_s[key][i]/(1-rf_s[key][i])))
            
            if i!=n_yr_run-1:
                x.append(i+season_frac)
                y.append(log10(rf_e[key][i]/(1-rf_e[key][i])))
                
            
                shapes.append(go.Scatter(x=[i+season_frac, i+1, i+1, i+season_frac],
                                y=[shape_min, shape_min, shape_max, shape_max],
                                fill="toself",
                                mode="lines",
                                showlegend=False,
                                line=dict(width=0, color="rgb(230,250,255)")))

        x_scat = list(range(len(y)))

        
        line = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        name=strain_attrs[key]['longname'],
                        line=dict(color=strain_attrs[key]['color'],
                                dash=strain_attrs[key]['dash'])
                        )
        
        scatter = go.Scatter(x=x_scat,
                        y=y_scat, 
                        name=strain_attrs[key]['longname'],
                        line=dict(color=strain_attrs[key]['color'],
                                dash=strain_attrs[key]['dash']),
                        showlegend=False,
                        mode="markers",
                        )
        
        traces.append(line)
        traces.append(scatter)

    traces = shapes + traces

    fig = go.Figure(data=traces, layout=standard_layout(True))
    fig.update_xaxes(title="Year")
    fig.update_yaxes(title="Frequency (logit scale)")
    
    fig.update_layout(legend=dict(x=0.07,
                        y=0.8,
                        font=dict(size=18),
                        bgcolor="rgba(255,255,255,0.5)"))

    fig.show()
    filename = conf_str.replace("/single/", "/single/strain_freqs/overtime")
    fig.write_image(filename)



# End of single Tactic

#----------------------------------------------------------------------------------------------
# * Changing dose

def SR_by_dose_plot(data, conf_str):
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
        
        line = go.Scatter(x=x, 
                y=y, 
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

    fig.show()
    filename = conf_str.replace("/single/", "/changing_dose/equal_ratio/")
    fig.write_image(filename)

# End of Changing dose

#----------------------------------------------------------------------------------------------
# * Changing fcide

def fcide_grid(x, y, z, filename, labels):
    traces = []

    trace = go.Heatmap(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = labels['cbar'],
            titleside = 'right',
        )
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=660, height=600)

    fig.update_xaxes(title=labels['x'])
    fig.update_yaxes(title=labels['y'])

    fig.show()
    fig.write_image(filename)

# End of Changing fcide

#----------------------------------------------------------------------------------------------
# * Grid of tactics

def dose_grid_heatmap(data, Config, to_plot, conf_str):
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

    fig.update_xaxes(title="Dose (fungicide A)")
    fig.update_yaxes(title="Dose (fungicide B)")

    fig.show()
    filename = conf_str.replace("/grid/", "/grid/dose_grid/")
    fig.write_image(filename)





def dose_sum_hobb_vs_me(data, Config, to_plot, conf_str):
    """
    Lines on heatmap, and then log ratio of RFs at break down
    """

    fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2)

    my_strat_traces = []
    heatmap_subplot = []
    hobb_strat_traces = []
    
    xheat = np.linspace(0, 1, Config.n_doses)
    yheat = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])

    for name_, clr, dash_ in zip(["Equal selection", "First year selection"],
                                ["blue", "black"],
                                ["dash", "dot"]):
        heatmap_subplot.append(go.Scatter(
            x=[-0.3,-0.31],
            y=[-0.3,-0.31],
            line=dict(color=clr, dash=dash_),
            mode="lines",
            name=name_
        ))

    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=[
            [0, NULL_HEATMAP_COLOUR],
            [1/np.amax(z), NULL_HEATMAP_COLOUR],
            [1/np.amax(z), "rgb(0, 0, 100)"],
            [1, "rgb(255, 255, 0)"],
        ],
        colorbar=dict(
            x=0.42, y=0.79, len=0.43,
            title = TITLE_MAP[to_plot], # title here
            titleside = 'right',
        )
    )

    heatmap_subplot.append(heatmap)
    

    # add lines on heatmap
    y_intrcpt = np.linspace(0,2,2*Config.n_doses-1)
    ind0 = 0
    mn = 0
    mx = 220

    colors = {}
    inds_list = []
    
    for ind in range(len(y_intrcpt)):

        dose_line_all_0 = True
        
        for i in range(ind):
            j = ind-i
            if i<z.shape[0] and j<z.shape[1] and z[i,j]>0:
                dose_line_all_0 = False

        if dose_line_all_0:
            ind0 = ind
            continue
        
        n_lines = 4
        n_interval = ceil((len(y_intrcpt) - ind0)/n_lines)
        # only actually want n_lines lines total
        if not (((ind - 1 -ind0) % n_interval == 0 and (len(y_intrcpt) - 1 - ind > n_interval))
                         or (ind==len(y_intrcpt)-1)):
            # only want a line every n_interval steps...
            # and want final point but then not another point v nearby
            continue

        clr = mn + (mx-mn)* (ind - 1 - ind0)/(len(y_intrcpt) - ind0 - 1)

        colors[ind] = f"rgba({255-clr},{0},{255-clr},0.8)"
        
        xx = y_intrcpt[:ind+1]
        xx = [x for x in xx if (x<=1 and x>=0)]
        yy = [y_intrcpt[ind] - x for x in xx]

        if len(yy)==Config.n_doses:
            # 1.001 since floating point error caused probs :/
            yy = [y for y in yy if (y<=1.001)]
            xx = xx[len(xx)-len(yy):]
        
        inds_list.append(ind)

        ds = round(y_intrcpt[ind], 2)
        heatmap_subplot.append(
            go.Scatter(
                x = xx,
                y = yy,
                line=dict(color=colors[ind]),
                name=f"Dose sum: {ds}",
            )
            )

    

    eq_sel = np.zeros(z.shape)
    eq_fy = np.zeros(z.shape)
    
    FY = data["FY"]
    for ind in range(len(y_intrcpt)):
        x_fy = []
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

                    eq_sel[j, i] = log10(s1/s2)

                    f1_y0 = data['start_freqs']['RS'][i,j,0]
                    f1_y1 = data['start_freqs']['RS'][i,j,1]
                    
                    f2_y0 = data['start_freqs']['SR'][i,j,0]
                    f2_y1 = data['start_freqs']['SR'][i,j,1]

                    s1_y1 = f1_y1/f1_y0
                    s2_y1 = f2_y1/f2_y0

                    # print(f1_y0, f1_y1, s1_y1)
                    # print(f2_y0, f2_y1, s2_y1)

                    x_fy.append(log10(s1_y1/s2_y1))
                    eq_fy[j, i] = log10(s1_y1/s2_y1)
                else:
                    eq_sel[j, i] = None
                    eq_fy[j, i] = None
                
        if not x or not (ind in inds_list):
            continue
        
        myline = go.Scatter(x=x,
                        y=y,
                        showlegend=False,
                        # name=f"Equal breakdown (dose sum): {round(y_intrcpt[ind],2)}",
                        line=dict(color=colors[ind], dash="solid"))
        my_strat_traces.append(myline)
        
        hobb_line_fy = go.Scatter(x=x_fy,
                        y=y,
                        showlegend=False,
                        line=dict(color=colors[ind], dash="dot"),
                        # name=f"First year selection (dose sum): {round(y_intrcpt[ind],2)}"
                        )
        hobb_strat_traces.append(hobb_line_fy)
    
    for trace in my_strat_traces:
        fig.add_trace(trace, row=1, col=2)
    
    for trace in hobb_strat_traces:
        fig.add_trace(trace, row=2, col=2)
    
    eq_contour = go.Contour(x=xheat,
                    y=yheat,
                    z=eq_sel,
                    contours=dict(start=0, end=0),
                    contours_coloring='lines',
                    colorscale=["blue", "blue"],
                    
                    line=dict(width=2, dash="dash"),

                    # hacky way to remove second colorbar
                    colorbar=dict(x=0.42, len=0.1, 
                            tickfont=dict(size=1,
                                color="rgba(0,0,0,0)"
                                )), 
                    )

    eq_fy_contour = go.Contour(x=xheat,
                    y=yheat,
                    z=eq_fy,
                    contours=dict(start=0, end=0),
                    contours_coloring='lines',
                    colorscale=["black", "black"],
                    
                    line=dict(width=2, dash="dot"),

                    # hacky way to remove second colorbar
                    colorbar=dict(x=0.42, len=0.1, 
                            tickfont=dict(size=1,
                                color="rgba(0,0,0,0)"
                                )), 
                    )

    heatmap_subplot.append(eq_contour)
    heatmap_subplot.append(eq_fy_contour)

    # col 1, row 1
    for trace in heatmap_subplot:
        fig.add_trace(trace, row=1, col=1)
    
    

    annotz = []

    for x_pos, y_pos, text, show_arr, arrow_length in zip(
            [0.58, 1.02, 0.96, 0.64],
            [-0.025, -0.025, -0.04, -0.04],
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
    fig.update_layout(width=1200, 
                        height=1150,
                        annotations=annotz,
                        legend=dict(
                                    x=0.25,
                                    # x=1.2,
                                    y=0.25,
                                    yanchor="middle",
                                    xanchor="center",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=18)
                                    )

    fig.update_xaxes(title="Log ratio of resistance<br>frequencies at breakdown", row=1, col=2, showgrid=False)
    fig.update_xaxes(title="Log ratio of selection<br>ratios after one year", row=2, col=2, showgrid=False)
    fig.update_yaxes(title="Failure year", row=1, col=2)
    fig.update_yaxes(title="Failure year", row=2, col=2)
    
    # if heatmap not contour use [0-dx, 1+dx] etc
    # is order correct? shape[0]/[1]
    # dx = 0.5*(1/(-1+z.shape[1]))
    # dy = 0.5*(1/(-1+z.shape[0]))
    dx = 0.01
    dy = 0.01
    
    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=1, showgrid=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=1, showgrid=False)

    fig.show()
    filename = conf_str.replace("/grid/", "/dose_space/dose_sum_hobb_vs_me/")
    fig.write_image(filename)








def dose_sum_LR(data, Config, to_plot, conf_str):
    """
    Lines on heatmap, and then log ratio of RFs at break down
    """

    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

    my_strat_traces = []
    heatmap_subplot = []
    
    xheat = np.linspace(0, 1, Config.n_doses)
    yheat = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])

    for name_, clr, dash_ in zip(["Equal selection", "First year selection"],
                                ["blue", "black"],
                                ["dash", "dot"]):
        heatmap_subplot.append(go.Scatter(
            x=[-0.3,-0.31],
            y=[-0.3,-0.31],
            line=dict(color=clr, dash=dash_),
            mode="lines",
            name=name_
        ))

    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=[
            [0, NULL_HEATMAP_COLOUR],
            [1/np.amax(z), NULL_HEATMAP_COLOUR],
            [1/np.amax(z), "rgb(0, 0, 100)"],
            [1, "rgb(255, 255, 0)"],
        ],
        colorbar=dict(
            title = TITLE_MAP[to_plot], # title here
            titleside = 'right',
        )
    )

    heatmap_subplot.append(heatmap)
    

    # add lines on heatmap
    y_intrcpt = np.linspace(0,2,2*Config.n_doses-1)
    ind0 = 0
    mn = 0
    mx = 220

    colors = {}
    inds_list = []
    
    for ind in range(len(y_intrcpt)):

        dose_line_all_0 = True
        
        for i in range(ind):
            j = ind-i
            if i<z.shape[0] and j<z.shape[1] and z[i,j]>0:
                dose_line_all_0 = False

        if dose_line_all_0:
            ind0 = ind
            continue
        
        n_lines = 4
        n_interval = ceil((len(y_intrcpt) - ind0)/n_lines)
        # only actually want n_lines lines total
        if not (((ind - 1 -ind0) % n_interval == 0 and (len(y_intrcpt) - 1 - ind > n_interval))
                         or (ind==len(y_intrcpt)-1)):
            # only want a line every n_interval steps...
            # and want final point but then not another point v nearby
            continue

        clr = mn + (mx-mn)* (ind - 1 - ind0)/(len(y_intrcpt) - ind0 - 1)

        colors[ind] = f"rgba({255-clr},{0},{255-clr},0.8)"
        
        xx = y_intrcpt[:ind+1]
        xx = [x for x in xx if (x<=1 and x>=0)]
        yy = [y_intrcpt[ind] - x for x in xx]

        if len(yy)==Config.n_doses:
            # 1.001 since floating point error caused probs :/
            yy = [y for y in yy if (y<=1.001)]
            xx = xx[len(xx)-len(yy):]
        
        inds_list.append(ind)

        ds = round(y_intrcpt[ind], 2)
        heatmap_subplot.append(
            go.Scatter(
                x = xx,
                y = yy,
                line=dict(color=colors[ind]),
                name=f"Dose sum: {ds}",
            )
            )

    

    eq_sel = np.zeros(z.shape)
    eq_fy = np.zeros(z.shape)
    
    FY = data["FY"]
    for ind in range(len(y_intrcpt)):
        x_fy = []
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

                    eq_sel[j, i] = log10(s1/s2)

                else:
                    eq_sel[j, i] = None
                    eq_fy[j, i] = None
                
        if not x or not (ind in inds_list):
            continue
        
        myline = go.Scatter(x=x,
                        y=y,
                        showlegend=False,
                        # name=f"Equal breakdown (dose sum): {round(y_intrcpt[ind],2)}",
                        line=dict(color=colors[ind], dash="solid"))
        my_strat_traces.append(myline)
    
    for trace in my_strat_traces:
        fig.add_trace(trace, row=1, col=1)
    
    eq_contour = go.Contour(x=xheat,
                    y=yheat,
                    z=eq_sel,
                    contours=dict(start=0, end=0),
                    contours_coloring='lines',
                    colorscale=["blue", "blue"],
                    
                    line=dict(width=2, dash="dash"),

                    # hacky way to remove second colorbar
                    colorbar=dict(x=0.42, len=0.1, 
                            tickfont=dict(size=1,
                                color="rgba(0,0,0,0)"
                                )), 
                    )

    heatmap_subplot.append(eq_contour)

    for trace in heatmap_subplot:
        fig.add_trace(trace, row=1, col=2)
    
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
    fig.update_layout(width=1360, 
                        height=650,
                        annotations=annotz,
                        legend=dict(
                                    x=0.36,
                                    y=0.95,
                                    yanchor="top",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20)
                                    )

    fig.update_xaxes(title="Log ratio of resistance<br>frequencies at breakdown", row=1, col=1, showgrid=False)
    fig.update_yaxes(title="Failure year", row=1, col=1)
    
    # if heatmap not contour use [0-dx, 1+dx] etc
    # is order correct? shape[0]/[1]
    # dx = 0.5*(1/(-1+z.shape[1]))
    # dy = 0.5*(1/(-1+z.shape[0]))
    dx = 0.01
    dy = 0.01
    
    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=2, showgrid=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=2, showgrid=False)

    fig.show()
    filename = conf_str.replace("/grid/", "/dose_space/dose_sum/")
    fig.write_image(filename)


# End of grid of tactics


#----------------------------------------------------------------------------------------------
# * Dose space

def dose_space_contour(data, to_plot, conf_str):
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

    fig.show()
    filename = conf_str.replace("/grid/", "/dose_space/yield_by_year/")
    fig.write_image(filename)


def radial(radial_data, grid_data, Config):
    
    fig = make_subplots(rows=1, cols=2, horizontal_spacing=0.2)

    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(grid_data['FY'])

    heatmap = go.Contour(
        x = x,
        y = y,
        z = z,
        colorscale=[
            [0, NULL_HEATMAP_COLOUR],
            [1/np.amax(z), NULL_HEATMAP_COLOUR],
            [1/np.amax(z), "rgb(0, 0, 100)"],
            [1, "rgb(255, 255, 0)"],
        ],
        colorbar=dict(
            title = TITLE_MAP['FY'],
            titleside = 'right',
        )
        )

    fig.add_trace(heatmap, row=1, col=2)

    angles = list(radial_data.angle.unique())
    
    for ind, angle in enumerate(angles):
        data_use = radial_data[(radial_data["angle"]==angle) & (radial_data['FY']>0)]
        
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
    fig.update_layout(width=1360,
                        height=650,
                        legend=dict(x=0.35,
                                    y=0.2,
                                    yanchor="middle",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20),
                        annotations=annotz
                        )


    fig.update_xaxes(title="Mixture strength<br>(radius measured from origin)", row=1, col=1)
    fig.update_yaxes(title="Failure year", row=1, col=1)

    dx = 0.01
    dy = 0.01

    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=2, showgrid=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=2, showgrid=False)

    fig.show()
    conf_str = Config.config_string_img
    filename = conf_str.replace("/grid/", "/dose_space/radial/")
    fig.write_image(filename)

# End of Dose space
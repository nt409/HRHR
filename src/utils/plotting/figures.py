
import plotly.colors as pltly_clrs
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import pandas as pd
from math import log2, floor, log10, pi
from scipy import stats
from PIL import Image



from ..params import PARAMS
from ..functions import EqualResFreqBreakdownArray, EqualSelectionArray

from .traces import get_RFB_diff_traces, get_eq_sel_traces, get_heatmap_lines, \
    get_strain_freq_traces, contour_at_0, get_multi_contour_traces, \
    get_MS_RFB_traces, contour_at_single_level

from .utils import get_text_annotation, get_arrow_annotation, standard_layout, \
    grey_colorscale, my_colorbar, get_big_text_annotation

from .consts import ATTRS_DICT, TITLE_MAP, PLOT_WIDTH, PLOT_HEIGHT, \
        FULL_PAGE_WIDTH



# TOC
# Single Tactic
# Changing dose
# Changing fcide
# Grid of tactics
# Dose space
# RF Ratio
# Paper Figs



#----------------------------------------------------------------------------------------------
# * Single Tactic

def yield_by_year(data, conf_str):
    traces = []
    
    y = data['yield_vec']
    x = list(range(1,1+len(y)))

    line = go.Scatter(
        x = x,
        y = y,
        line = dict(color="green"),
    )

    y_low = y[-1]-2
    
    annotz = [dict(
        xref="x1",
        yref="y1",
        x=10,
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
    
    traces.append(line)
    traces.append(shape)

    fig = go.Figure(data=traces, layout=standard_layout(False))
    fig.update_layout(annotations=annotz)

    fig.update_xaxes(title="Year")
    fig.update_yaxes(title="Yield<br>(% of disease free)")

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
            name = ATTRS_DICT[str(ind)]['name'],
            line = dict(color=ATTRS_DICT[str(ind)]['colour'], dash=ATTRS_DICT[str(ind)]['dash'])
        )

        traces.append(line)

    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_xaxes(title="Time (degree-days)")
    fig.update_yaxes(title="Amount")
    
    if 15 in indices:
        fig.update_layout(legend=dict(x=0.5, y=0.85))
    elif 0 in indices:
        fig.update_layout(legend=dict(x=0.2, y=0.55))
    else:
        fig.update_layout(legend=dict(x=0.01, y=0.95))
    
    fig.update_layout(width=PLOT_WIDTH*(2/3), height=PLOT_HEIGHT)

    fig.show()
    filename = conf_str.replace("/single/", "/single/within_season/plot" + "".join(str(e) for e in indices))
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

    fig.update_xaxes(title="Year", row=2, col=1, showgrid=False, zeroline=False)
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
    
    rf_s = data['start_of_season']
    rf_e = data['end_of_season']


    traces = get_strain_freq_traces(rf_s, rf_e)

    fig = go.Figure(data=traces, layout=standard_layout(True))
    fig.update_xaxes(title="Year")
    fig.update_yaxes(title="Frequency (logit scale)")
    
    fig.update_layout(legend=dict(x=0.07,
                        y=1,
                        orientation="h",
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
    
    

    text = get_text_annotation(1.05, 1.1, 'Increasing<br>resistance<br>frequency')
    arrow = get_arrow_annotation(1.05, 0.78, 0, 150)
    
    annotz = [text, arrow]

    fig = go.Figure(data=traces, layout=standard_layout(True))
                    
    fig.update_layout(legend=dict(x=0.02, y=1.15, 
                        bgcolor='rgba(255,255,255,0.5)',
                        font=dict(size=14)),
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
        colorbar=my_colorbar(labels['cbar'])
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH - 50)

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
        colorbar=my_colorbar(TITLE_MAP[to_plot])
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH - 50)

    fig.update_xaxes(title="Dose (fungicide A)")
    fig.update_yaxes(title="Dose (fungicide B)")

    fig.show()
    filename = conf_str.replace("/grid/", f"/grid/dose_grid/{to_plot}")
    fig.write_image(filename)



def dose_grid_RA_heatmap(data, Config, conf_str, yr):
    traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data["res_arrays"]["f2"][:,:,yr])

    trace = go.Contour(
        x = x,
        y = y,
        z = z,
        colorbar=my_colorbar("Res. f. B"),
        # colorscale= [(0,"yellow"), (1,"red")],
        # range_color=[0,1]
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH - 50)

    fig.update_xaxes(title="Dose (fungicide A)")
    fig.update_yaxes(title="Dose (fungicide B)")

    fig.show()
    filename = conf_str.replace("/grid/", f"/grid/dose_grid/res_array_{yr}")
    fig.write_image(filename)





def dose_sum_hobb_vs_me(data, Config, to_plot, conf_str):
    """
    Lines on heatmap, and then log ratio of RFs at break down
    """

    fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2)

    heatmap_subplot = []
    
    xheat = np.linspace(0, 1, Config.n_doses)
    yheat = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data[to_plot])
    
    # for legend
    for name_, clr, dash_ in zip(["Equal resistance<br>at breakdown", 
                                        "First year selection"],
                                ["blue", "black"],
                                ["dash", "dot"]):
        heatmap_subplot.append(go.Scatter(
            x=[-0.3,-0.31],
            y=[-0.3,-0.31],
            line=dict(color=clr, dash=dash_),
            mode="lines",
            name=name_
        ))
    
    clrbar = my_colorbar(TITLE_MAP[to_plot])
    clrbar.update(dict(x=0.42, y=0.79, len=0.43))

    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=grey_colorscale(z),
        colorbar=clrbar
    )

    heatmap_subplot.append(heatmap)
    

    # add lines on heatmap
    N_y_int = 2*Config.n_doses-1
    y_intrcpt = np.linspace(0,2,N_y_int)
    trc_out, colors, inds_list = get_heatmap_lines(Config, z, y_intrcpt)

    heatmap_subplot += trc_out

    my_strat_traces, RBF_diff = get_RFB_diff_traces(data, z, N_y_int, inds_list, colors)
    hobb_strat_traces, eq_fy = get_eq_sel_traces(data, z, N_y_int, inds_list, colors)
    
    for trace in my_strat_traces:
        fig.add_trace(trace, row=1, col=2)
    
    for trace in hobb_strat_traces:
        fig.add_trace(trace, row=2, col=2)
    
    eq_contour = contour_at_0(xheat, yheat, RBF_diff, "blue", "dash")
    eq_fy_contour = contour_at_0(xheat, yheat, eq_fy, "black", "dot")


    heatmap_subplot.append(eq_contour)
    heatmap_subplot.append(eq_fy_contour)

    # col 1, row 1
    for trace in heatmap_subplot:
        fig.add_trace(trace, row=1, col=1)
    
    

    annotz = []

    annot1 = get_text_annotation(0.58, -0.025, 'More selection<br>for f. B')
    annot2 = get_text_annotation(1.02, -0.025, 'More selection<br>for f. A')
    arrow1 = get_arrow_annotation(0.96, -0.04, -200, 0)
    arrow2 = get_arrow_annotation(0.64, -0.04,  200, 0)
    
    annot1b = get_text_annotation(0.58,  0.55, 'More resistance<br>to f. B')
    annot2b = get_text_annotation(1.02,  0.55, 'More resistance<br>to f. A')
    arrow1b = get_arrow_annotation(0.96, 0.535, -200, 0)
    arrow2b = get_arrow_annotation(0.64, 0.535,  200, 0)
    
    annotz += [annot1, annot2, arrow1, arrow2, annot1b, annot2b, arrow1b, arrow2b]


    fig.update_layout(standard_layout(True))
    fig.update_layout(width=2*PLOT_WIDTH, 
                        height=2*PLOT_WIDTH - 50,
                        annotations=annotz,
                        legend=dict(
                                    x=0.25,
                                    y=0.25,
                                    yanchor="middle",
                                    xanchor="center",
                                    font=dict(size=18)
                                    ),
                        font=dict(size=18)
                                    )

    fig.update_xaxes(title="Difference in logit of<br>res. freqs. at breakdown", row=1, col=2, showgrid=False, zeroline=False)
    fig.update_xaxes(title="Difference of log of<br>selection ratios after one year", row=2, col=2, showgrid=False, zeroline=False)
    fig.update_yaxes(title="Effective life", row=1, col=2)
    fig.update_yaxes(title="Effective life", row=2, col=2)
    
    # if heatmap not contour use [0-dx, 1+dx] etc
    # is order correct? shape[0]/[1]
    # dx = 0.5*(1/(-1+z.shape[1]))
    # dy = 0.5*(1/(-1+z.shape[0]))
    dx = 0.01
    dy = 0.01
    
    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=1, showgrid=False, zeroline=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=1, showgrid=False, zeroline=False)

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
    
    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=grey_colorscale(z),
        colorbar=my_colorbar(TITLE_MAP[to_plot])
    )

    heatmap_subplot.append(heatmap)

    for name_, clr, dash_ in zip(["Equal resistance<br>at breakdown"],
                                ["blue"],
                                ["dash"]):
        heatmap_subplot.append(go.Scatter(
            x=[-0.3,-0.31],
            y=[-0.3,-0.31],
            line=dict(color=clr, dash=dash_),
            mode="lines",
            name=name_
        ))


    

    # add lines on heatmap
    N_y_int = 2*Config.n_doses-1
    y_intrcpt = np.linspace(0,2,N_y_int)
    trc_out, colors, inds_list = get_heatmap_lines(Config, z, y_intrcpt)

    heatmap_subplot += trc_out

    
    my_strat_traces, RBF_diff = get_RFB_diff_traces(data, z, N_y_int, inds_list, colors)
    
    
    for trace in my_strat_traces:
        fig.add_trace(trace, row=1, col=1)
    
    eq_contour = contour_at_0(xheat, yheat, RBF_diff, "blue", "dash")

    heatmap_subplot.append(eq_contour)

    for trace in heatmap_subplot:
        fig.add_trace(trace, row=1, col=2)
    
    annotz = []

    # for x_pos, y_pos, text, show_arr, arrow_length in zip(
    #         [0, 0.4, 0.34, 0.06],
    #         [-0.06, -0.06, -0.08, -0.08],
    #         ['More resistance<br>to f. B', 'More resistance<br>to f. A', '', ''],
    #         [False, False, True, True],
    #         [None, None, -200, 200]):
    #     annotz.append(dict(
    #         x=x_pos,
    #         y=y_pos,
    #         text=text,

    #         showarrow=show_arr,
    #         arrowcolor=LABEL_COLOR,
    #         arrowsize=2,
    #         arrowwidth=1,
    #         arrowhead=2,
            
    #         ax=arrow_length,
    #         ay=0,
                
    #         xref='paper',
    #         yref='paper',

    #         xanchor="center",
    #         yanchor="top",

    #         font=dict(
    #                 size=14,
    #                 color=LABEL_COLOR,
    #             ),
    #     ))

    annot1 = get_text_annotation(0, -0.06, 'More resistance<br>to f. B')
    annot2 = get_text_annotation(0.4, -0.06, 'More resistance<br>to f. A')
    arrow1 = get_arrow_annotation(0.34, -0.08, -200, 0)
    arrow2 = get_arrow_annotation(0.06, -0.08,  200, 0)
    
    annotz += [annot1, annot2, arrow1, arrow2]


    fig.update_layout(standard_layout(True))
    fig.update_layout(width=2*PLOT_WIDTH, 
                        height=PLOT_WIDTH - 50,
                        annotations=annotz,
                        legend=dict(
                                    x=0.36,
                                    y=0.95,
                                    yanchor="top",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20)
                                    )

    fig.update_xaxes(title="Difference in logit of<br>resistance frequencies<br>at breakdown", row=1, col=1, showgrid=False, zeroline=False)
    fig.update_yaxes(title="Effective life", row=1, col=1)
    
    # if heatmap not contour use [0-dx, 1+dx] etc
    # is order correct? shape[0]/[1]
    # dx = 0.5*(1/(-1+z.shape[1]))
    # dy = 0.5*(1/(-1+z.shape[0]))
    dx = 0.01
    dy = 0.01
    
    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=2, showgrid=False, zeroline=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=2, showgrid=False, zeroline=False)

    fig.show()
    filename = conf_str.replace("/grid/", f"/dose_space/dose_sum/{to_plot}")
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
        colorbar=my_colorbar(TITLE_MAP[to_plot])
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH-50)

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
        colorscale=grey_colorscale(z),
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

    # for x_pos, y_pos, text, show_arr, arrow_length in zip(
    #         [0.4, 0.34],
    #         [-0.055, -0.08],
    #         ['Stronger fungicide<br>mixture', ''],
    #         [False, True],
    #         [None, -200]):
    #     annotz.append(dict(
    #         x=x_pos,
    #         y=y_pos,
    #         text=text,

    #         showarrow=show_arr,
    #         arrowcolor=LABEL_COLOR,
    #         arrowsize=2,
    #         arrowwidth=1,
    #         arrowhead=2,
            
    #         ax=arrow_length,
    #         ay=0,
                
    #         xref='paper',
    #         yref='paper',

    #         xanchor="center",
    #         yanchor="top",

    #         font=dict(
    #                 size=14,
    #                 color=LABEL_COLOR,
    #             ),
    #     ))
    
    annot1 = get_text_annotation(0.4, -0.055, 'Stronger fungicide<br>mixture')
    arrow1 = get_arrow_annotation(0.34, -0.08, -200, 0)
    
    annotz += [annot1, arrow1]


    fig.update_layout(standard_layout(True))
    fig.update_layout(width=2*PLOT_WIDTH,
                        height=2*PLOT_WIDTH - 50,
                        legend=dict(x=0.35,
                                    y=0.2,
                                    yanchor="middle",
                                    font=dict(size=14)
                                    ),
                        font=dict(size=20),
                        annotations=annotz
                        )


    fig.update_xaxes(title="Mixture strength<br>(radius measured from origin)", row=1, col=1)
    fig.update_yaxes(title="Effective life", row=1, col=1)

    dx = 0.01
    dy = 0.01

    fig.update_xaxes(title="Dose (fungicide A)", range=[0-dx,1+dx], row=1, col=2, showgrid=False)
    fig.update_yaxes(title="Dose (fungicide B)", range=[0-dy,1+dy], row=1, col=2, showgrid=False)

    fig.show()
    conf_str = Config.config_string_img
    filename = conf_str.replace("/grid/", "/dose_space/radial/")
    fig.write_image(filename)




def first_year_yield(data, Config):
    
    traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data['yield_array'][:,:,0])

    trace = go.Contour(
        x = x,
        y = y,
        z = z,
        colorbar=dict(
            title = 'Yield',
            titleside = 'right',
        )
    )

    traces.append(trace)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH-50)

    fig.update_xaxes(title="Dose (fungicide A)",
        showgrid=False)

    fig.update_yaxes(title="Dose (fungicide B)",
        showgrid=False,
        )

    fig.show()
    conf_str = Config.config_string_img
    filename = conf_str.replace("/grid/", "/dose_space/first_year_yield/")
    fig.write_image(filename)


def eq_RFB_contours(data, Config, title=None):
    
    traces = get_multi_contour_traces(data, Config)

    fig = go.Figure(data=traces, layout=standard_layout(False))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH-50)
    # dy = 0.5*(1/(-1+z.shape[0]))
    dx = 0.01
    dy = 0.01

    fig.update_xaxes(title="Dose (fungicide A)",
        range=[0-dx,1+dx],
        showgrid=False)

    fig.update_yaxes(title="Dose (fungicide B)",
        range=[0-dy,1+dy],
        showgrid=False)

    fig.update_layout(title=title)

    fig.show()
    conf_str = Config.config_string_img
    filename = conf_str.replace("/grid/", "/dose_space/eq_RFB_contours/")
    fig.write_image(filename)



def MS_RFB_scatter_plot(data, Config):

    traces = get_MS_RFB_traces(data)

    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_layout(width=PLOT_WIDTH, height=PLOT_WIDTH-50)

    fig.update_layout(legend=dict(x=1.05,
                            y=1.1,
                            yanchor="top",
                            xanchor="right",
                            font=dict(size=16),
                            # bgcolor="rgba(255,255,255,0.5)"
                            ))
    

    fig.update_xaxes(title=u"\u0394" + "<i><sub>RFB</sub></i>")
    fig.update_yaxes(title="Dose sum", range=[0.3, 2.04])

    fig.show()
    conf_str = Config.config_string_img
    filename = conf_str.replace("/grid/", "/dose_space/MS_RFB_scatter_plot/")
    fig.write_image(filename)




# End of Dose space


# * RF Ratio


def dose_grid_heatmap_with_contours(data, Config, contours, conf_str):
    traces = []
    
    x = np.linspace(0, 1, Config.n_doses)
    y = np.linspace(0, 1, Config.n_doses)

    z = np.transpose(data["FY"])

    trace = go.Heatmap(
        x = x,
        y = y,
        z = z,
        colorbar=my_colorbar(TITLE_MAP["FY"])
    )

    traces.append(trace)
    
    names = ["Equal R.Fs<br>at breakdown", "Equal selection<br>in first year"]
    for cont_df, name in zip(contours, names):
        xx = cont_df.f1
        yy = cont_df.f2

        scatter = go.Scatter(x=xx, y=yy, name=name)

        traces.append(scatter)

    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_layout(width=PLOT_WIDTH, 
                    height=PLOT_WIDTH-50,
                    legend=dict(x=1.25,
                                y=1,
                                yanchor="top",
                                xanchor="left",
                                font=dict(size=16)
                                ))

    fig.update_xaxes(title="Dose (fungicide A)")
    fig.update_yaxes(title="Dose (fungicide B)")

    fig.show()
    filename = conf_str.replace("/grid/", "/rf_ratio/dose_grid_and_contours/")
    fig.write_image(filename)


def outcomes_by_ratio(data, conf_str):
    traces = []
    
    x = [log10(d) for d in data.ratio]
    yEqS = data.EqS
    yRFB = data.RFB

    y_list = [yRFB, yEqS]
    names = ["STRATEGY: Equal resistance frequencies at breakdown",
            "STRATEGY: Equal selection in first year"
            ]

    for yy, name in zip(y_list, names):
        line = go.Scatter(
            x = x,
            y = yy,
            name=name
            )

        traces.append(line)
    
    
    fig = go.Figure(data=traces, layout=standard_layout(True))

    fig.update_layout(width=2*PLOT_WIDTH,
                    height=PLOT_WIDTH,
                    legend=dict(x=0.3,
                                y=1,
                                yanchor="top",
                                xanchor="left",
                                font=dict(size=16)
                                ))
    xticks = list(x[::2])
    
    if 0 not in xticks:
        xticks.append(0)


    xtick_text = [10**(x) for x in xticks]

    fig.update_xaxes(title="Ratio of initial resistance frequencies",
            tickvals=xticks,
            ticktext=xtick_text
            )
    fig.update_yaxes(title="Effective lifetime")

    fig.show()
    
    x_str = [str(int(i)) for i in x]
    ratio_string = "".join(x_str)
    filename = conf_str.replace("/grid/", f"/rf_ratio/outcomes/{ratio_string}")
    fig.write_image(filename)

# End of RF Ratio



# Paper Figs
    
class DiseaseProgressCurvesAll:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 650

        self.xx = data['t_vec']
        
        self.array = data['sol_array']

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)




    def _generate_figure(self):

        self.config = self.get_config()

        traces_dict = self.get_model_output_overview_traces()

        ugly_fig = self.add_traces_to_layout_model_output_overview(traces_dict)

        fig = self._sort_layout(ugly_fig)

        return fig
        
    


    @staticmethod
    def get_config():

        E_cols = pltly_clrs.n_colors("rgb(0,0,255)", "rgb(255,255,255)", 5, colortype="rgb")[:4]
        
        I_cols = pltly_clrs.n_colors("rgb(255,0,0)", "rgb(255,255,255)", 5, colortype="rgb")[:4]


        return {
            'S': dict(color='limegreen', dash='solid', name='Susceptible', ind=PARAMS.S_ind),
            'R': dict(color='rgb(100,100,100)', dash='solid', name='Removed', ind=PARAMS.R_ind),

            'ERR': dict(color=E_cols[0], dash='dot', name='E (rr)', ind=PARAMS.ER_ind),
            'ERS': dict(color=E_cols[1], dash='dash', name='E (rs)', ind=PARAMS.ERS_ind),
            'ESR': dict(color=E_cols[2], dash='dashdot', name='E (sr)', ind=PARAMS.ESR_ind),
            'ESS': dict(color=E_cols[3], dash='solid', name='E (ss)', ind=PARAMS.ES_ind),

            'IRR': dict(color=I_cols[0], dash='dot', name='I (rr)', ind=PARAMS.IR_ind),
            'IRS': dict(color=I_cols[1], dash='dash', name='I (rs)', ind=PARAMS.IRS_ind),
            'ISR': dict(color=I_cols[2], dash='dashdot', name='I (sr)', ind=PARAMS.ISR_ind),
            'ISS': dict(color=I_cols[3], dash='solid', name='I (ss)', ind=PARAMS.IS_ind),

            'F1': dict(color='turquoise', dash='solid', name='Fungicide A', ind=PARAMS.Fung1_ind),
            'F2': dict(color='magenta', dash='dot', name='Fungicide B', ind=PARAMS.Fung2_ind),
        }


    def get_model_output_overview_traces(self):
        S_R_traces = self.get_S_R_traces()
        E_traces = self.get_E_traces()
        I_traces = self.get_I_traces()
        F_traces = self.get_F_traces()

        out = dict(S_R = S_R_traces,
                    E = E_traces,
                    I = I_traces,
                    F = F_traces)

        return out
    



    def get_S_R_traces(self):
        
        out = []

        for key in ['S', 'R']:
            out.append(self.get_DPC_trace(key))
        return out
    
    

    def get_E_traces(self):
        
        out = []

        for key in ['ERR', 'ERS', 'ESR', 'ESS']:
            out.append(self.get_DPC_trace(key))

        return out
    


    
    def get_I_traces(self):
        
        out = []

        for key in ['IRR', 'IRS', 'ISR', 'ISS']:
            out.append(self.get_DPC_trace(key))

        return out
    


    def get_F_traces(self):

        out = []

        for key in ['F1', 'F2']:
            out.append(self.get_DPC_trace(key))

        return out



    def get_DPC_trace(self, key):
        clr = self.config[key]['color']
        dash = self.config[key]['dash']
        name = self.config[key]['name']
        ind = self.config[key]['ind']

        return go.Scatter(x=self.xx,
                    y=self.array[:, ind, 0],
                    line=dict(color=clr, dash=dash),
                    name=name
                    )






    def add_traces_to_layout_model_output_overview(self, data_dict):
        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2)

        self.add_traces(fig, data_dict['S_R'], 1, 1)
        self.add_traces(fig, data_dict['E'], 1, 2)
        self.add_traces(fig, data_dict['I'], 2, 1)
        self.add_traces(fig, data_dict['F'], 2, 2)

        return fig



    def _sort_layout(self, fig):
        fig = self.update_axes(fig)


        fig.update_layout(standard_layout(False, self.width, self.height))
        
        fig = self.add_corner_text_labels(fig)

        fig = self.sort_legend(fig)

        return fig


    def sort_legend(self, fig):
        fig.update_layout(showlegend=True, legend=dict(font=dict(size=16)))
        return fig



    def add_corner_text_labels(self, fig):
        top_row = 1.08
        bottom_row = 0.5
        
        left = 0.01
        middle = 0.58

        c1 = get_big_text_annotation(left, top_row, 'A')
        c2 = get_big_text_annotation(middle, top_row, 'B')
        c3 = get_big_text_annotation(left, bottom_row, 'C')
        c4 = get_big_text_annotation(middle, bottom_row, 'D')
        
        annotz = [c1, c2, c3, c4]

        fig.update_layout(annotations=annotz)
        return fig




    def add_traces(self, fig, traces, row, col):
        for trace in traces:
            fig.add_trace(trace, row, col)

        return fig




    def update_axes(self, fig):
        fig.update_xaxes(row=1, col=1, showgrid=False)
        fig.update_xaxes(row=1, col=2, showgrid=False)
        
        fig.update_xaxes(title="Time (degree-days)", row=2, col=1, showgrid=False)
        fig.update_xaxes(title="Time (degree-days)", row=2, col=2, showgrid=False)
        
        fig.update_yaxes(title="L.A.I.", row=1, col=1)
        fig.update_yaxes(title="L.A.I.", row=1, col=2)
        fig.update_yaxes(title="L.A.I.", row=2, col=1)
        fig.update_yaxes(title="Concentration", row=2, col=2)

        return fig
    


    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/single/", "/paper_figs/model_overview_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)










class CombinedModelPlot:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 950

        self.xx = data['t_vec']
        
        self.array = data['sol_array']

        self.data = data

        self.DPC_year = 5

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)




    def _generate_figure(self):

        self.config = self.get_config()

        m_o_trcs_dict = self.get_model_output_overview_traces()

        y_rf_trcs_dict = self.get_yield_RF_traces()

        traces_dict = {**m_o_trcs_dict, **y_rf_trcs_dict}

        ugly_fig = self.add_traces_to_layout(traces_dict)

        fig = self._sort_layout(ugly_fig)

        return fig
        
    


    @staticmethod
    def get_config():

        E_cols = pltly_clrs.n_colors("rgb(0,0,255)", "rgb(255,255,255)", 5, colortype="rgb")[:4]
        
        I_cols = pltly_clrs.n_colors("rgb(255,0,0)", "rgb(255,255,255)", 5, colortype="rgb")[:4]


        return {
            'S': dict(color='limegreen', dash='solid', name='Susceptible', ind=PARAMS.S_ind),
            'R': dict(color='rgb(100,100,100)', dash='solid', name='Removed', ind=PARAMS.R_ind),

            'ERR': dict(color=E_cols[0], dash='dot', name='E (rr)', ind=PARAMS.ER_ind),
            'ERS': dict(color=E_cols[1], dash='dash', name='E (rs)', ind=PARAMS.ERS_ind),
            'ESR': dict(color=E_cols[2], dash='dashdot', name='E (sr)', ind=PARAMS.ESR_ind),
            'ESS': dict(color=E_cols[3], dash='solid', name='E (ss)', ind=PARAMS.ES_ind),

            'IRR': dict(color=I_cols[0], dash='dot', name='I (rr)', ind=PARAMS.IR_ind),
            'IRS': dict(color=I_cols[1], dash='dash', name='I (rs)', ind=PARAMS.IRS_ind),
            'ISR': dict(color=I_cols[2], dash='dashdot', name='I (sr)', ind=PARAMS.ISR_ind),
            'ISS': dict(color=I_cols[3], dash='solid', name='I (ss)', ind=PARAMS.IS_ind),

            'F1': dict(color='turquoise', dash='solid', name='Fungicide A', ind=PARAMS.Fung1_ind),
            'F2': dict(color='magenta', dash='dot', name='Fungicide B', ind=PARAMS.Fung2_ind),
        }





    def get_model_output_overview_traces(self):
        S_R_traces = self.get_S_R_traces()
        E_traces = self.get_E_traces()
        I_traces = self.get_I_traces()
        F_traces = self.get_F_traces()

        out = dict(S_R = S_R_traces,
                    E = E_traces,
                    I = I_traces,
                    F = F_traces)

        return out
    




    def get_S_R_traces(self):
        
        out = []

        for key in ['S', 'R']:
            out.append(self.get_DPC_trace(key))
        return out
    
    

    def get_E_traces(self):
        
        out = []

        for key in ['ERR', 'ERS', 'ESR', 'ESS']:
            out.append(self.get_DPC_trace(key))

        return out
    


    
    def get_I_traces(self):
        
        out = []

        for key in ['IRR', 'IRS', 'ISR', 'ISS']:
            out.append(self.get_DPC_trace(key))

        return out
    


    def get_F_traces(self):

        out = []

        for key in ['F1', 'F2']:
            out.append(self.get_DPC_trace(key))

        return out



    def get_DPC_trace(self, key):
        clr = self.config[key]['color']
        dash = self.config[key]['dash']
        name = self.config[key]['name']
        ind = self.config[key]['ind']

        return go.Scatter(x=self.xx,
                    y=self.array[:, ind, self.DPC_year],
                    line=dict(color=clr, dash=dash),
                    legendgroup="DPC",
                    name=name
                    )


    def get_yield_RF_traces(self):        

        yield_traces = self.get_yield_traces()
        
        RF_traces = self.get_RF_traces()
        
        traces = {"yield": yield_traces, "RF": RF_traces}

        return traces
    
    
    def get_yield_traces(self):
        out = []

        yy = self.data['yield_vec']
        
        xx = list(range(1,1+len(yy)))

        line = go.Scatter(x=xx, y=yy, name="Yield", legendgroup="YRF")
        
        Y_LOW = yy[-1]-2

        self.yield_lower_lim = Y_LOW
        
        X_END = 0.5 + xx[-1]
                
        shape = go.Scatter(x=[0, 0, X_END, X_END],
                            y=[Y_LOW, 95, 95, Y_LOW],
                            fill="toself",
                            mode="lines",
                            showlegend=False,
                            line=dict(width=0, color="rgb(150,150,150)"))
        
        out.append(shape)
        out.append(line)

        return out
    




    def get_RF_traces(self):
        out = []
        
        y1 = self.data['res_vec_dict']['f1']
        y2 = self.data['res_vec_dict']['f2']
        
        xx = list(range(len(y1)))

        for data, name, dash, col in zip([y1, y2], ['A', 'B'], ['solid', 'dot'], ['red', 'blue']):
            line = go.Scatter(x=xx, 
                y=data,
                name=f"R.F. (fung. {name})",
                legendgroup="YRF",
                line=dict(dash=dash, color=col))
            out.append(line)

        return out



    def add_traces_to_layout(self, data_dict):
        fig = make_subplots(rows=2, cols=3, 
                            horizontal_spacing=0.2,
                            shared_xaxes=True,
                            # row_heights=[0.3, 0.3]
                            )

        self.add_traces(fig, data_dict['S_R'], 1, 1)
        self.add_traces(fig, data_dict['E'], 1, 2)
        self.add_traces(fig, data_dict['I'], 2, 1)
        self.add_traces(fig, data_dict['F'], 2, 2)
        
        self.add_traces(fig, data_dict['yield'], 1, 3)
        self.add_traces(fig, data_dict['RF'], 2, 3)

        return fig



    def _sort_layout(self, fig):
        fig = self.update_axes(fig)

        fig.update_layout(standard_layout(False, self.width, self.height))
        
        annotz = self.get_corner_text_labels_comb()
        annotz += self.get_unacceptable_yield_annotation()
        
        fig.update_layout(annotations=annotz)

        fig = self.add_diagram(fig)

        fig.update_layout(margin=dict(t=450))

        fig = self.sort_legend(fig)

        return fig


    def add_diagram(self, fig):
        
        img = Image.open("runHRHR/img/diagram.png")

        fig.add_layout_image(
            dict(
                source=img,
                
                xref="paper", yref="paper",
                x=0,
                y=1.15,
                sizex=1.5,
                sizey=1.5,
                xanchor="left",
                yanchor="bottom",
                opacity=1,
                ))
        return fig


    def sort_legend(self, fig):
        fig.update_layout(showlegend=True, legend=dict(font=dict(size=16)))
        return fig


    @staticmethod
    def get_corner_text_labels_comb():

        very_top = 1.8
        
        top_row = 1.12
        bottom_row = 0.54
        
        left = -0.04
        middle = 0.33
        right = 0.72

        annotz = [
            get_big_text_annotation(left, very_top, 'A'),
            get_big_text_annotation(left, top_row, 'B'),
            get_big_text_annotation(left, bottom_row, 'C'),
            get_big_text_annotation(middle, top_row, 'D'),
            get_big_text_annotation(middle, bottom_row, 'E'),
            get_big_text_annotation(right, top_row, 'F'),
            get_big_text_annotation(right, bottom_row, 'G'),
            ]

        return annotz



    def get_unacceptable_yield_annotation(self):
        # 0.5*(95+self.yield_lower_lim)

        return [dict(
            xref="x3",
            # yref="y1",
            xanchor="left",
            # yanchor="top",
            x= 0.5,
            y= 0.55,
            text="Unacc.<br>yield",
            showarrow=False,
            )]



    def add_traces(self, fig, traces, row, col):
        for trace in traces:
            fig.add_trace(trace, row, col)

        return fig




    def update_axes(self, fig):

        fig.update_xaxes(row=1, col=1, showgrid=False)
        fig.update_xaxes(row=1, col=2, showgrid=False)
        
        fig.update_xaxes(title="Time<br>(degree-days)", row=2, col=1, showgrid=False, zeroline=False)
        fig.update_xaxes(title="Time<br>(degree-days)", row=2, col=2, showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Time (years)", row=2, col=3, showgrid=False, zeroline=False)
        fig.update_xaxes(row=1, col=3, showgrid=False, zeroline=False)
        
        
        fig.update_yaxes(title="L.A.I.", row=1, col=1, showgrid=False, zeroline=False)
        fig.update_yaxes(title="L.A.I.", row=1, col=2, showgrid=False, zeroline=False)
        fig.update_yaxes(title="L.A.I.", row=2, col=1, showgrid=False, zeroline=False)
        fig.update_yaxes(title="Concentration", row=2, col=2, showgrid=False, zeroline=False)

        fig.update_yaxes(title="Yield", row=1, col=3, showgrid=False, zeroline=False)
        fig.update_yaxes(title="R.F.", row=2, col=3, showgrid=False, zeroline=False)

        return fig
    


    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/single/", "/paper_figs/model_overview_combined_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)


# End of CombinedModelPlot














class DoseSpaceScenariosPlot:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 650

        self.data = data

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)



    def _generate_figure(self):
        traces = self._get_traces()

        ugly_fig = self._add_traces_to_figure(traces)

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        traces = []

        traces.append(self.get_ERFB_legend_entry())
        traces.append(self.get_EqSel_legend_entry())

        traces.append(self.get_FY_heatmap())
        
        traces.append(self.get_ERFB_contour())
        traces.append(self.get_EqSel_contour())
        
        traces.append(self.get_full_dose_point())
        traces.append(self.get_min_dose_point())
        

        return traces


    def get_EqSel_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    name="ES contour"
                    )


    def get_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="solid"),
                    name="ERFB contour"
                    )



    def get_FY_heatmap(self):
        FY = self.data['FY']

        xheat = np.linspace(0, 1, FY.shape[0])
        yheat = np.linspace(0, 1, FY.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = np.transpose(FY),
            colorscale = grey_colorscale(FY),
            colorbar = my_colorbar(TITLE_MAP["FY"])
            )

        return heatmap


    def get_ERFB_contour(self):
        z = EqualResFreqBreakdownArray(self.data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'solid')
        out['name'] = "Delta RFB"

        return out



    def get_EqSel_contour(self):
        z = EqualSelectionArray(self.data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"

        return out



    @staticmethod
    def get_full_dose_point():
        out = go.Scatter(x=[1],
                y=[1],
                marker=dict(color='red', size=16),
                marker_symbol='circle',
                mode='markers',
                name="Full dose",
                )

        return out



    def get_min_dose_point(self):
        FYs = np.transpose(self.data['FY'])

        x = np.linspace(0, 1, FYs.shape[0])
        
        minEqDoseELVec = np.asarray([float(FYs[i, i]) for i in range(FYs.shape[0])])

        ind = np.where(minEqDoseELVec>0)

        min_val = x[int(ind[0][0])]

        out = go.Scatter(x=[min_val],
                y=[min_val],
                marker=dict(color='green', size=16),
                mode='markers',
                marker_symbol='square',
                name="Min. dose",
                )

        return out




    def _add_traces_to_figure(self, traces):
        fig = go.Figure(data=traces, layout=standard_layout(True, self.width, self.height))
        return fig



    def _sort_layout(self, fig):
        fig = self._update_axes(fig)
        fig = self._update_legend(fig)
        return fig
    


    @staticmethod
    def _update_axes(fig):
        eps = 0.04
        fig.update_xaxes(title="Dose (fungicide A)", range=[0-eps,1+eps], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide B)", range=[0-eps,1+eps], showgrid=False, zeroline=False)
        return fig
    



    @staticmethod
    def _update_legend(fig):
        fig.update_layout(legend=dict(x=1.2, y=0.95))
        return fig
    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/grid/", "/paper_figs/dose_space_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)





class DosesScatterPlot:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 620

        self.data = data

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)



    def _generate_figure(self):
        traces = self._get_traces()

        ugly_fig = self._add_traces_to_figure(traces)

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        traces = []
        
        line = go.Scatter(x=[0,0],
                    y=[0,12.5],
                    line=dict(color='rgb(50,50,50)', dash='dot'),
                    mode="lines"
                    )
        
        traces.append(line)

        z = EqualResFreqBreakdownArray(self.data).array
        FYs = self.data['FY']

        x = np.asarray(z).flatten()
        y = np.asarray(FYs).flatten()

        dose_sum_cols = self.get_dose_sum_vec(z.shape, y)
        

        scatter = go.Scatter(
                x=x,
                y=y,
                mode="markers",
                text=dose_sum_cols,
                marker=dict(color=dose_sum_cols,
                    size=6,
                    line=dict(width=0.2,
                            color='black'),

                    colorbar=dict(title="Sum of doses"),
                    colorscale='Viridis',
                    showscale=True)
                )

        traces.append(scatter)

        return traces



    def get_dose_sum_vec(self, matrix_shape, FYs_flat):
        array = np.zeros(matrix_shape)

        for i in range(matrix_shape[0]):
            for j in range(matrix_shape[1]):
                array[i,j] = i+j
        
        array = array*(2/array[-1,-1])

        
        ds_cols = array.flatten()

        dose_sum_cols = [ds_cols[i] if FYs_flat[i] else "grey" for i in range(len(ds_cols))]

        return dose_sum_cols



    def _add_traces_to_figure(self, traces):
        fig = go.Figure(data=traces, layout=standard_layout(False, self.width, self.height))
        return fig



    def _sort_layout(self, fig):
        fig = self._update_axes(fig)
        return fig
    


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title=r"$\Delta_{RFB}$")
        fig.update_yaxes(title="Effective life")
        return fig

    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/grid/", "/paper_figs/doses_scatter_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)





class YieldAndRfPlot:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 620

        self.data = data

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)



    def _generate_figure(self):
        trace_dict = self.get_trace_dict()

        ugly_fig = self._add_traces_to_figure(trace_dict)

        fig = self._sort_layout(ugly_fig)

        return fig



    def get_trace_dict(self):        
        

        yield_traces = self.get_yield_traces()
        
        RF_traces = self.get_RF_traces()
        
        traces = {"yield": yield_traces, "RF": RF_traces}

        return traces





    def get_yield_traces(self):
        out = []

        yy = self.data['yield_vec']
        
        xx = list(range(1,1+len(yy)))

        line = go.Scatter(x=xx, y=yy, name="Yield")
        
        Y_LOW = yy[-1]-2

        self.yield_lower_lim = Y_LOW
        
        X_END = 0.5 + xx[-1]
                
        shape = go.Scatter(x=[0, 0, X_END, X_END],
                            y=[Y_LOW, 95, 95, Y_LOW],
                            fill="toself",
                            mode="lines",
                            showlegend=False,
                            line=dict(width=0, color="rgb(150,150,150)"))
        
        out.append(shape)
        out.append(line)

        return out
    




    def get_RF_traces(self):
        out = []
        
        y1 = self.data['res_vec_dict']['f1']
        y2 = self.data['res_vec_dict']['f2']
        
        xx = list(range(len(y1)))

        for data, name, dash in zip([y1, y2], ['A', 'B'], ['dot', 'solid']):
            line = go.Scatter(x=xx, 
                y=data,
                name=f"<br>Resistance<br>frequency<br>(fungicide {name})",
                line=dict(dash=dash))
            out.append(line)

        return out


    def _add_traces_to_figure(self, trace_dict):
        fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
        
        fig = self._add_yield_traces(fig, trace_dict['yield'])

        fig = self._add_RF_traces(fig, trace_dict['RF'])
        
        return fig



    @staticmethod
    def _add_yield_traces(fig, traces):
        for trace in traces:
            fig.add_trace(trace, row=1, col=1)
        return fig
    


    @staticmethod
    def _add_RF_traces(fig, traces):
        for trace in traces:
            fig.add_trace(trace, row=2, col=1)
        return fig






    def _sort_layout(self, fig):
        fig = self._update_layout(fig)
        fig = self._update_axes(fig)
        return fig
    
    
    
    def _update_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))
        
        text = self.get_unacceptable_yield_annotation()
        corners = self._get_corner_text_labels()
        annotz = corners + text
        
        fig.update_layout(annotations=annotz)
        return fig


    def get_unacceptable_yield_annotation(self):
        return [dict(
            xref="x1",
            yref="y1",
            xanchor="left",
            x= 1,
            y= 2 + 0.5*(95+self.yield_lower_lim),
            text="Unacceptable yield",
            showarrow=False,
            )]



    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title="Time (years)", row=2, col=1, showgrid=False, zeroline=False)
        
        fig.update_yaxes(title="Yield<br>(% of disease free)", row=1, col=1)
        
        fig.update_yaxes(title="Resistance<br>frequency", row=2, col=1)

        return fig



    def _get_corner_text_labels(self):
        top_row = 1.08
        bottom_row = 0.52
        
        left = 0.02

        cA = get_big_text_annotation(left, top_row, 'A')
        cB = get_big_text_annotation(left, bottom_row, 'B')
        
        out = [cA, cB]
        return out
    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/single/", "/paper_figs/yield_rf_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)







class ParamScanPlotMeVsHobb:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 800

        self.data = self._process_data(data)

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)
    



    def _process_data(self, data):

        req_cols = data[['run', 'RFB_maxContEL', 'EqSel_maxContEL',
                                 "best_value_RFB", "best_value_EqSel",
                                 "RS", "SR",
                                 "RR", "sr_prop",
                                #  "delta_1", "delta_2",
                                 "omega_1", "omega_2",
                                 ]]
        
        complete = req_cols.dropna()

        complete['successMetric'] = self._get_success_metric(complete)

        complete['IRFMetric'] = complete.apply(lambda x: log10(x['RS']) - log10(x['SR']), axis=1)
        
        complete['AsympMetric'] = complete.apply(lambda x: x['omega_1']/(x['omega_1']+x['omega_2']), axis=1)

        return complete



    @staticmethod
    def _get_success_metric(x):

        bv_ES = np.asarray(x['best_value_EqSel'])
        ES_EL = np.asarray(x['EqSel_maxContEL'])
        bv_RFB = np.asarray(x['best_value_RFB'])
        RBF_EL = np.asarray(x['RFB_maxContEL'])

        out = [100*max(bv_ES[i], ES_EL[i]) / max(bv_RFB[i], RBF_EL[i]) for i in range(len(bv_ES))]

        return out



        


    def _generate_figure(self):
        trace_dict = self._get_traces()

        ugly_fig = self._add_traces_to_figure(trace_dict)

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        out = {}

        data = self.data
        
        data['logRR'] = [log10(x) for x in data['RR']]        

        for key in ['IRFMetric', 'AsympMetric', 'logRR', 'sr_prop']:
            
            scatter = go.Scatter(x=data[key],
                y=data['successMetric'],
                mode='markers',
                marker=dict(opacity=0.2))

            out[key] = scatter

        return out




    def _add_traces_to_figure(self, trace_dict):
        fig = make_subplots(rows=2, cols=2, vertical_spacing=0.3)

        fig.add_trace(trace_dict['IRFMetric'], row=1, col=1)
        fig.add_trace(trace_dict['AsympMetric'], row=2, col=1)
        fig.add_trace(trace_dict['logRR'], row=1, col=2)
        fig.add_trace(trace_dict['sr_prop'], row=2, col=2)

        return fig



    def _sort_layout(self, fig):
        fig = self._update_axes(fig)
        fig = self._update_layout(fig)
        fig = self._add_corner_text_labels(fig)
        return fig


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title="Log difference in <br>single resistance frequencies", row=1, col=1)
        fig.update_xaxes(title="Asymptote metric", row=2, col=1)
        fig.update_xaxes(title="Double resistant<br>frequency (log scale)", row=1, col=2)
        fig.update_xaxes(title="Proportion of sex", row=2, col=2)
        
        fig.update_yaxes(title="Success metric", row=1, col=1)
        fig.update_yaxes(title="Success metric", row=2, col=1)

        return fig

    
    def _update_layout(self, fig):
        fig.update_layout(standard_layout(False, self.width, self.height))
        return fig


    def _add_corner_text_labels(self, fig):
        top_row = 1.06
        bottom_row = 0.42
        
        left = -0.01
        middle = 0.56

        c1 = get_big_text_annotation(left, top_row, 'A')
        c2 = get_big_text_annotation(middle, top_row, 'B')
        c3 = get_big_text_annotation(left, bottom_row, 'C')
        c4 = get_big_text_annotation(middle, bottom_row, 'D')
        
        annotz = [c1, c2, c3, c4]

        fig.update_layout(annotations=annotz)
        return fig
    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = f"../outputs/figures/paper_figs/param_scan_hobb_me_{conf_str}"
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)






class ParamScanPlotHighLowDose:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 800

        self.data = self._process_data(data)

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)
    



    @staticmethod
    def _process_data(data):

        req_cols = data[['run', 'RFB_minDS', 'RFB_maxDS', 
                            'max_dose_sums',
                            'min_dose_sums',
                            "RS", "SR", "RR", "sr_prop", "delta_1", "delta_2"]]
        
        complete = req_cols.dropna()

        complete['minMixStrength'] = complete.apply(lambda x: (
            100*(x['RFB_minDS']- x['min_dose_sums'])) / (x['max_dose_sums']
             - x['min_dose_sums'])
            , axis=1)
        
        complete['maxMixStrength'] = complete.apply(lambda x: (
            100*(x['RFB_maxDS']- x['min_dose_sums'])) / (x['max_dose_sums']
             - x['min_dose_sums'])
            , axis=1)
        
        
        complete['meanMixStrength'] = complete.apply(lambda x: (
            100*(x['RFB_minDS'] - x['min_dose_sums'] + 0.5*(x['RFB_maxDS'] - x['RFB_minDS']))) / (x['max_dose_sums']
             - x['min_dose_sums'])
            , axis=1)
        
        print(complete['meanMixStrength'])

        complete['IRFMetric'] = complete.apply(lambda x: log10(x['RS']) - log10(x['SR']), axis=1)
        
        complete['DecayMetric'] = complete.apply(lambda x: x['delta_1']/(x['delta_1']+x['delta_2']), axis=1)

        return complete



        


    def _generate_figure(self):
        trace_dict = self._get_traces()

        ugly_fig = self._add_traces_to_figure(trace_dict)

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        out = {}

        data = self.data
        
        # y = data['minMixStrength']
        # y = data['maxMixStrength']
        y = data['meanMixStrength']
        
        x1 = data['IRFMetric']
        # x1 = [log10(x) for x in data['SR']]

        x2 = data['DecayMetric']

        x3 = [log10(x) for x in data['RR']]
        
        x4 = data['sr_prop']
        
        out['IRFMetric'] = [self._get_scatter(x1, y)]

        out['DecayMetric'] = [self._get_scatter(x2, y)]
        
        out['RR'] = [self._get_scatter(x3, y)]

        out['sr_prop'] = [self._get_scatter(x4, y),
                            self.get_best_fit(x4, y, "blue")
                            ]

        return out




    def _get_scatter(self, x, y):
        return go.Scatter(x=x,
            y=y,
            mode='markers',
            marker=dict(opacity=0.2))
        


    @staticmethod
    def get_best_fit(xx, y, col):
        slope, intercept, _, _, _ = stats.linregress(xx,y)

        line = slope*np.asarray(xx) + intercept
        
        best_fit = go.Scatter(
                    x=xx,
                    y=line,
                    mode='lines',
                    marker=dict(color=col),
                    name='Fit'
                    )

        return best_fit



    def _add_traces_to_figure(self, trace_dict):
        fig = make_subplots(rows=2, cols=2, vertical_spacing=0.3)
        
        self.add_traces_each_subplot(fig, 1, 1, trace_dict['IRFMetric'])
        self.add_traces_each_subplot(fig, 2, 1, trace_dict['DecayMetric'])
        self.add_traces_each_subplot(fig, 1, 2, trace_dict['RR'])
        self.add_traces_each_subplot(fig, 2, 2, trace_dict['sr_prop'])
        
        return fig





    @staticmethod
    def add_traces_each_subplot(fig, row, col, traces):
        for trace in traces:
            fig.add_trace(trace, row=row, col=col)
        return fig

    def _sort_layout(self, fig):
        fig = self._update_axes(fig)
        fig = self._update_layout(fig)
        fig = self._add_corner_text_labels(fig)
        return fig
    


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title="Log difference in <br>single resistance frequencies", row=1, col=1)
        fig.update_xaxes(title="Decay rate metric", row=2, col=1)
        fig.update_xaxes(title="Double resistant<br>frequency (log scale)", row=1, col=2)
        fig.update_xaxes(title="Proportion of sex", row=2, col=2)
        
        fig.update_yaxes(title="Optimal distance<br>along contour (%)", row=1, col=1)
        fig.update_yaxes(title="Optimal distance<br>along contour (%)", row=2, col=1)

        return fig

    
    def _update_layout(self, fig):
        fig.update_layout(standard_layout(False, self.width, self.height))
        return fig


    def _add_corner_text_labels(self, fig):
        top_row = 1.06
        bottom_row = 0.42
        
        left = -0.01
        middle = 0.56

        c1 = get_big_text_annotation(left, top_row, 'A')
        c2 = get_big_text_annotation(middle, top_row, 'B')
        c3 = get_big_text_annotation(left, bottom_row, 'C')
        c4 = get_big_text_annotation(middle, bottom_row, 'D')
        
        annotz = [c1, c2, c3, c4]

        fig.update_layout(annotations=annotz)
        return fig
    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = f"../outputs/figures/paper_figs/param_scan_high_low_{conf_str}"
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)









class FigTemplateClass:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 620

        self.data = data

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)



    def _generate_figure(self):
        traces = self._get_traces()

        ugly_fig = self._add_traces_to_figure(traces)

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        traces = []

        return traces


    def _add_traces_to_figure(self, traces):
        fig = go.Figure(data=traces, layout=standard_layout(False, self.width, self.height))
        return fig



    def _sort_layout(self, fig):
        fig = self._update_axes(fig)
        return fig
    


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="y")
        return fig

    
    
    
    
    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/grid/", "/paper_figs/doses_scatter_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)
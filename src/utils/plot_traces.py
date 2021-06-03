
from math import ceil
import numpy as np
import plotly.graph_objects as go

from .plot_consts import STRAIN_ATTRS, TITLE_MAP
from .plot_utils import invisible_colorbar, my_colorbar, grey_colorscale
from .functions import logit10, log10_difference, logit10_difference

# * RFB

def _update_RFB_x_y(ind, data, RFB_diff):
    x = []
    y = []

    FY = data["FY"]

    for i in range(ind+1):
        j = ind-i

        if i<FY.shape[0] and j<FY.shape[1]:
            fy = int(FY[i,j])
            
            if fy>0:
                rr1 = data['res_arrays']['f1'][i,j,fy]
                rr2 = data['res_arrays']['f2'][i,j,fy]

                rf_diff_breakdown = logit10_difference(rr1, rr2)

                x.append(rf_diff_breakdown)
                y.append(fy)

                RFB_diff[j, i] = rf_diff_breakdown

            else:
                RFB_diff[j, i] = None
    
    return x, y, RFB_diff


def _get_line(x, y, color, dash_):
    return go.Scatter(x=x,
                        y=y,
                        showlegend=False,
                        line=dict(color=color, dash=dash_))
    

def get_RFB_diff_traces(data, z, N_y_int, inds_list, colors):
    RFB_diff = np.zeros(z.shape)
    traces_RFB = []
    
    for ind in range(N_y_int):
        
        x, y, RFB_diff = _update_RFB_x_y(ind, data, RFB_diff)
        
        if not x or not (ind in inds_list):
            continue
        
        myline = _get_line(x, y, colors[ind], "solid")
        traces_RFB.append(myline)
        
    return traces_RFB, RFB_diff

# * End of RFB


# * FY selection


def _get_sel_one_strain(data, strain, i, j):
    y0 = data['start_freqs'][strain][i,j,0]
    y1 = data['start_freqs'][strain][i,j,1]
    return y1/y0


def _get_fy_sel_diff(data, i, j):
    s1_y1 = _get_sel_one_strain(data, 'RS', i, j)
    s2_y1 = _get_sel_one_strain(data, 'SR', i, j)
    return log10_difference(s1_y1, s2_y1)



def _update_x_y_fy_sel(ind, data, FY, fy_sel):
    x = []
    y = []

    for i in range(ind+1):
        
        j = ind-i

        if i<FY.shape[0] and j<FY.shape[1]:
            fy = int(FY[i,j])
            
            if fy>0:
                fy_selection = _get_fy_sel_diff(data, i, j)

                x.append(fy_selection)
                y.append(fy)

                fy_sel[j, i] = fy_selection
            else:
                fy_sel[j, i] = None
            
    return x, y, fy_sel



def get_eq_sel_traces(data, z, N_y_int, inds_list, colors):
    
    fy_sel = np.zeros(z.shape)
    traces_sel = []
    
    FY = data["FY"]

    for ind in range(N_y_int):
        
        x, y, fy_sel = _update_x_y_fy_sel(ind, data, FY, fy_sel)
        
        if not x or not (ind in inds_list):
            continue
        
        line = _get_line(x, y, colors[ind], "dot")
        traces_sel.append(line)
    
    return traces_sel, fy_sel

# * End of FY selection


# * Heatmap lines

def _get_hm_line_col(ind, min_ind, max_ind):
    mn = 0
    mx = 220

    clr = mn + (mx-mn)* (ind - min_ind - 1)/(max_ind - min_ind - 1)
    
    return f"rgba({255-clr},{0},{255-clr},0.8)"

def _get_color_range(inds_list):
    colors = {}

    ind_min = min(inds_list)
    ind_max = max(inds_list)
    
    for ind in inds_list:

        colors[ind] = _get_hm_line_col(ind, ind_min, ind_max)
    
    return colors

def _get_inds_list(N_y_int, ind0):

    n_lines = 4
    inds_list = []

    for ind in range(ind0+1, N_y_int):

        n_interval = ceil((N_y_int - ind0)/n_lines)
        # only actually want n_lines lines total
        if not (((ind - 1 -ind0) % n_interval == 0 and
                        (N_y_int - 1 - ind > n_interval))
                    or (ind==N_y_int-1)):
            # only want a line every n_interval steps...
            # and want final point but then not another point v nearby
            continue

        inds_list.append(ind)

    return inds_list


def _get_ind0(N_y_int, z):
    ind0 = 0
    for ind in range(N_y_int):

        dose_line_all_0 = True
        
        for i in range(ind):
            j = ind-i
            if i<z.shape[0] and j<z.shape[1] and z[i,j]>0:
                dose_line_all_0 = False

        if dose_line_all_0:
            ind0 = ind
            continue
    
    return ind0



def _get_hm_lines(y_intrcpt, inds_list, colors, N_d):
    out = []
    for ind in inds_list:
        
        xx = y_intrcpt[:ind+1]
        xx = [x for x in xx if (x<=1 and x>=0)]
        yy = [y_intrcpt[ind] - x for x in xx]

        if len(yy)==N_d:
            # 1.001 since floating point error caused probs :/
            yy = [y for y in yy if (y<=1.001)]
            xx = xx[len(xx)-len(yy):]

        ds = round(y_intrcpt[ind], 2)
        
        scat = go.Scatter(
                x = xx,
                y = yy,
                line=dict(color=colors[ind]),
                name=f"Dose sum: {ds}",
            )
        
        out.append(scat)

    return out



def get_heatmap_lines(Config, z, y_intrcpt):

    N_y_int = len(y_intrcpt)
    ind0 = _get_ind0(N_y_int, z)
    inds_list = _get_inds_list(N_y_int, ind0)

    colors = _get_color_range(inds_list)
    traces = _get_hm_lines(y_intrcpt, inds_list, colors, Config.n_doses)

    return traces, colors, inds_list


# * End of Heatmap lines


# * Strain freq traces

def _get_strain_freq_line(n_yr_run, key, rf_s, rf_e, season_frac):
    x = []
    y = []

    for i in range(n_yr_run):
        x.append(i)
        y.append(logit10(rf_s[key][i]))

        if i!=n_yr_run-1:
            x.append(i+season_frac)
            y.append(logit10(rf_e[key][i]))
    
    line = go.Scatter(x=x,
                    y=y,
                    mode="lines",
                    name=STRAIN_ATTRS[key]['abbrv'],
                    line=dict(color=STRAIN_ATTRS[key]['color'],
                            dash=STRAIN_ATTRS[key]['dash'])
                    )
    return line


def _get_strain_freq_scat(n_yr_run, key, rf_s):
    y_scat = []
    
    for i in range(n_yr_run):
        y_scat.append(logit10(rf_s[key][i]))

    x_scat = list(range(len(y_scat)))
    
    scatter = go.Scatter(x=x_scat,
                    y=y_scat, 
                    line=dict(color=STRAIN_ATTRS[key]['color']),
                    showlegend=False,
                    mode="markers",
                    )
    return scatter


def _get_strain_freq_shape(n_yr_run, season_frac, shape_min, shape_max):
    shapes = []
    for i in range(n_yr_run-1):
        shapes.append(go.Scatter(x=[i+season_frac, i+1, i+1, i+season_frac],
                        y=[shape_min, shape_min, shape_max, shape_max],
                        fill="toself",
                        mode="lines",
                        showlegend=False,
                        line=dict(width=0, color="rgb(200,200,255)")))
    return shapes


def _get_shape_min_max_y(rf_s, rf_e, n_yr_run, keys):
    min_ = 1
    max_ = 0
    for key in keys:
        min_ = min(np.amin(rf_s[key][:n_yr_run]), np.amin(rf_e[key][:n_yr_run]), min_)
        max_ = max(np.amax(rf_s[key][:n_yr_run]), np.amax(rf_e[key][:n_yr_run]), max_)
    
    shape_min = logit10(min_) - 0.2
    shape_max = logit10(max_) + 0.2
    
    return shape_min, shape_max



def get_strain_freq_traces(rf_s, rf_e):
    n_yr_run = len(rf_s['RR']) - 1
    
    keys = list(rf_s.keys())
    keys.reverse()

    shape_min, shape_max = _get_shape_min_max_y(rf_s, rf_e, n_yr_run, keys)
    
    season_frac = 0.75
    
    traces = []

    for key in keys:
        line = _get_strain_freq_line(n_yr_run, key, rf_s, rf_e, season_frac)
        shapes = _get_strain_freq_shape(n_yr_run, season_frac, shape_min, shape_max)
        scatter = _get_strain_freq_scat(n_yr_run, key, rf_s)
        
        traces.append(line)
        traces.append(scatter)
    
    traces = shapes + traces
    
    return traces

# * End of Strain freq traces


def contour_at_0(x, y, z, color, dash):
    return go.Contour(x=x,
                    y=y,
                    z=z,
                    contours=dict(start=0, end=0),
                    contours_coloring='lines',
                    colorscale=[color]*2,
                    line=dict(width=2, dash=dash),
                    colorbar=invisible_colorbar(0.42),
                    )


def contour_at_single_level(x, y, z, level, color, dash):
    return go.Contour(x=x,
                    y=y,
                    z=z,
                    contours=dict(start=level, end=level),
                    contours_coloring='lines',
                    colorscale=[color]*2,
                    line=dict(width=2, dash=dash),
                    colorbar=invisible_colorbar(0.42),
                    )


# Eq RFB contours

def _multi_contours(x, y, z, cont_list, colors):
    out = []

    for i in range(len(cont_list)):
        color = colors[i]

        out.append(go.Contour(x=x,
                    y=y,
                    z=z,
                    contours=dict(start=cont_list[i],
                        end=cont_list[i],
                        showlabels = True,
                        ),
                    contours_coloring='lines',
                    colorscale=[color]*2,
                    line=dict(width=2, dash="solid"),
                    colorbar=invisible_colorbar(0.42),
                    ))
    return out


def get_multi_contour_traces(data, Config):
    # cont_list = [-4, -1, -0.1, 0, 0.1, 1]
    # cont_list = [-1, 0, 1]
    cont_list = [0]
    
    xheat = np.linspace(0, 1, Config.n_doses)
    yheat = np.linspace(0, 1, Config.n_doses)
    z = np.transpose(data["FY"])
    
    traces = []
    
    clrbar = my_colorbar(TITLE_MAP["FY"])
    
    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=grey_colorscale(z),
        colorbar=clrbar
    )
    
    traces.append(heatmap)

    colors = _get_color_range(list(range(len(cont_list))))

    N_y_int = 2*Config.n_doses-1
    _, RBF_diff = get_RFB_diff_traces(data, z, N_y_int, [], colors)
    
    traces += _multi_contours(xheat, yheat, RBF_diff, cont_list, colors)
    
    return traces

# End of RFB contours

# MRFB contours

def _add_contour_lines(contours, names, colors):
    
    out = []
    for i, cont in enumerate(contours):
        line = go.Scatter(x=cont[0],
                    y=cont[1],
                    name=names[i],
                    line=dict(color=colors[i]),
                    mode="lines+markers")
        out.append(line)
    return out


def get_MRFB_contour_traces(grid, contours, names, Config):
    
    xheat = np.linspace(0, 1, Config.n_doses)
    yheat = np.linspace(0, 1, Config.n_doses)
    z = np.transpose(grid["FY"])
    
    clrbar = my_colorbar(TITLE_MAP["FY"])
    
    heatmap = go.Heatmap(
        x = xheat,
        y = yheat,
        z = z,
        colorscale=grey_colorscale(z),
        colorbar=clrbar
        )
    
    traces = [heatmap]

    colors = _get_color_range(list(range(len(names))))
    
    traces += _add_contour_lines(contours, names, colors)

    return traces

# End of MRFB contours


# MS_RFB_scatter
def _get_MS_RFB_FY_df(data, ind):
    rfb_list = []
    ms_list = []
    fy_list = []

    FY = data["FY"]
            
    for i in range(ind):
        j = ind - i - 1

        if i<FY.shape[0] and j<FY.shape[1]:
            fy = int(FY[i,j])
            
            if fy>0:
                rr1 = data['res_arrays']['f1'][i,j,fy]
                rr2 = data['res_arrays']['f2'][i,j,fy]

                rf_diff_breakdown = logit10_difference(rr1, rr2)

                rfb_list.append(rf_diff_breakdown)
                fy_list.append(fy)
                ms_list.append(2*ind/(2*FY.shape[0]-1))

    return rfb_list, ms_list, fy_list




def _get_clr(x, limits, colors):
    if x>=limits[-1]:
        out = colors[-1]
    elif x>=limits[-2] and x<limits[-1]:
        out = colors[-2]
    elif x>=limits[-3] and x<limits[-2]:
        out = colors[-3]
    elif x>=limits[-4] and x<limits[-3]:
        out = colors[-4]
    else:
        out = colors[-5]
    return out


def _color_map(lst, limits, colors):
    out = [""]*len(lst)
    
    for i in range(len(lst)):
        out[i] = _get_clr(lst[i], limits, colors)

    
    return out


def _get_legend_trace(col, name):
    return go.Scatter(x=[0],
                        y=[2.5],
                        marker=dict(color=col),
                        name=name,
                        mode="markers",
                        )

def _get_legend_traces(colors, limits):
    
    traces = []
    for i in range(1,len(limits)+1):
        
        # if i==0:
            # name = f"EL<{str(int(limits[i]))}"
        if i==len(limits):
            name = "EL" + u"\u2265" + f"{str(int(limits[i-1]))}"
        else:
            name = f"{str(int(limits[i-1]))}" + u"\u2264" +  f"EL<{str(int(limits[i]))}"

        traces.append(
        _get_legend_trace(colors[i], name)
        )
        
    return traces


def get_MS_RFB_traces(data):
    
    FY = data["FY"]
    N_lim = ceil(0.5*(1+np.amax(FY)-np.amin(FY[FY>0])))
    limits = np.linspace(np.amin(FY[FY>0]), np.amax(FY), N_lim)
    
    # how many points?
    # pos = np.sum(np.array(FY) > 0, axis=0)
    # print(sum(pos))

    
    # colors = ["rgb(0,0,0)",
    #         "rgb(100,100,100)",
    #         "rgb(200,200,200)",
    #         "orange",
    #         "red",
    #         "blue"]

    clrs = _get_color_range(list(range(len(limits)+1)))
    colors = [clrs[key] for key in clrs.keys()]

    colors[-3:] = ['orange', 'red', 'blue']

    traces = []
    
    N = 2*FY.shape[0]

    for i in range(N):
        rfb_list, ms_list, fy_list = _get_MS_RFB_FY_df(data, i)
        if rfb_list:
            traces.append(go.Scatter(x=rfb_list,
                    y=ms_list,
                    marker=dict(color=_color_map(fy_list, limits, colors)),
                    mode="markers",
                    showlegend=False,
                    ))
    
    traces += _get_legend_traces(colors, limits)

    return traces
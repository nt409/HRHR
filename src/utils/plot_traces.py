
from math import log10, ceil
import numpy as np
import plotly.graph_objects as go

from .plot_consts import STRAIN_ATTRS
from .plot_utils import invisible_colorbar

def logit10(x):
    if x>0 and x<1:
        return log10(x/(1-x))
    else:
        raise Exception(f"x={x} - invalid value")

def log10_difference(x1, x2):
    return log10(x1) - log10(x2)

def logit10_difference(x1, x2):
    return logit10(x1) - logit10(x2)

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

def _get_hm_line_col(ind, y_int_N, ind0):
    mn = 0
    mx = 220

    clr = mn + (mx-mn)* (ind - 1 - ind0)/(y_int_N - ind0 - 1)
    
    return f"rgba({255-clr},{0},{255-clr},0.8)"

def _get_hm_colors(inds_list, ind0, N_y_int):
    colors = {}
    
    for ind in inds_list:

        colors[ind] = _get_hm_line_col(ind, N_y_int, ind0)
    
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

    colors = _get_hm_colors(inds_list, ind0, N_y_int)
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
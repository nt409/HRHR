"""Figures for the paper."""

from __future__ import annotations
from abc import ABC, abstractmethod
import copy
from numpy.core.fromnumeric import size
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
# import pandas as pd
from math import floor, log10
from scipy import stats
from PIL import Image
import plotly.express as px
import random


from model.utils import logit10
from model.strategy_arrays import EqualResFreqBreakdownArray, EqualSelectionArray

from plotting.traces import contour_at_0, contour_at_single_level

from plotting.utils import divergent_color_scale, get_shape_annotation, get_text_annotation, get_arrow_annotation,\
    grey_colorscale_discrete, grey_colorscale_discrete_N, \
    invisible_colorbar, my_colorbar_subplot, standard_layout, \
    grey_colorscale, my_colorbar, get_big_text_annotation

from plotting.consts import ATTRS_DICT, LABEL_COLOR, LIGHT_GREY_TEXT, NULL_HEATMAP_COLOUR, TITLE_MAP, PLOT_WIDTH, PLOT_HEIGHT, \
        FULL_PAGE_WIDTH


# Paper Figs
class BasicFig(ABC):
    def __init__(self) -> None:
        pass

    @abstractmethod 
    def _generate_figure(self):
        pass

    @abstractmethod 
    def _sort_layout(self):
        pass

    def _save_and_show(self, fig):
        fig.show()
        
        print(f"saving figure to: \n{self.filename}")

        fig.write_image(self.filename)












class CombinedModelPlot(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 950
        
        self.states_list = data.states_list

        self.data = data

        self.DPC_year = 5

        self.xx = data.states_list[self.DPC_year].t

        fig = self._generate_figure()
        
        self.filename = conf_str.replace("/single/", "/paper_figs/model_overview_combined_")

        self._save_and_show(fig)




    def _generate_figure(self):

        m_o_trcs_dict = self.get_model_output_overview_traces()

        y_rf_trcs_dict = self.get_yield_RF_traces()

        traces_dict = {**m_o_trcs_dict, **y_rf_trcs_dict}

        ugly_fig = self.add_traces_to_layout(traces_dict)

        fig = self._sort_layout(ugly_fig)

        return fig
        
    




    def get_model_output_overview_traces(self):
        S_traces = self.get_DPC_traces(['S'])
        E_traces = self.get_DPC_traces(['ERR', 'ERS', 'ESR', 'ESS'])
        I_traces = self.get_DPC_traces(['IRR', 'IRS', 'ISR', 'ISS'])
        R_traces = self.get_DPC_traces(['R'])

        out = dict(S = S_traces,
                    E = E_traces,
                    I = I_traces,
                    R = R_traces)

        return out
    


    
    def get_DPC_traces(self, keys):
        out = []
        for key in keys:
            out.append(self.get_DPC_trace(key))
        return out
    




    def get_DPC_trace(self, key):
        clr = ATTRS_DICT[key]['color']
        dash = ATTRS_DICT[key]['dash']
        name = ATTRS_DICT[key]['name']
        group = ATTRS_DICT[key]['legendgrp']

        xx = (1e-3)*self.xx
        yy = vars(self.states_list[self.DPC_year])[key]
        
        if group in ["E", "I"]:
            yy = (1e3)*yy

        return go.Scatter(x=xx,
                    y=yy,
                    line=dict(color=clr, dash=dash),
                    legendgroup=group,
                    name=name
                    )


    def get_yield_RF_traces(self):        

        yield_traces = self.get_yield_traces()
        
        RF_traces = self.get_RF_traces()
        
        traces = {"yield": yield_traces, "RF": RF_traces}

        return traces
    
    
    def get_yield_traces(self):
        out = []

        yy = self.data.yield_vec
        
        xx = list(range(1,1+len(yy)))

        line = go.Scatter(x=xx, y=yy, name="Yield", line=dict(color="darkgreen"), legendgroup="Y")
        
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
        
        y1 = self.data.res_vec_dict['f1']
        y2 = self.data.res_vec_dict['f2']
        
        xx = list(range(len(y1)))

        for data, name, dash, col in zip([y1, y2], ['A', 'B'], ['solid', 'dot'], ['red', 'blue']):
            line = go.Scatter(x=xx, 
                y=data,
                name=f"Resistance<br>frequency<br>(fungicide <i>{name}</i>)",
                legendgroup="RF",
                line=dict(dash=dash, color=col))
            out.append(line)

        return out



    def add_traces_to_layout(self, data_dict):
        fig = make_subplots(rows=2, cols=3, 
                            horizontal_spacing=0.25,
                            shared_xaxes=True,
                            # row_heights=[0.3, 0.3]
                            )

        self.add_traces(fig, data_dict['S'], 1, 1)
        self.add_traces(fig, data_dict['R'], 2, 1)
        self.add_traces(fig, data_dict['E'], 1, 2)
        self.add_traces(fig, data_dict['I'], 2, 2)
        
        self.add_traces(fig, data_dict['yield'], 1, 3)
        self.add_traces(fig, data_dict['RF'], 2, 3)

        return fig



    def _sort_layout(self, fig):
        fig = self.update_axes(fig)

        fig.update_layout(standard_layout(True, self.width, self.height))
        
        annotz = self.get_corner_text_labels_comb()
        annotz += self.get_unacceptable_yield_annotation()
        
        fig.update_layout(annotations=annotz)

        fig = self.add_diagram(fig)

        fig.update_layout(margin=dict(t=450))

        fig.update_layout(legend=dict(
            font=dict(size=14),
            x=1.03,
            y=1.10,
            xanchor="left",
            yanchor="top",
            ))

        return fig


    def add_diagram(self, fig):
        
        img = Image.open("figs/img/diagram.png")

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




    @staticmethod
    def get_corner_text_labels_comb():

        very_top = 1.8
        
        top_row = 1.11
        bottom_row = 0.53
        
        left = -0.04
        middle = 0.35
        right = 0.76

        annotz = [
            get_big_text_annotation(left, very_top, 'A'),
            
            get_big_text_annotation(left, top_row, 'B'),
            get_big_text_annotation(middle, top_row, 'C'),
            get_big_text_annotation(left, bottom_row, 'D'),
            get_big_text_annotation(middle, bottom_row, 'E'),

            get_big_text_annotation(right, top_row, 'F'),
            get_big_text_annotation(right, bottom_row, 'G'),
            ]

        return annotz



    def get_unacceptable_yield_annotation(self):
        return [dict(
            xref="x3",
            xanchor="left",
            font=dict(size=16),
            textangle=-90,
            x= 0.35,
            y= 0.4,
            text="Unacceptable<br>yield",
            align="left",
            showarrow=False,
            )]



    def add_traces(self, fig, traces, row, col):
        for trace in traces:
            fig.add_trace(trace, row, col)

        return fig




    def update_axes(self, fig):

        fig.update_xaxes(row=1, col=1, showgrid=False)
        fig.update_xaxes(row=1, col=2, showgrid=False)
        fig.update_xaxes(row=1, col=3, showline=True, showgrid=False)
        
        fig.update_xaxes(title="Time<br>(degree-days x10<sup>3</sup>)", row=2, col=1, showgrid=False)
        fig.update_xaxes(title="Time<br>(degree-days x10<sup>3</sup>)", row=2, col=2, showgrid=False)
        fig.update_xaxes(title="Time (years)", row=2, col=3, showgrid=False)
        
        
        fig.update_yaxes(title="Leaf proportion", row=1, col=1, showline=True, showgrid=False)
        fig.update_yaxes(title="Leaf proportion<br>(x10<sup>-3</sup>)", row=1, col=2, showline=True, showgrid=False)

        fig.update_yaxes(title="Leaf proportion", row=2, col=1, showline=True, showgrid=False)
        fig.update_yaxes(title="Leaf proportion<br>(x10<sup>-3</sup>)", row=2, col=2, showline=True, showgrid=False)

        fig.update_yaxes(title="Yield", row=1, col=3, showgrid=False)
        fig.update_yaxes(title="Resistance<br>frequency", row=2, col=3, showgrid=False)

        return fig
    


# End of CombinedModelPlot














class DoseSpaceOverview(BasicFig):
    def __init__(self, data, contour_data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 700

        self.data = data

        self.contour_data = contour_data

        self.opt_point_index = dict(ERFB = int(5), ESFY = int(8))

        self.cbar_attrs = dict(
                        x=[0.375, 1.005],
                        y=[0.21, 0.79],
                        # intercepts are where the blue/black lines cross 0.5/0 resp.
                        intercepts=[0.238, 0.765],
                        len=0.46)

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_overview_")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.28)

        fig.add_traces(self._get_EL_traces(), rows=1, cols=1)
        fig.add_traces(self._get_ESFY_traces(), rows=1, cols=2)

        fig.add_traces(self._get_ERFB_traces(), rows=2, cols=1)
        fig.add_traces(self._get_contour_traces(), rows=2, cols=2)

        fig = self._sort_layout(fig)

        return fig


    # A
    def _get_EL_traces(self):        
        traces = [

            self._get_ERFB_legend_entry(),
            self._get_ESFY_legend_entry(),

            self._get_FY_trace(),
            
            self._get_ERFB_contour_single(),
            self._get_ESFY_contour_single(),
            ]

        traces += self._get_ERFB_cont_scatter_trcs_heatmap_plot()
        traces += self._get_ESFY_cont_scatter_trcs_heatmap_plot()

        return traces
    
    # B
    def _get_ESFY_traces(self):
        traces = [
            self.get_grey_heatmap(),
            self.get_ESFY_trace(),
            self._get_ESFY_contour_single(),
            ]
        traces += self._get_ESFY_cont_scatter_trcs_heatmap_plot()
        return traces
    
    # C
    def _get_ERFB_traces(self):
        traces = [
            self.get_grey_heatmap(),
            self.get_ERFB_trace(),
            self._get_ERFB_contour_single(),
            ]
        traces += self._get_ERFB_cont_scatter_trcs_heatmap_plot()
        
        return traces    
    
    # D
    def _get_contour_traces(self):
        traces = [
            self._get_ERFB_cont_line(),
            self._get_ESFY_cont_line(),
            ]

        traces += self._get_ESFY_cont_scatter_trcs_DS()
        traces += self._get_ERFB_cont_scatter_trcs_DS()
        return traces







    def _get_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    legendgroup="cont",
                    name=u"\u03A9<sub>SFY</sub>=0.5 contour"
                    )


    def _get_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="dash"),
                    legendgroup="cont",
                    name=u"\u0394<sub>RFB</sub>=0 contour"
                    )



    def _get_FY_trace(self):
        FYs = np.transpose(self.data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            colorscale = grey_colorscale_discrete(FYs),
            colorbar = my_colorbar_subplot("E.L.", self.cbar_attrs['x'][0], self.cbar_attrs['y'][1], self.cbar_attrs['len'])
            )

        return heatmap
    

    def get_ERFB_trace(self):
        z = EqualResFreqBreakdownArray(self.data).array
        z_transpose = np.transpose(z)

        xheat = np.linspace(0, 1, z.shape[0])
        yheat = np.linspace(0, 1, z.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = z_transpose,
            colorscale = divergent_color_scale(z_transpose, 0),
            colorbar = my_colorbar_subplot(u"\u0394<sub>RFB", self.cbar_attrs['x'][0], self.cbar_attrs['y'][0], self.cbar_attrs['len'])
            )

        return heatmap


    def get_ESFY_trace(self):
        z = EqualSelectionArray(self.data).array
        z_transpose = np.transpose(z)

        xheat = np.linspace(0, 1, z.shape[0])
        yheat = np.linspace(0, 1, z.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = z_transpose,
            colorscale = divergent_color_scale(z_transpose, 0.5),
            colorbar = my_colorbar_subplot(u"\u03A9<sub>SFY", self.cbar_attrs['x'][1], self.cbar_attrs['y'][1], self.cbar_attrs['len'])
            )

        return heatmap


    def _get_ERFB_cont_line(self):
        x = self.contour_data['ERFB']['DS']
        y = self.contour_data['ERFB']['FY']
        return go.Scatter(x=x, y=y, mode="lines",
                        line=dict(color="black", dash="dash"),
                        showlegend=False)

    def _get_ESFY_cont_line(self):
        x = self.contour_data['ESFY']['DS']
        y = self.contour_data['ESFY']['FY']
        return go.Scatter(x=x, y=y, mode="lines", 
                        line=dict(color="blue", dash="dot"),
                        showlegend=False)
    


    def _get_ERFB_cont_scatter_trcs_DS(self):
        x = self.contour_data['ERFB']['DS']
        y = self.contour_data['ERFB']['FY']

        opt_ind = self.opt_point_index['ERFB']

        names = ["ERFB: best doses", "ERFB: highest doses"]
        
        xs = [x[opt_ind], x[-1]]
        ys = [y[opt_ind], y[-1]]
        return self.scatter_point_pair(xs, ys, "black", names, True)

        
    
    def _get_ESFY_cont_scatter_trcs_DS(self):
        x = self.contour_data['ESFY']['DS']
        y = self.contour_data['ESFY']['FY']

        opt_ind = self.opt_point_index['ESFY']
        
        names = ["ESFY: best doses", "ESFY: highest doses"]

        xs = [x[opt_ind], x[-1]]
        ys = [y[opt_ind], y[-1]]
        return self.scatter_point_pair(xs, ys, "blue", names, True)
    


    def _get_ERFB_cont_scatter_trcs_heatmap_plot(self):
        x = self.contour_data['ERFB']['x']
        y = self.contour_data['ERFB']['y']

        opt_ind = self.opt_point_index['ERFB']

        xs = [x[opt_ind], x[-1]]
        ys = [y[opt_ind], y[-1]]
        return self.scatter_point_pair(xs, ys, "black")


    def _get_ESFY_cont_scatter_trcs_heatmap_plot(self):
        x = self.contour_data['ESFY']['x']
        y = self.contour_data['ESFY']['y']

        opt_ind = self.opt_point_index['ESFY']

        xs = [x[opt_ind], x[-1]]
        ys = [y[opt_ind], y[-1]]
        return self.scatter_point_pair(xs, ys, "blue")
    



    def scatter_point_pair(self, xs, ys, col, names=None, showlegend=False):
        # marker_size = 10

        if names is None:
            names = [None, None]

        trcs = [
            go.Scatter(x=[xs[0]], y=[ys[0]], mode="markers",
                marker=dict(color=col, size=11, symbol="star"),
                name=names[0],
                legendgroup=col,
                showlegend=showlegend),

            go.Scatter(x=[xs[1]], y=[ys[1]], mode="markers",
                marker=dict(color=col, size=9),
                name=names[1],
                legendgroup=col,
                showlegend=showlegend),
            ]
        return trcs

    def get_grey_heatmap(self):
        z = [[0,0], [0,0]]
        x = np.linspace(0,1,2)
        y = np.linspace(0,1,2)

        grey_map = go.Heatmap(
            x = x,
            y = y,
            z = z,
            colorscale = [[0, NULL_HEATMAP_COLOUR],[1, NULL_HEATMAP_COLOUR]],
            colorbar = invisible_colorbar(self.cbar_attrs['x'][0], self.cbar_attrs['y'][0]),
            )
        return grey_map



    def _get_ERFB_contour_single(self):
        z = EqualResFreqBreakdownArray(self.data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'dash')
        out['name'] = "Delta RFB"

        return out



    def _get_ESFY_contour_single(self):
        z = EqualSelectionArray(self.data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"

        return out





    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))
        

        fig.layout.coloraxis.showscale = False

        fig.update_layout(legend=dict(
                        x=0,
                        y=1.07,
                        xanchor="left", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))

        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", row=1, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", row=1, col=2, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="",                   row=1, col=2, range=[0,1], showgrid=False, zeroline=False, showticklabels=False)

        fig.update_xaxes(title="Dose sum",       row=2, col=2, range=[0,2.05], showgrid=False, showline=True)
        fig.update_yaxes(title="Effective life", row=2, col=2, showgrid=False)
        

        top_row = 1.08
        bottom_row = 0.51
        
        left = -0.04
        middle = 0.61

        cbarbox_x = 0.43
        cbarbox_y = 0.58


        annotz = [  
            dict(text="", x=self.cbar_attrs['x'][0]+0.014, 
                        y=self.cbar_attrs['intercepts'][0],
                        xref="paper", yref="paper",
                        ay=0, ax=30, arrowsize=2, arrowwidth=4,
                        arrowhead=0, arrowcolor="black"),
            dict(text="", x=self.cbar_attrs['x'][1]+0.014, 
                        y=self.cbar_attrs['intercepts'][1],
                        xref="paper", yref="paper",
                        ay=0, ax=30, arrowsize=2, arrowwidth=4,
                        arrowhead=0, arrowcolor="blue"),
            get_big_text_annotation(left, top_row, 'A'),
            get_big_text_annotation(middle, top_row, 'B'),
            get_big_text_annotation(left, bottom_row, 'C'),
            get_big_text_annotation(middle, bottom_row, 'D'),
            get_shape_annotation(cbarbox_x, cbarbox_y),
            ]

        fig.update_layout(annotations=annotz)


        return fig
    
    










class DoseSpaceScenariosPlot(BasicFig):
    def __init__(self, data_SS, data_SD, data_DS, data_DD, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 750

        self.data_SS = data_SS
        self.data_SD = data_SD
        self.data_DS = data_DS
        self.data_DD = data_DD

        self.cbar_attrs = dict(x=[0.43, 1.005], y=[0.21, 0.79], len=0.46)

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_scenarios_")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=2, cols=2, 
                    shared_xaxes=True,
                    shared_yaxes=True,
                    horizontal_spacing=0.15
                    )
        
        self.max_FY = max(
                np.amax(self.data_SS.FY),
                np.amax(self.data_SD.FY),
                np.amax(self.data_DS.FY),
                np.amax(self.data_DD.FY),
                )

        EL_SS_traces = self._get_DSS_EL_traces(self.data_SS, 0, 1, showlegend=True)
        EL_DS_traces = self._get_DSS_EL_traces(self.data_DS, 1, 1)
        EL_SD_traces = self._get_DSS_EL_traces(self.data_SD, 0, 0)
        EL_DD_traces = self._get_DSS_EL_traces(self.data_DD, 1, 0)

        fig.add_traces(EL_SS_traces, rows=1, cols=1)
        fig.add_traces(EL_DS_traces, rows=1, cols=2)
        
        fig.add_traces(EL_SD_traces, rows=2, cols=1)
        fig.add_traces(EL_DD_traces, rows=2, cols=2)

        fig = self._sort_layout(fig)

        return fig



    def _get_DSS_EL_traces(self, data, cbar_x, cbar_y, showlegend=False):  
        traces = []

        if showlegend:
            traces.append(self._get_DSS_ERFB_legend_entry())
            traces.append(self._get_DSS_ESFY_legend_entry())

        traces.append(self._get_DSS_FY_trace(data, cbar_x, cbar_y))
        
        traces.append(self._get_DSS_ERFB_contour_single(data))
        traces.append(self._get_DSS_ESFY_contour_single(data))

        return traces


    def _get_DSS_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    name=u"\u03A9<sub>SFY</sub>=0.5 contour"
                    )


    def _get_DSS_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="dash"),
                    name=u"\u0394<sub>RFB</sub>=0 contour"
                    )



    def _get_DSS_FY_trace(self, data, cbar_x, cbar_y):
        FYs = np.transpose(data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            colorscale = grey_colorscale_discrete(FYs),
            colorbar = my_colorbar_subplot("E.L.",
                            self.cbar_attrs['x'][cbar_x],
                            self.cbar_attrs['y'][cbar_y],
                            self.cbar_attrs['len']),
            )

        return heatmap



    def _get_DSS_ERFB_contour_single(self, data):
        z = EqualResFreqBreakdownArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'dash')
        out['name'] = "Delta RFB"
        out['coloraxis'] = "coloraxis"

        return out

    def _get_DSS_ESFY_contour_single(self, data):
        z = EqualSelectionArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"
        out['coloraxis'] = "coloraxis"

        return out





    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_layout(legend=dict(
                        x=0,
                        y=1.065,
                        xanchor="left", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))

        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", row=1, col=1, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", row=2, col=2, range=[0,1], showgrid=False, zeroline=False)
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)

        fig.layout.coloraxis.showscale = False
        
        top_row = 1.075
        bottom_row = 0.495
        
        left = 0.02
        middle = 0.59

        annotz = [
                  get_big_text_annotation(left, top_row, 'A'),
                  get_big_text_annotation(middle, top_row, 'B'),
                  get_big_text_annotation(left, bottom_row, 'C'),
                  get_big_text_annotation(middle, bottom_row, 'D'),
                    ]

        fig.update_layout(annotations=annotz)

        return fig
    
    




















class ParamScanPlotMeVsHobb(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 800

        self.use_pc = False

        self.data = self._process_data(data)
        
        self.data_high_low = self._process_data_high_low(data)

        fig = self._generate_figure()
        
        self.filename = f"../outputs/figures/paper_figs/par_scan_{conf_str}"

        self._save_and_show(fig)
    



    def _process_data(self, data):

        data_use = data.loc[:, ["run", 
                                "c_E_lowDoseMaxEL",
                                "c_R_maxContEL",
                                "c_E_maxContEL",
                                "I_R_best_value",
                                "I_E_best_value",
                                "RS",
                                "SR",
                                "RR",
                                "sr_prop",
                                "max_grid_EL",
                                "omega_1",
                                "omega_2",
                                 ]]
        

        print("\nNB NOW COMPARING LOW DOSES")

        if self.use_pc:
            # data_use['successMetric'] = 100*data_use['c_E_maxContEL']/data_use['c_R_maxContEL']
            data_use['successMetric'] = 100*data_use['c_E_lowDoseMaxEL']/data_use['c_R_maxContEL']

            bad_df = data_use[data_use['successMetric']>100]
            
            print(f"\n CURRENTLY REPLACING {bad_df.shape[0]} ROWS GREATER THAN 100 WITH 100... will need to check/justify!! \n")
            print(bad_df.loc[:, ["run", "c_R_maxContEL", "c_E_maxContEL", "c_E_lowDoseMaxEL", "max_grid_EL"]])

            # data_use.loc[data_use['successMetric']>100, ['successMetric']] = 100
            data_use.loc[data_use['successMetric']>100, ['successMetric']] = None
        
        else:
            # data_use['successMetric'] = data_use['c_R_maxContEL'] - data_use['c_E_maxContEL']
            data_use['successMetric'] = data_use['c_R_maxContEL'] - data_use['c_E_lowDoseMaxEL']
            

            bad_df = data_use[data_use['successMetric']<0]
            
            print(f"\n CURRENTLY REPLACING {bad_df.shape[0]} ROWS GREATER THAN 100 WITH 100... will need to justify!! \n")
            print(bad_df.loc[:, ["run", "c_R_maxContEL", "c_E_maxContEL", "c_E_lowDoseMaxEL", "max_grid_EL"]])

            # data_use.loc[data_use['successMetric']<0, ['successMetric']] = 0
            data_use.loc[data_use['successMetric']<0, ['successMetric']] = None



        data_use['IRFMetric'] = data_use.apply(lambda x: log10(x['RS']) - log10(x['SR']), axis=1)
        
        data_use['AsympMetric'] = data_use.apply(lambda x: x['omega_1']/(x['omega_1']+x['omega_2']), axis=1)
        
        # data_use['DecMetric'] = data_use.apply(lambda x: x['delta_1']/(x['delta_2']+x['delta_1']), axis=1)

        return data_use


    @staticmethod
    def _process_data_high_low(data):

        data_use = data.loc[:, [
                            # 'run', 
                            # 'c_R_minDS',
                            # 'c_R_maxDS', 
                            # 'max_dose_sums',
                            # 'min_dose_sums',
                            "RS",
                            "SR",
                            "RR", 
                            "sr_prop",
                            "c_R_lowDoseMaxEL",
                            "c_R_medDoseMaxEL",
                            "c_R_highDoseMaxEL",
                            "delta_1",
                            "delta_2",
                            "omega_1",
                            "omega_2",
                            "theta_1",
                            "theta_2",
                            ]]


        # srs = list(data_use["RS"])
        # rss = list(data_use["SR"])
        # data_use["singRatio"] = [logit10(rss[ii]) - logit10(srs[ii]) for ii in range(len(srs))]
        # data_use["RRlower"] = data_use["RR"] < data_use["SR"] * data_use["RS"]

        data_use['diff'] = data_use['c_R_highDoseMaxEL'] - data_use['c_R_lowDoseMaxEL']
        
        # data_use.sort_values(by="diff", ascending=False, inplace=True)
        
        # data_use = data_use.loc[(data_use.omega_1>0.8) & (data_use.omega_2>0.8) 
        #                     & (data_use.theta_1>8) & (data_use.theta_2>8)
        #                     # & (data_use.delta_1>1e-2) & (data_use.delta_2>1e-2)
        #                     ]

        # print("HIGH EFFICACY:")
        # print(data_use)


        
        data_use['highLow'] = data_use.apply(lambda x: x['c_R_highDoseMaxEL']/(x['c_R_lowDoseMaxEL'] + x['c_R_highDoseMaxEL']), axis=1)
        
        # print("High doses pref (all data):")
        # data['highLow'] = data.apply(lambda x: x['c_R_highDoseMaxEL']/(x['c_R_lowDoseMaxEL'] + x['c_R_highDoseMaxEL']), axis=1)
        # print(data.loc[data.highLow>0.5])


        return data_use










    def _generate_figure(self):

        fig = make_subplots(rows=2, cols=2, vertical_spacing=0.3, horizontal_spacing=0.16)

        trace_dict = self._get_trace_dict()

        fig.add_traces(trace_dict['IRFMetric'],   rows=1, cols=1)
        fig.add_traces(trace_dict['logRR'],       rows=1, cols=2)
        fig.add_traces(trace_dict['AsympMetric'], rows=2, cols=1)
        # fig.add_trace(trace_dict['DecMetric'], row=2, col=1)
        
        fig.add_traces(self._get_traces_high_low(), rows=2, cols=2)

        return self._sort_layout(fig)



    def _get_trace_dict(self):        
        out = {}

        data = self.data
        
        data['logRR'] = [log10(x) for x in data['RR']]        
        data['IRFMetric'] = [abs(x) for x in data['IRFMetric']]
        # DecMetric

        for key in ['IRFMetric', 'AsympMetric', 'logRR', 'sr_prop']:

            LIMIT = 0

            xblack = data[data['successMetric']<=LIMIT][key]
            yblack = data[data['successMetric']<=LIMIT]['successMetric']
            
            xblue = data[data['successMetric']>LIMIT][key]
            yblue = data[data['successMetric']>LIMIT]['successMetric']
            
            scatter_blue = go.Scatter(x=xblue,
                y=yblue,
                mode='markers',
                marker=dict(opacity=0.2, color="blue"))

            scatter_grey = go.Scatter(x=xblack,
                y=yblack,
                mode='markers',
                marker=dict(opacity=0.1, color="black"))

            out[key] = [scatter_blue, scatter_grey]

        return out



    def _get_traces_high_low(self):        
        data = self.data_high_low
        
        xblack = data[data['highLow']<0.5]['sr_prop']
        yblack = data[data['highLow']<0.5]['highLow']
        
        xblue = data[data['highLow']>=0.5]['sr_prop']
        yblue = data[data['highLow']>=0.5]['highLow']


        traces = [
                go.Scatter(x=xblue,
                            y=yblue,
                            mode='markers',
                            marker=dict(opacity=0.2, color="blue")),
                go.Scatter(x=xblack,
                            y=yblack,
                            mode='markers',
                            marker=dict(opacity=0.1, color="black"))
                    ]

        return traces




    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(False, self.width, self.height))
        fig = self._update_axes(fig)
        fig = self._add_corner_text_labels(fig)
        return fig


    def _update_axes(self, fig):

        fig.update_xaxes(title="Log difference in single resistance<br>frequencies (absolute value)", row=1, col=1, showgrid=False, showline=True)
        
        fig.update_xaxes(title="Double resistant<br>frequency (log scale)", row=1, col=2, showgrid=False, showline=True)

        fig.update_xaxes(title="Asymptote metric<br>  ", row=2, col=1, showgrid=False, showline=True)
        
        fig.update_xaxes(title="Proportion of between-season<br>sexual reproduction", row=2, col=2, showgrid=False, showline=True)

        
        if self.use_pc:
            ylab = "ESFY vs ERFB (max EL, %)"
        else:
            ylab = "ERFB EL - ESFY EL"

        fig.update_yaxes(title=ylab, row=1, col=1, showgrid=False, zeroline=False)
        fig.update_yaxes(            row=1, col=2, showgrid=False, showline=True, zeroline=False)
        fig.update_yaxes(title=ylab, row=2, col=1, showgrid=False, showline=True, zeroline=False)
        
        fig.update_yaxes(title="High/low dose metric", row=2, col=2, showgrid=False)

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
    
    










class SRPlot(BasicFig):
    def __init__(self, def_eff_data, low_eff_data, sf_ratio, grid_output, filestr) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 820

        self.def_eff_data = def_eff_data
        self.low_eff_data = low_eff_data

        self.sf_ratio = sf_ratio

        self.grid_output = grid_output

        fig = self._generate_figure()

        self.filename = f"../outputs/figures/paper_figs/sr_grid_{filestr}.png"

        self._save_and_show(fig)



    def _generate_figure(self):
        
        trcs1 = self.get_heatmap_and_contour(self.def_eff_data)
        trcs2 = self.get_heatmap_and_contour(self.low_eff_data)
        trcs3 = self.get_sfr_trace()
        trcs4 = self.get_SR_grid_EL_traces()
        
        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.2, vertical_spacing=0.2)
        
        fig.add_traces(trcs1, 1, 1)
        fig.add_traces(trcs2, 1, 2)
        fig.add_traces(trcs3, 2, 1)
        fig.add_traces(trcs4, 2, 2)
        
        self._sort_layout(fig)
        
        return fig

    def get_sfr_trace(self):
        data = self.sf_ratio
        traces = []

        # dashes = ["solid", "dash"]
        names = ["Low", "High"]
        cols = ["turquoise", "red"]

        for dd, name, col in zip(data.RR.unique(), names, cols):
            d1 = data.loc[data.RR==dd]
            
            xx = d1.d
            y1 = d1.sf_ratio_sing
            y2 = d1.sf_ratio_doub

            traces += [
                go.Scatter(x=xx, y=y2, name=f"{name} resistance: double resistant strain", 
                        legendgroup=name,
                        line=dict(color=col, dash="dot"), mode="lines"),
                go.Scatter(x=xx, y=y1, name=f"{name} resistance: single resistant strains", 
                        legendgroup=name,
                        line=dict(color=col, dash="solid"), mode="lines"),
                ]


        traces += [go.Scatter(x=[0.15,1], y=[1,1], 
                        line=dict(color="rgba(0,0,0,0.5)", dash="dash"), 
                        mode="lines",
                        name="No selection",
                        )]
        return traces




    def get_heatmap_and_contour(self, data):

        data = data.sort_values(['bs', 'ws'])
        
        x = np.array(data['bs'].unique())
        y = np.array(data['ws'].unique())

        z = np.array(data['Z_metric'])
        z = z.reshape((len(x), len(y)))
        z_transpose = np.transpose(z)

        heatmap = go.Heatmap(
            x = x,
            y = y,
            z = z_transpose,
            coloraxis = "coloraxis",
            )
        

        contour = contour_at_single_level(x, y, z_transpose, 0.5, 'darkred', 'dot')
        
        traces = [heatmap, contour]

        return traces













    def get_SR_grid_EL_traces(self):  
        traces = []

        traces.append(self._get_SR_grid_ERFB_legend_entry())
        traces.append(self._get_SR_grid_ESFY_legend_entry())

        traces.append(self._get_SR_grid_FY_trace())
        
        traces.append(self._get_SR_grid_ERFB_contour_single())
        traces.append(self._get_SR_grid_ESFY_contour_single())

        return traces


    def _get_SR_grid_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    legendgroup="cont",
                    name=u"\u03A9<sub>SFY</sub>=0.5 contour"
                    )


    def _get_SR_grid_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="dash"),
                    legendgroup="cont",
                    name=u"\u0394<sub>RFB</sub>=0 contour"
                    )



    def _get_SR_grid_FY_trace(self):
        data = self.grid_output
        FYs = np.transpose(data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            coloraxis="coloraxis2",
            )

        return heatmap



    def _get_SR_grid_ERFB_contour_single(self):
        data = self.grid_output
        z = EqualResFreqBreakdownArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'dash')
        out['name'] = "Delta RFB"

        return out

    def _get_SR_grid_ESFY_contour_single(self):
        data = self.grid_output
        z = EqualSelectionArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"

        return out












    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_layout(template="plotly_white")

        fig.update_layout(

            coloraxis=dict(colorbar=dict(
                    title = "EL<sub>full</sub> / (EL<sub>low</sub> + EL<sub>full</sub>)",
                    titleside = 'right',
                    len = 0.43,
                    y = 1.01,
                    yanchor="top",
                ),
                colorscale=divergent_color_scale(np.array(
                                        self.def_eff_data['Z_metric']), 0.5)
                ),

            coloraxis2=dict(colorbar=dict(
                    title = "Effective life",
                    titleside = 'right',
                    len = 0.43,
                    y = 0,
                    yanchor="bottom",
                ),
                colorscale=grey_colorscale_discrete(self.grid_output.FY)
                ),

                )

        #axes

        fig.update_xaxes(range=[0,1], row=1, col=1, showline=False, zeroline=False)
        fig.update_xaxes(range=[0,1], row=1, col=2, showline=False, zeroline=False)
        
        fig.update_xaxes(range=[0,1], row=2, col=1,
                title=dict(text="Dose",
                                font=dict(size=18, color="black")),            
                            showgrid=False)
        
        fig.update_xaxes(range=[0,1], row=2, col=2,
                title=dict(text="Dose (fungicide <i>A</i>)",
                                font=dict(size=18, color="black")),            
                            showgrid=False, zeroline=False)



        fig.update_yaxes(range=[0,1], row=1, col=1,
                title=dict(text="Within-season sexual<br>reproduction proportion (<i>p<sub>W</sub></i>)",
                                font=dict(size=18, color="black")),
                            showline=False, zeroline=False)

        fig.update_yaxes(range=[0,1], row=1, col=2,
                            showline=False, zeroline=False)
        
        fig.update_yaxes(range=[0,10], row=2, col=1,
                title=dict(text="Ratio of frequencies<br>year 2 vs year 1",
                                font=dict(size=18, color="black")),
                            showline=True, zeroline=False)
        
        fig.update_yaxes(range=[0,1], row=2, col=2,
                title=dict(text="Dose (fungicide <i>B</i>)",
                                font=dict(size=18, color="black")),            
                            showline=False, zeroline=False)
        

        # annotations
        top_row = 1.05
        bottom_row = 0.46
        
        left = 0
        middle = 0.6

        c1 = get_big_text_annotation(left, top_row, 'A: Default efficacy', xanchor="left")
        c2 = get_big_text_annotation(middle, top_row, 'B: Lower efficacy', xanchor="left")
        c3 = get_big_text_annotation(left, bottom_row, 'C', xanchor="left")
        c4 = get_big_text_annotation(middle, bottom_row, 'D', xanchor="left")

        x_lab = get_big_text_annotation(0.5, 0.55, 'Between-season sexual reproduction proportion (<i>p<sub>B</sub></i>)')

        txt_x = 1.03

        t1 = get_text_annotation(txt_x, 0.99, "High dose best", xanchor="left", yanchor="bottom")
        t2 = get_text_annotation(txt_x, 0.59, "Low dose best", xanchor="left", yanchor="top")
        
        c1['font']['size'] = 16
        c2['font']['size'] = 16
        c3['font']['size'] = 16
        c4['font']['size'] = 16
        
        t1['font']['color'] = LIGHT_GREY_TEXT
        t2['font']['color'] = LIGHT_GREY_TEXT

        x_lab['font'] = dict(size=18, color="black")
        
        annotz = [c1,
                    c2,
                    c3, 
                    c4,
                    x_lab,
                    t1,
                    t2,
                    ]

        fig.update_layout(annotations=annotz)

        fig.update_layout(legend=dict(
                            x=-0.05, y=-0.13,
                            orientation="h",
                            font=dict(size=14),
                            yanchor="top", 
                            xanchor="left")
                            )

        return fig






class DoseResponse(BasicFig):
    def __init__(self, dr_data, eff_data, filename) -> None:
        
        self.width = FULL_PAGE_WIDTH
        self.height = 700

        self.dr_data = dr_data
        self.eff_data = eff_data

        fig = self._generate_figure()

        self.filename = filename
        
        self._save_and_show(fig)



    def _generate_figure(self):
        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.25, shared_xaxes=True) # vertical_spacing=0.18,

        traces = self._get_dr_traces_1()
        fig.add_traces(traces, rows=1, cols=1)
        
        traces = self._get_dr_traces_2()
        fig.add_traces(traces, rows=2, cols=1)
        
        traces = self._get_conc_traces()
        fig.add_traces(traces, rows=1, cols=2)

        traces = self._get_eff_traces()
        fig.add_traces(traces, rows=2, cols=2)

        fig = self._sort_layout(fig)
        return fig
    


    def _get_eff_traces(self):
        xs = self.eff_data['x']
        ys = self.eff_data['y']
        names = self.eff_data['name']
        cols = self.eff_data['cols']
        dashes = self.eff_data['dashes']

        traces = []
 
        for x, y, name, col, dash in zip(xs, ys, names, cols, dashes):
            trc = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        showlegend=False,
                        name=name,
                        line=dict(color=col, dash=dash)
                        )
            traces.append(trc)
        
        traces.append(trc)

        return traces

    def _get_conc_traces(self):
        xs = self.eff_data['x']
        ys = self.eff_data['concs']
        names = self.eff_data['name']
        cols = self.eff_data['cols']
        dashes = self.eff_data['dashes']

        traces = []
 
        for x, y, name, col, dash in zip(xs, ys, names, cols, dashes):
            trc = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        showlegend=False,
                        name=name,
                        line=dict(color=col, dash=dash)
                        )
            traces.append(trc)
        
        traces.append(trc)

        return traces

    def _get_dr_traces_1(self):
        x = self.dr_data['x'][1]
        y = self.dr_data['y'][1]
        name = self.dr_data['name'][1]
        col = self.dr_data['cols'][1]
        dash = self.dr_data['dashes'][1]

        traces = []

        trc = go.Scatter(x=x,
                    y=y,
                    mode="lines",
                    name=name,
                    line=dict(color=col, dash=dash)
                    )

        traces.append(trc)


        return traces
    
    def _get_dr_traces_2(self):
        xs = self.dr_data['x']
        ys = self.dr_data['y']

        traces = []

        asymp_trc = go.Scatter(x=[xs[0][0], xs[0][-1]],
                        y=[0.2,0.2],
                        mode="lines",
                        showlegend=False,
                        line=dict(dash="dot", 
                                color="rgba(0,0,0,0.5)"),
                        )

        traces.append(asymp_trc)

        x = xs[0]
        y = ys[0]
        name = self.dr_data['name'][0]
        col = self.dr_data['cols'][0]
        dash = self.dr_data['dashes'][0]
        
        trc = go.Scatter(x=x,
                    y=y,
                    mode="lines",
                    name=name,
                    line=dict(color=col, dash=dash)
                    )

        traces.append(trc)


        return traces



    def _sort_layout(self, fig):

        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_xaxes(showgrid=False, row=1, col=1)
        fig.update_xaxes(title="Concentration (C)", range=[-0.01,1], showgrid=False, row=2, col=1)
        fig.update_xaxes(showgrid=False, row=1, col=2)
        fig.update_xaxes(title="Time (degree days)", showgrid=False, row=2, col=2)
        
        
        fig.update_yaxes(title=u"Growth rate<br>factor (\u03B4<sub>F,s</sub>)", range=[-0.01,1.01], showgrid=False, row=1, col=1)
        fig.update_yaxes(title=u"Growth rate<br>factor (\u03B4<sub>F,s</sub>)", range=[-0.01,1.01], showgrid=False, row=2, col=1)
        fig.update_yaxes(title="Concentration (C)", range=[-0.01,1.2], showgrid=False, showline=True, row=1, col=2)
        fig.update_yaxes(title=u"Growth rate<br>factor (\u03B4<sub>F,s</sub>)", range=[-0.01,1.01], showgrid=False, showline=True, row=2, col=2)
        
        
        fig.update_layout(legend=dict(
                        orientation="h",
                        x=0, y=1.1,
                        xanchor="left", 
                        yanchor="bottom",
                        font=dict(size=14),
                        ),
                        margin=dict(l=120)
                        )
        
        
        left = 0
        middle = 0.62
        
        top_row = 1.10
        bottom_row = 0.5

        annotz = [
            get_text_annotation(-0.035, 0.6, "Max. effect", "right", "bottom"),
            get_text_annotation(-0.035, 1.02, "No effect", "right", "top"),
            get_big_text_annotation(left, top_row, 'A', xanchor="left"),
            get_big_text_annotation(left,   bottom_row, 'B', xanchor="left"),
            get_big_text_annotation(middle, top_row, 'C', xanchor="left"),
            get_big_text_annotation(middle, bottom_row, 'D', xanchor="left"),
            ]

        fig.update_layout(annotations=annotz)

        return fig














class CurveDose(BasicFig):
    def __init__(self, dr_data, eff_data, filename) -> None:
        
        self.width = FULL_PAGE_WIDTH
        self.height = 450

        self.dr_data = dr_data
        self.eff_data = eff_data

        fig = self._generate_figure()

        self.filename = filename
        
        self._save_and_show(fig)



    def _generate_figure(self):
        fig = make_subplots(rows=1, cols=2, shared_yaxes=True)

        traces = self._get_dr_traces()
        fig.add_traces(traces, rows=1, cols=1)
        
        traces = self._get_eff_traces()
        fig.add_traces(traces, rows=1, cols=2)

        fig = self._sort_layout(fig)
        return fig
    


    def _get_eff_traces(self):
        xs = self.eff_data['x']
        ys = self.eff_data['y']
        names = self.eff_data['name']
        cols = self.eff_data['cols']
        dashes = self.eff_data['dashes']

        traces = []
 
        for x, y, name, col, dash in zip(xs, ys, names, cols, dashes):
            trc = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        showlegend=False,
                        name=name,
                        line=dict(color=col, dash=dash)
                        )
            traces.append(trc)
        
        traces.append(trc)

        return traces

    def _get_dr_traces(self):
        xs = self.dr_data['x']
        ys = self.dr_data['y']
        names = self.dr_data['name']
        cols = self.dr_data['cols']
        dashes = self.dr_data['dashes']

        traces = []

        for x, y, name, col, dash in zip(xs, ys, names, cols, dashes):
            trc = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        name=name,
                        line=dict(color=col, dash=dash)
                        )

            traces.append(trc)


        return traces



    def _sort_layout(self, fig):

        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_xaxes(title="Concentration (C)", range=[-0.01,2], showgrid=False, row=1, col=1)
        fig.update_yaxes(title=u"Growth rate<br>factor (\u03B4<sub>F,s</sub>)", range=[0,1], showgrid=False, row=1, col=1)
        
        fig.update_xaxes(title="Time (degree days)", showgrid=False, row=1, col=2)
        fig.update_yaxes(showgrid=False, row=1, col=2)
        
        
        fig.update_layout(legend=dict(
                        orientation="h",
                        x=0, y=1.12,
                        xanchor="left", 
                        yanchor="bottom",
                        font=dict(size=14),
                        ),
                        margin=dict(l=120)
                        )
        
        
        left = 0
        middle = 0.54
        top_row = 1.15

        annotz = [
            get_text_annotation(-0.05, 0.05, "Max. effect", "right", "bottom"),
            get_text_annotation(-0.05, 1.05, "No effect", "right", "top"),
            get_big_text_annotation(left, top_row, 'A', xanchor="left"),
            get_big_text_annotation(middle, top_row, 'B', xanchor="left"),
            ]

        fig.update_layout(annotations=annotz)

        return fig









class DeltaCurve(BasicFig):
    def __init__(self, data, filename) -> None:
        
        self.width = 0.7*FULL_PAGE_WIDTH
        self.height = 450

        self.data = data

        fig = self._generate_figure()

        self.filename = filename
        
        self._save_and_show(fig)

    def _generate_figure(self):
        traces = self._get_trace()
        fig = go.Figure(data=traces, layout=standard_layout(True, self.width, self.height))
        fig = self._sort_layout(fig)
        return fig
    
    def _get_trace(self):
        xs = self.data['x']
        ys = self.data['y']
        names = self.data['name']
        dashes = self.data['dash']

        traces = []

        for x,y,name,dash in zip(xs,ys,names,dashes):
            trc = go.Scatter(x=x,
                        y=y,
                        mode="lines",
                        name=name,
                        line=dict(dash=dash),
                        )
            traces.append(trc)
        # fix order
        traces.reverse()
        return traces



    def _sort_layout(self, fig):
        fig.update_xaxes(title=u"Fungicide effect (\u03B4)", showgrid=False)
        fig.update_yaxes(title=u"Selection metric", showgrid=False)
  
        fig.update_layout(legend=dict(
                        # orientation="h", 
                        x=1, y=1.05,
                        # font=dict(size=14),
                        xanchor="right", 
                        yanchor="top"),
                        )
        annotz = [
            get_text_annotation(0, -0.07, "Max. effect", "center", "top"),
            get_text_annotation(1, -0.07, "No effect", "center", "top"),
            ]

        fig.update_layout(annotations=annotz)
        return fig












class DoseSpaceScenarioSingle(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 750

        self.data = data

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_scenario_single")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = go.Figure()

        EL_SS_traces = self._get_DSS_EL_traces(self.data, showlegend=True)

        fig.add_traces(EL_SS_traces)

        fig = self._sort_layout(fig)

        return fig



    def _get_DSS_EL_traces(self, data, showlegend=False):  
        traces = []

        if showlegend:
            traces.append(self._get_DSS_ERFB_legend_entry())
            traces.append(self._get_DSS_ESFY_legend_entry())

        traces.append(self._get_DSS_FY_trace(data))
        
        traces.append(self._get_DSS_ERFB_contour_single(data))
        traces.append(self._get_DSS_ESFY_contour_single(data))

        return traces


    def _get_DSS_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    name=u"\u03A9<sub>SFY</sub>=0.5 contour"
                    )


    def _get_DSS_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="dash"),
                    name=u"\u0394<sub>RFB</sub>=0 contour"
                    )



    def _get_DSS_FY_trace(self, data):
        FYs = np.transpose(data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            colorscale = grey_colorscale_discrete(FYs),
            colorbar = my_colorbar("E.L."),
            )

        return heatmap



    def _get_DSS_ERFB_contour_single(self, data):
        z = EqualResFreqBreakdownArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'dash')
        out['name'] = "Delta RFB"

        return out

    def _get_DSS_ESFY_contour_single(self, data):
        z = EqualSelectionArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"

        return out





    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))
        fig.update_layout(template="plotly_white")

        fig.update_layout(legend=dict(
                        x=0,
                        y=1,
                        xanchor="left", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))

        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", range=[0,1], showgrid=False, zeroline=False)
        
        return fig
    
    




class DoseSpaceScenarioDouble(BasicFig):
    def __init__(self, data_def, data_pert, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 450

        self.data_def = data_def
        self.data_pert = data_pert

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_scenario_single")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=1, cols=2, shared_yaxes=True)

        EL_SS_traces_def = self._get_DSS_EL_traces(self.data_def, showlegend=True)
        EL_SS_traces_pert = self._get_DSS_EL_traces(self.data_pert, showlegend=False)

        fig.add_traces(EL_SS_traces_def, rows=1, cols=1)
        fig.add_traces(EL_SS_traces_pert, rows=1, cols=2)

        fig = self._sort_layout(fig)

        return fig



    def _get_DSS_EL_traces(self, data, showlegend=False):  
        traces = []

        if showlegend:
            traces.append(self._get_DSS_ERFB_legend_entry())

        traces.append(self._get_DSS_FY_trace(data))
        
        traces.append(self._get_DSS_ERFB_contour_single(data))

        return traces


    def _get_DSS_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    name=u"\u03A9<sub>SFY</sub>=0.5 contour"
                    )


    def _get_DSS_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="dash"),
                    name=u"\u0394<sub>RFB</sub>=0 contour"
                    )



    def _get_DSS_FY_trace(self, data):
        FYs = np.transpose(data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            coloraxis = "coloraxis",
            )

        return heatmap



    def _get_DSS_ERFB_contour_single(self, data):
        z = EqualResFreqBreakdownArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'dash')
        out['name'] = "Delta RFB"

        return out

    def _get_DSS_ESFY_contour_single(self, data):
        z = EqualSelectionArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_single_level(x, y, z_transpose, 0.5, 'blue', 'dot')
        out['name'] = "Equal Selection"

        return out



    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_layout(legend=dict(
                        x=1,
                        y=1,
                        xanchor="right", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))
        
        fig.update_layout(
            coloraxis = dict(colorbar=dict(
                    title = "E.L.",
                    titleside = 'right',
                ),
                colorscale = grey_colorscale_discrete(self.data_def.FY),
                )
            )
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", range=[0,1], showgrid=False, zeroline=False, row=1, col=1)
        fig.update_xaxes(title="Dose (fungicide <i>A</i>)", range=[0,1], showgrid=False, zeroline=False, row=1, col=2)
        fig.update_yaxes(title="Dose (fungicide <i>B</i>)", range=[0,1], showgrid=False, zeroline=False, row=1, col=1)

        
        left = 0
        middle = 0.54
        top_row = 1.15
        c1 = get_big_text_annotation(left, top_row, 'A', xanchor="left")
        c2 = get_big_text_annotation(middle, top_row, 'B', xanchor="left")

        annotz = [c1,c2]

        fig.update_layout(annotations=annotz)



        
        return fig
    
    



class SREffectAppendix(BasicFig):
    def __init__(self, out_l, out_m, out_h, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 900

        self.outputs_l = out_l
        self.outputs_m = out_m
        self.outputs_h = out_h

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/sr_app_")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=3, cols=3, 
                shared_xaxes=True,
                )

        cols3 = ["maroon", "cornflowerblue", "lime"]

        trcs_l_rr, trcs_l_rs, trcs_l_sr = self.get_rf_traces(self.outputs_l, [0, 0.2, 1], cols3, True)
        trcs_m_rr, trcs_m_rs, trcs_m_sr = self.get_rf_traces(self.outputs_m, [0, 0.2, 1], cols3)
        trcs_h_rr, trcs_h_rs, trcs_h_sr = self.get_rf_traces(self.outputs_h, [0, 0.2, 1], cols3)
        
        fig.add_traces(trcs_l_rr, rows=1, cols=1)
        fig.add_traces(trcs_m_rr, rows=1, cols=2)
        fig.add_traces(trcs_h_rr, rows=1, cols=3)
        
        fig.add_traces(trcs_l_rs, rows=2, cols=1)
        fig.add_traces(trcs_m_rs, rows=2, cols=2)
        fig.add_traces(trcs_h_rs, rows=2, cols=3)
        
        fig.add_traces(trcs_l_sr, rows=3, cols=1)
        fig.add_traces(trcs_m_sr, rows=3, cols=2)
        fig.add_traces(trcs_h_sr, rows=3, cols=3)

        fig = self._sort_layout(fig)

        return fig






    def get_rf_traces(self, outputs, bss, cols, showledge=False):
        
        traces_rr = []
        traces_rs = []
        traces_sr = []
        
        for data, bs, col in zip(outputs, bss, cols):
            traces = self.get_rf_trace(data, bs, col, showledge)

            traces_rr += [traces[0]]
            traces_rs += [traces[1]]
            traces_sr += [traces[2]]


        return traces_rr, traces_rs, traces_sr




    def get_rf_trace(self, data, bs, col, showledge):

        fy = int(np.amax(data.FY))

        fRR = data.start_freqs_DA["RR"][4, 19, :(fy+1)]
        fRS = data.start_freqs_DA["RS"][4, 19, :(fy+1)]
        fSR = data.start_freqs_DA["SR"][4, 19, :(fy+1)]
        # fSS = data.start_freqs_DA["SS"][-1,-1,:]

        # print(data.FY[5, 10:20])
        print(fy)

        x = list(range(len(fRR)))

        traces = []

        for ff, dash, name in zip([fRR, fRS, fSR],
                            ["solid", "dash", "dot"],
                            [f"Sexual proportion={bs}, strain <i>rr</i>",
                                f"Sexual proportion={bs}, strain <i>rs</i>", 
                                f"Sexual proportion={bs}, strain <i>sr</i>"],
                            ):

            ff = [None if ee==0 else log10(ee/(1-ee)) for ee in ff]

            showlegend = showledge if dash=="solid" else False

            name_use = name.split(",")[0]

            trc = dict(x=x, y=ff,
                        line=dict(dash=dash, color=col),
                        name=name_use,
                        # legendgroup=col,
                        showlegend=showlegend,
                        mode="lines"
                        )
            
            traces.append(trc)

        return traces






    def _sort_layout(self, fig):

        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_yaxes(title="Frequency <i>rr</i><br>(logit scale)", row=1, col=1)
        fig.update_yaxes(title="Frequency <i>rs</i><br>(logit scale)", row=2, col=1)
        fig.update_yaxes(title="Frequency <i>sr</i><br>(logit scale)", row=3, col=1)                
        
        fig.update_xaxes(range=[-0.5, 17.5], row=3, col=1)
        fig.update_xaxes(range=[-0.5, 17.5], row=3, col=2)
        fig.update_xaxes(range=[-0.5, 17.5], row=3, col=3)

        left = -0.015
        middle = 0.34
        right = 0.695

        top_row = 1.05
        middle_row = 0.68
        bottom_row = 0.32

        annotz = [
            get_big_text_annotation(left,   top_row, 'A: low <i>rr</i>', xanchor="left"),
            get_big_text_annotation(middle, top_row, 'B: medium <i>rr</i>', xanchor="left"),
            get_big_text_annotation(right,  top_row, 'C: high <i>rr</i>', xanchor="left"),
            
            get_big_text_annotation(left,   middle_row, 'D', xanchor="left"),
            get_big_text_annotation(middle, middle_row, 'E', xanchor="left"),
            get_big_text_annotation(right,  middle_row, 'F', xanchor="left"),

            get_big_text_annotation(left,   bottom_row, 'H', xanchor="left"),
            get_big_text_annotation(middle, bottom_row, 'I', xanchor="left"),
            get_big_text_annotation(right,  bottom_row, 'J', xanchor="left"),
            ]


        for cc in annotz:
            cc['font'] = dict(size=20, color=LIGHT_GREY_TEXT)
        
        
        x_lab2 = get_big_text_annotation(0.5, -0.08, 'Time (years)')
        x_lab2['font'] = dict(size=18, color="black")

        annotz += [x_lab2]

        fig.update_layout(annotations=annotz)

        fig.update_layout(legend=dict(x=-0.1, y=-0.072, 
                    # orientation="h",
                    font=dict(size=18),
                    xanchor="left", yanchor="top"))

        
        return fig
    
    





class SREffectResultsOld(BasicFig):
    def __init__(self, data, out_l, out_m, out_h, n_its, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 900

        self.data = data
        
        self.outputs_l = out_l
        self.outputs_m = out_m
        self.outputs_h = out_h


        self.n_its = n_its

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/sr_effect_old")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=3, cols=3, vertical_spacing=0.18, 
                row_heights=[0.5, 0.25, 0.25])

        cs = ["maroon", "turquoise", "lime"]

        step = 1

        cols2 = cs[0::(2*step)]
        cols3 = cs[0::step]

        print(cols2, cols3)
        
        fig.add_traces(self.get_EL_traces(1, 6), rows=1, cols=1)
        fig.add_traces(self.get_EL_traces(2, 10), rows=1, cols=2)
        fig.add_traces(self.get_EL_traces(3, 24), rows=1, cols=3)

        trcs_l_rr, trcs_l_sing = self.get_rf_traces(self.outputs_l, [0, 0.5, 1], cols3, True)
        trcs_m_rr, trcs_m_sing = self.get_rf_traces(self.outputs_m, [0, 1], cols2)
        trcs_h_rr, trcs_h_sing = self.get_rf_traces(self.outputs_h, [0, 1], cols2)
        
        fig.add_traces(trcs_l_rr, rows=2, cols=1)
        fig.add_traces(trcs_m_rr, rows=2, cols=2)
        fig.add_traces(trcs_h_rr, rows=2, cols=3)
        
        fig.add_traces(trcs_l_sing, rows=3, cols=1)
        fig.add_traces(trcs_m_sing, rows=3, cols=2)
        fig.add_traces(trcs_h_sing, rows=3, cols=3)

        fig = self._sort_layout(fig)

        return fig



    def get_EL_traces(self, col, index):

        df = copy.copy(self.data)

        # e.g. 0-9, 10-19, 20-29
        df = df.loc[((df["run"]<col*self.n_its) 
                    & (df["run"]>=(col-1)*self.n_its))]        


        traces = []

        for rr in df.run.unique():
            xx = df.loc[df["run"]==rr].bs_sex_prop
            yy = df.loc[df["run"]==rr].maxEL

            trc = dict(x=xx, 
                    y=yy,
                    showlegend=False,
                    opacity = 0.6,
                    name=rr,
                    )
            
            traces.append(trc)
        
        if col>1:
            index = index % ((col-1)*self.n_its)
        
        index = int(index)


        highlighted_trace = traces[index]
        highlighted_trace["opacity"] = 1
        highlighted_trace["line"] = dict(color="black", dash="dot")

        traces = traces[:(index)] + traces[(index+1):] + [highlighted_trace]

        return traces






    def get_rf_traces(self, outputs, bss, cols, showlegend=False):
        
        traces_rr = []
        traces_sing = []
        
        for data, bs, col in zip(outputs, bss, cols):
            traces = self.get_rf_trace(data, bs, col, showlegend)

            traces_rr += [traces[0]]
            traces_sing += traces[1:]


        return traces_rr, traces_sing
    

    def get_rf_trace(self, data, bs, col, showlegend):

        fy = int(np.amax(data.FY))
        
        fRR = data.start_freqs_DA["RR"][-1,-1,:(fy+1)]
        fRS = data.start_freqs_DA["RS"][-1,-1,:(fy+1)]
        fSR = data.start_freqs_DA["SR"][-1,-1,:(fy+1)]
        # fSS = data.start_freqs_DA["SS"][-1,-1,:]

        x = list(range(len(fRR)))

        traces = []

        for ff, dash, name in zip([fRR, fRS, fSR],
                            ["solid", "dash", "dot"],
                            [f"<i>rr</i>, sexual proportion={bs}", f"<i>rs</i>, sexual proportion={bs}", f"<i>sr</i>, sexual proportion={bs}"],
                            ):

            ff = [None if ee==0 else log10(ee/(1-ee)) for ee in ff]

            trc = dict(x=x, y=ff,
                        line=dict(dash=dash, color=col),
                        name=name,
                        legendgroup=name[:5],
                        showlegend=showlegend,
                        mode="lines"
                        )
            traces.append(trc)
        return traces






    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_yaxes(range=[5,22.5], title="Maximum<br>effective life", row=1, col=1)
        fig.update_yaxes(range=[5,22.5], row=1, col=2)
        fig.update_yaxes(range=[5,22.5], row=1, col=3)
        
        fig.update_yaxes(title="Double resistant<br>frequency<br>(logit scale)", row=2, col=1)
        fig.update_yaxes(row=2, col=2)
        fig.update_yaxes(row=2, col=3)
        
        fig.update_yaxes(title="Single resistant<br>frequencies<br>(logit scale)", row=3, col=1)
        fig.update_yaxes(row=3, col=2)
        fig.update_yaxes(row=3, col=3)


        fig.update_xaxes(range=[-0.5,16.5], row=2, col=1)
        fig.update_xaxes(range=[-0.5,16.5], row=3, col=1)
        
        fig.update_xaxes(range=[-0.5,8.5],  row=2, col=2)
        fig.update_xaxes(range=[-0.5,8.5],  row=3, col=2)

        fig.update_xaxes(range=[-0.5,13.5], row=2, col=3)
        fig.update_xaxes(range=[-0.5,13.5], row=3, col=3)
        

        left = -0.015
        middle = 0.34
        right = 0.695

        top_row = 1.05
        middle_row = 0.56
        bottom_row = 0.23

        annotz = [
            get_big_text_annotation(left,   top_row, 'A: low <i>rr</i>', xanchor="left"),
            get_big_text_annotation(middle, top_row, 'B: medium <i>rr</i>', xanchor="left"),
            get_big_text_annotation(right,  top_row, 'C: high <i>rr</i>', xanchor="left"),
            
            get_big_text_annotation(left,   middle_row, 'D', xanchor="left"),
            get_big_text_annotation(middle, middle_row, 'E', xanchor="left"),
            get_big_text_annotation(right,  middle_row, 'F', xanchor="left"),

            get_big_text_annotation(left,   bottom_row, 'H', xanchor="left"),
            get_big_text_annotation(middle, bottom_row, 'I', xanchor="left"),
            get_big_text_annotation(right,  bottom_row, 'J', xanchor="left"),
            ]


        for cc in annotz:
            cc['font'] = dict(size=20, color=LIGHT_GREY_TEXT)
        
        x_lab = get_big_text_annotation(0.5, 0.61, 'Between-season sexual reproduction proportion (<i>p<sub>B</sub></i>)')
        x_lab['font'] = dict(size=18, color="black")
        
        x_lab2 = get_big_text_annotation(0.5, -0.08, 'Time (years)')
        x_lab2['font'] = dict(size=18, color="black")

        annotz += [x_lab, x_lab2]

        fig.update_layout(annotations=annotz)

        fig.update_layout(legend=dict(x=-0.15, y=-0.15, orientation="h",
                    xanchor="left", yanchor="top",
                    font=dict(size=14),
                    ))

        
        return fig







class SREffectResults3Panel(BasicFig):
    def __init__(self, data, double_freqs, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 450

        self.data = data

        self.double_freqs = double_freqs

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/sr_effect")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=1, cols=3, shared_yaxes=True)

        fig.add_traces(self.get_traces(1), rows=1, cols=1)
        fig.add_traces(self.get_traces(2), rows=1, cols=2)
        fig.add_traces(self.get_traces(3), rows=1, cols=3)

        fig = self._sort_layout(fig)

        return fig



    def get_traces(self, col):

        df = copy.copy(self.data)

        double_freq_factor = self.double_freqs[col-1]

        df = df.loc[df["double_freq_factor"]==double_freq_factor]

        traces = []

        colors = list(px.colors.qualitative.Bold)


        for rr in df.run.unique():
            xx = df.loc[df["run"]==rr].bs_sex_prop
            yy = df.loc[df["run"]==rr].maxEL

            showlegend = True if col==1 else False

            colr = colors[int(rr)]
            
            opacity = 0.8 if rr in [0,1,2] else 0.15

            trc = dict(x=xx, 
                    y=yy,
                    name=f"Scenario {int(rr+1)}",
                    showlegend=showlegend,
                    line=dict(color=colr),
                    opacity=opacity,
                    )

            traces.append(trc)

        return traces


    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))

        fig.update_yaxes(range=[6.5,17.5], title="Maximum effective life", row=1, col=1)

        left = -0.02
        middle = 0.33
        right = 0.69

        top_row = 1.12

        c1 = get_big_text_annotation(left,   top_row, 'A: low <i>rr</i>', xanchor="left")
        c2 = get_big_text_annotation(middle, top_row, 'B: expected <i>rr</i>', xanchor="left")
        c3 = get_big_text_annotation(right,  top_row, 'C: high <i>rr</i>', xanchor="left")
        c1['font'] = dict(size=18, color=LIGHT_GREY_TEXT)
        c2['font'] = dict(size=18, color=LIGHT_GREY_TEXT)
        c3['font'] = dict(size=18, color=LIGHT_GREY_TEXT)

        x_lab = get_big_text_annotation(0.5, -0.14, 'Between-season sexual reproduction proportion (<i>p<sub>B</sub></i>)')
        x_lab['font'] = dict(size=18, color="black")

        annotz = [c1,c2,c3,x_lab]

        fig.update_layout(annotations=annotz)
        fig.update_layout(legend=dict(x=-0.05, y=1.15, 
            font=dict(size=16),
            yanchor="bottom",
            xanchor="left",
            orientation="h",
            ))

        return fig

"""Figures for the paper."""

from __future__ import annotations
from abc import ABC, abstractmethod
from model.utils import logit10
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
# import pandas as pd
from math import log10
from scipy import stats
from PIL import Image


from model.strategy_arrays import EqualResFreqBreakdownArray, EqualSelectionArray

from plotting.traces import contour_at_0, contour_at_single_level

from plotting.utils import get_text_annotation, get_arrow_annotation, grey_colorscale_discrete, grey_colorscale_discrete_N, invisible_colorbar, my_colorbar_subplot, standard_layout, \
    grey_colorscale, my_colorbar, get_big_text_annotation

from plotting.consts import ATTRS_DICT, LABEL_COLOR, NULL_HEATMAP_COLOUR, TITLE_MAP, PLOT_WIDTH, PLOT_HEIGHT, \
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

        for key in ['fung_1', 'fung_2']:
            out.append(self.get_DPC_trace(key))

        return out



    def get_DPC_trace(self, key):
        clr = ATTRS_DICT[key]['color']
        dash = ATTRS_DICT[key]['dash']
        name = ATTRS_DICT[key]['name']

        return go.Scatter(x=self.xx,
                    y=vars(self.states_list[self.DPC_year])[key],
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

        yy = self.data.yield_vec
        
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
        
        y1 = self.data.res_vec_dict['f1']
        y2 = self.data.res_vec_dict['f2']
        
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
        
        img = Image.open("create_figs/img/diagram.png")

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
    


# End of CombinedModelPlot














class DoseSpaceOverview(BasicFig):
    def __init__(self, data, contour_data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 720

        self.data = data

        self.contour_data = contour_data

        self.cbar_attrs = dict(x=[0.375, 1.005], y=[0.21, 0.79], len=0.46)

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_overview_")

        self._save_and_show(fig)



    def _generate_figure(self):

        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.25)

        EL_traces = self._get_EL_traces()
        
        ESFY_traces = self._get_ESFY_traces()
        ERFB_traces = self._get_ERFB_traces()
        cont_traces = self._get_contour_traces()

        fig.add_traces(EL_traces, rows=1, cols=1)
        fig.add_traces(ESFY_traces, rows=1, cols=2)
        fig.add_traces(ERFB_traces, rows=2, cols=1)
        fig.add_traces(cont_traces, rows=2, cols=2)

        fig = self._sort_layout(fig)

        return fig



    def _get_EL_traces(self):        
        traces = []

        traces.append(self._get_ERFB_legend_entry())
        traces.append(self._get_ESFY_legend_entry())

        traces.append(self._get_FY_trace())
        
        traces.append(self._get_ERFB_contour_single())
        traces.append(self._get_ESFY_contour_single())

        return traces

    def _get_ERFB_traces(self):
        traces = []
        traces.append(self.get_grey_heatmap())
        traces.append(self.get_ERFB_trace())
        traces.append(self._get_ERFB_contour_single())
        return traces    
    
    def _get_ESFY_traces(self):
        traces = []
        traces.append(self.get_grey_heatmap())
        traces.append(self.get_ESFY_trace())
        traces.append(self._get_ESFY_contour_single())
        return traces
    
    def _get_contour_traces(self):
        traces = []
        traces.append(self.get_ERFB_cont_trace())
        traces.append(self.get_ESFY_cont_trace())
        return traces

    def get_ERFB_cont_trace(self):
        x = self.contour_data['ERFB']['DS']
        y = self.contour_data['ERFB']['FY']
        return go.Scatter(x=x, y=y, mode="lines", line=dict(color="black"), showlegend=False)

    def get_ESFY_cont_trace(self):
        x = self.contour_data['ESFY']['DS']
        y = self.contour_data['ESFY']['FY']
        return go.Scatter(x=x, y=y, mode="lines", line=dict(color="blue", dash="dot"), showlegend=False)






    def _get_ESFY_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="blue", dash="dot"),
                    name=u"\u0394<sub>SFY</sub>=0.5 contour"
                    )


    def _get_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="solid"),
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
            # colorscale = grey_colorscale_discrete(FYs),
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
            # colorscale = grey_colorscale_discrete(FYs),
            colorbar = my_colorbar_subplot(u"\u0394<sub>SFY", self.cbar_attrs['x'][1], self.cbar_attrs['y'][1], self.cbar_attrs['len'])
            )

        return heatmap


    def get_grey_heatmap(self):
        z = [[0,0], [0,0]]
        x = np.linspace(0,1,2)
        y = np.linspace(0,1,2)

        grey_map = go.Heatmap(
            x = x,
            y = y,
            z = z,
            colorscale = [[0, NULL_HEATMAP_COLOUR],[1, NULL_HEATMAP_COLOUR]],
            colorbar = invisible_colorbar(self.cbar_attrs['x'][0], self.cbar_attrs['y'][0])
            )
        return grey_map

    def _get_ERFB_contour_single(self):
        z = EqualResFreqBreakdownArray(self.data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'solid')
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

        fig.update_layout(legend=dict(
                        x=0,
                        y=1.05,
                        xanchor="left", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))

        fig.update_yaxes(title="Dose (fungicide B)", row=1, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide A)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide B)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide A)", row=1, col=2, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="",                   row=1, col=2, range=[0,1], showgrid=False, zeroline=False, showticklabels=False)

        fig.update_xaxes(title="Dose sum",       row=2, col=2, range=[0.3,2], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Effective life", row=2, col=2, showgrid=False, zeroline=False)
        

        top_row = 1.07
        bottom_row = 0.49
        
        left = -0.04
        middle = 0.59


        annotz = [dict(text="", x=self.cbar_attrs['x'][0]+0.014, y=0.24, xref="paper", yref="paper",
                        ay=0, ax=30, arrowsize=2, arrowwidth=4, arrowhead=0, arrowcolor="black"),
                  dict(text="", x=self.cbar_attrs['x'][1]+0.014, y=0.77, xref="paper", yref="paper",
                        ay=0, ax=30, arrowsize=2, arrowwidth=4, arrowhead=0, arrowcolor="blue"),
                  get_big_text_annotation(left, top_row, 'A'),
                  get_big_text_annotation(middle, top_row, 'B'),
                  get_big_text_annotation(left, bottom_row, 'C'),
                  get_big_text_annotation(middle, bottom_row, 'D'),
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
                    name=u"\u0394<sub>SFY</sub>=0.5 contour"
                    )


    def _get_DSS_ERFB_legend_entry(self):
        return go.Scatter(x=[1], 
                    y=[1],
                    mode="lines",
                    line=dict(color="black", dash="solid"),
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
            colorbar = my_colorbar_subplot("EL",
                            self.cbar_attrs['x'][cbar_x],
                            self.cbar_attrs['y'][cbar_y],
                            self.cbar_attrs['len']),
            # coloraxis = "coloraxis",
            )

        return heatmap



    def _get_DSS_ERFB_contour_single(self, data):
        z = EqualResFreqBreakdownArray(data).array

        x = np.linspace(0, 1, z.shape[0])
        y = np.linspace(0, 1, z.shape[1])

        z_transpose = np.transpose(z)

        out = contour_at_0(x, y, z_transpose, 'black', 'solid')
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
                        x=0,
                        y=1.065,
                        xanchor="left", 
                        yanchor="bottom",
                        
                        orientation="h",
                        font=dict(
                            color=LABEL_COLOR
                        )
                        ))

        fig.update_yaxes(title="Dose (fungicide B)", row=1, col=1, range=[0,1], showgrid=False, zeroline=False)
        fig.update_yaxes(title="Dose (fungicide B)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        
        fig.update_xaxes(title="Dose (fungicide A)", row=2, col=2, range=[0,1], showgrid=False, zeroline=False)
        fig.update_xaxes(title="Dose (fungicide A)", row=2, col=1, range=[0,1], showgrid=False, zeroline=False)
        
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
    
    










class DosesScatterPlot(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 620

        self.data = data

        fig = self._generate_figure()
        
        self.filename = conf_str.replace("/grid/", "/paper_figs/doses_scatter_")

        self._save_and_show(fig)



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
        FYs = self.data.FY

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
        fig.update_xaxes(title=r"$\Delta_{RFB}$")
        fig.update_yaxes(title="Effective life")
        return fig
    












class ParamScanPlotMeVsHobb(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 800

        # data.set_index("run", inplace=True)
        # par_data.set_index("run", inplace=True, drop=False)

        # data = pd.concat([data, par_data], axis=1)

        self.data = self._process_data(data)
        
        self.data_high_low = self._process_data_high_low(data)

        fig = self._generate_figure()
        
        self.filename = f"../outputs/figures/paper_figs/par_scan_{conf_str}"

        self._save_and_show(fig)
    



    def _process_data(self, data):

        data_use = data.loc[:, ["run", 
                                "c_R_maxContEL",
                                "c_E_maxContEL",
                                "I_R_best_value",
                                "I_E_best_value",
                                "RS",
                                "SR",
                                "RR",
                                "sr_prop",
                                "max_grid_EL",
                                # "delta_1",
                                # "delta_2",
                                "omega_1",
                                "omega_2",
                                 ]]
        

        data_use['successMetric'] = self._get_success_metric(data_use)
        
        bad_df = data_use[data_use['successMetric']>100]
        
        print(f"\n CURRENTLY FILTERING OUT {bad_df.shape[0]} ROWS GREATER THAN 100... will need to justify!! \n")
        print(bad_df.loc[:, ["run", "c_R_maxContEL", "c_E_maxContEL", "max_grid_EL"]])

        data_use.loc[data_use['successMetric']>100, ['successMetric']] = None

        data_use['IRFMetric'] = data_use.apply(lambda x: logit10(x['RS']) - logit10(x['SR']), axis=1)
        
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
                            "omega_1",
                            "omega_2",
                            ]]

        srs = list(data_use["RS"])
        rss = list(data_use["SR"])

        data_use["singRatio"] = [logit10(rss[ii]) - logit10(srs[ii]) for ii in range(len(srs))]
        
        data_use["RRlower"] = data_use["RR"] < data_use["SR"] * data_use["RS"]
        

        data_use = data_use.loc[(data_use.omega_1>0.8) & (data_use.omega_2>0.8)]
        # data_use = data_use.loc[(data_use.omega_1>0.8) & (data_use.omega_2>0.8) & (data_use.sr_prop>0.5)]
        
        # data_use = data_use.loc[(data_use.omega_1>0.92) & (data_use.omega_2>0.92) & (data_use.sr_prop>0.6)]
        # data_use = data_use.loc[(data_use.omega_1>0.8) & (data_use.omega_2>0.8) & (data_use.RRlower)]

        
        # data_use['highLow'] = data_use.apply(lambda x: x['c_R_highDoseMaxEL']-x['c_R_lowDoseMaxEL'],
        #                                                         axis=1)
        data_use['highLow'] = data_use.apply(lambda x: x['c_R_highDoseMaxEL']/(x['c_R_lowDoseMaxEL'] + x['c_R_highDoseMaxEL']),
                                                                axis=1)
        
        
        print("High doses pref:")
        print(data_use.loc[data_use.highLow>0.5])
        # print(data_use.loc[data_use.highLow>=0.5])
        # print(data_use.loc[data_use.highLow<0.5])
        # print(data_use.loc[data_use.highLow==1].describe())
        # exit()

        # print(sum(data_use['highLow']))

        return data_use



    def _get_traces_high_low(self):        
        data = self.data_high_low
        
        y = data['highLow']
        # y = data['maxMixStrength']
        # y = data['meanMixStrength']
        
        x4 = data['sr_prop']


        traces = [go.Scatter(x=x4,
                            y=y,
                            mode='markers',
                            marker=dict(opacity=0.2, color="black"))]

        return traces






    @staticmethod
    def _get_success_metric(x):

        # bv_ES = np.asarray(x['I_E_best_value'])
        # bv_RFB = np.asarray(x['I_R_best_value'])

        c_E_EL = np.asarray(x['c_E_maxContEL'])
        c_R_EL = np.asarray(x['c_R_maxContEL'])

        # out = [100*max(bv_ES[i], c_E_EL[i]) / max(bv_RFB[i], c_R_EL[i]) for i in range(len(bv_ES))]
        out = [100*c_E_EL[i] / c_R_EL[i] for i in range(len(c_R_EL))]

        return out



        


    def _generate_figure(self):
        

        ugly_fig = self._add_traces_to_figure()

        fig = self._sort_layout(ugly_fig)

        return fig



    def _get_traces(self):        
        out = {}

        data = self.data
        
        data['logRR'] = [log10(x) for x in data['RR']]        
        data['IRFMetric'] = [abs(x) for x in data['IRFMetric']]
        # DecMetric

        for key in ['IRFMetric', 'AsympMetric', 'logRR', 'sr_prop']:
            
            scatter = go.Scatter(x=data[key],
                y=data['successMetric'],
                mode='markers',
                marker=dict(opacity=0.2, color="blue"))

            out[key] = scatter

        return out




    def _add_traces_to_figure(self):
        fig = make_subplots(rows=2, cols=2, vertical_spacing=0.3, horizontal_spacing=0.2)

        trace_dict = self._get_traces()

        fig.add_trace(trace_dict['IRFMetric'], row=1, col=1)
        fig.add_trace(trace_dict['logRR'], row=1, col=2)
        fig.add_trace(trace_dict['AsympMetric'], row=2, col=1)
        # fig.add_trace(trace_dict['DecMetric'], row=2, col=1)
        
        fig.add_traces(self._get_traces_high_low(), rows=2, cols=2)

        return fig



    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(False, self.width, self.height))
        fig = self._update_axes(fig)
        fig = self._add_corner_text_labels(fig)
        return fig


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title="Log diff. in single res.<br>freqs. (absolute value)", row=1, col=1)
        
        fig.update_xaxes(title="Double resistant<br>frequency (log scale)", row=1, col=2)

        fig.update_xaxes(title="Asymptote metric", row=2, col=1)
        
        fig.update_xaxes(title="Proportion of between-season<br>sexual reproduction", row=2, col=2)

        
        fig.update_yaxes(title="ESFY vs ERFB (max EL, %)", row=1, col=1)

        fig.update_yaxes(title="ESFY vs ERFB (max EL, %)", row=2, col=1)
        fig.update_yaxes(title="High/low dose metric<br>(high efficacy fungicides)", row=2, col=2)

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
    
    










class SR_grid(BasicFig):
    def __init__(self, def_eff_data, low_eff_data, sf_ratio, filestr) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 720

        self.def_eff_data = def_eff_data
        self.low_eff_data = low_eff_data

        self.sf_ratio = sf_ratio

        fig = self._generate_figure()

        self.filename = f"../outputs/figures/paper_figs/sr_grid_{filestr}.png"

        self._save_and_show(fig)



    def _generate_figure(self):
        
        trcs1 = self.get_trace(self.def_eff_data)
        trcs2 = self.get_trace(self.low_eff_data)
        trcs3 = self.get_sfr_trace()
        
        fig = make_subplots(rows=2, cols=2, horizontal_spacing=0.1, vertical_spacing=0.18)
        
        fig.add_traces(trcs1, 1, 1)
        fig.add_traces(trcs2, 1, 2)
        fig.add_traces(trcs3, 2, 1)
        
        self._sort_layout(fig)
        
        return fig

    def get_sfr_trace(self):
        data = self.sf_ratio
        traces = []

        dashes = ["solid", "dash"]
        # cols1 = ["blue", "blue"]
        # cols2 = ["blue", "blue"]
        names = ["Low", "High"]

        for dd, dash, name in zip(data.RR.unique(), dashes, names):
            d1 = data.loc[data.RR==dd]
            
            xx = d1.d
            y1 = d1.sf_ratio_sing
            y2 = d1.sf_ratio_doub

            traces += [
                go.Scatter(x=xx, y=y2, name=f"{name} initial freq.: double resistant", 
                        legendgroup=name,
                        line=dict(color="black", dash=dash), mode="lines"),
                go.Scatter(x=xx, y=y1, name=f"{name} initial freq.: single resistant", 
                        legendgroup=name,
                        line=dict(color="blue", dash=dash), mode="lines"),
                ]

        # reverse order for legend
        # traces.reverse()

        traces += [go.Scatter(x=[0.15,1], y=[1,1], 
                        line=dict(color="rgba(0,0,0,0.5)", dash="dot"), 
                        # showlegend=False,
                        mode="lines",
                        name="No selection",
                        )]
        return traces

    def get_trace(self, data):

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
        

        contour = contour_at_single_level(x, y, z_transpose, 0.5, 'red', 'solid')
        
        traces = [heatmap, contour]

        return traces



    def _sort_layout(self, fig):
        fig.update_layout(standard_layout(True, self.width, self.height))
        fig.update_layout(
            coloraxis=dict(colorbar=dict(
                    title = "EL<sub>full</sub> / (EL<sub>low</sub> + EL<sub>full</sub>)",
                    titleside = 'right',
                    len = 0.43,
                    y = 1.01,
                    yanchor="top",
                )))

        fig.update_xaxes(range=[0,1], row=1, col=1)
        fig.update_xaxes(range=[0,1], row=1, col=2)
        
        fig.update_xaxes(title="Dose", range=[0.15,1], row=2, col=1, showgrid=False)

        fig.update_yaxes(range=[0,1], row=1, col=1,
                title=dict(text="Within-season<br>sexual reproduction proportion",
                            font=dict(size=18, color="black")),
                            )
        
        fig.update_yaxes(title="Ratio of frequencies<br>year 2 vs year 1",
                            range=[0,9], row=2, col=1)

        
        top_row = 1.05
        bottom_row = 0.45
        
        left = 0
        middle = 0.55

        c1 = get_big_text_annotation(left, top_row, 'A: Default', xanchor="left")
        c2 = get_big_text_annotation(middle, top_row, 'B: Lower efficacy', xanchor="left")
        c3 = get_big_text_annotation(left, bottom_row, 'C: ratio of initial frequencies', xanchor="left")
        # c4 = get_big_text_annotation(middle, bottom_row, 'D', xanchor="left")

        text = get_big_text_annotation(0.5, 0.56, 'Between-season sexual reproduction proportion')
        
        c1['font']['size'] = 16
        c2['font']['size'] = 16
        c3['font']['size'] = 16
        # c4['font']['size'] = 16
        text['font'] = dict(size=18, color="black")
        
        annotz = [c1,
                    c2,
                    c3, 
                    # c4,
                    text]

        fig.update_layout(annotations=annotz)
        fig.update_layout(legend=dict(x=0.5, y=0, yanchor="bottom", xanchor="left"))
        return fig




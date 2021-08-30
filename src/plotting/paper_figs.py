"""Figures for the paper."""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
# import pandas as pd
from math import log10
from scipy import stats
from PIL import Image


from model.strategy_arrays import EqualResFreqBreakdownArray, EqualSelectionArray

from plotting.traces import contour_at_0, contour_at_single_level

from plotting.utils import get_text_annotation, get_arrow_annotation, standard_layout, \
    grey_colorscale, my_colorbar, get_big_text_annotation

from plotting.consts import ATTRS_DICT, TITLE_MAP, PLOT_WIDTH, PLOT_HEIGHT, \
        FULL_PAGE_WIDTH


# Paper Figs

class BasicFig:
    def __init__(self) -> None:
        pass

    def _save_and_show(self, fig):
        fig.show()
        
        print(f"""
        saving figure to:
        {self.filename}
        """)

        fig.write_image(self.filename)





class DiseaseProgressCurvesAll(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 650

        self.xx = data.states_list[0].t
        
        self.states_list = data.states_list

        fig = self._generate_figure()

        self.filename = conf_str.replace("/single/", "/paper_figs/model_overview_")

        self._save_and_show(fig)




    def _generate_figure(self):

        traces_dict = self.get_model_output_overview_traces()

        ugly_fig = self.add_traces_to_layout_model_output_overview(traces_dict)

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
        # ind = ATTRS_DICT[key]['ind']

        yy = vars(self.states_list[0])[key]

        return go.Scatter(x=self.xx,
                    y=yy,
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
    












class CombinedModelPlot:
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 950
        
        self.states_list = data.states_list

        self.data = data

        self.DPC_year = 5

        self.xx = data.states_list[self.DPC_year].t

        fig = self._generate_figure()

        self._save_and_show(fig, conf_str)




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
    


    def _save_and_show(self, fig, conf_str):
        fig.show()
        filename = conf_str.replace("/single/", "/paper_figs/model_overview_combined_")
        
        print("saving figure to: \n", filename)

        fig.write_image(filename)


# End of CombinedModelPlot














class DoseSpaceScenariosPlot(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 650

        self.data = data

        fig = self._generate_figure()

        self.filename = conf_str.replace("/grid/", "/paper_figs/dose_space_")

        self._save_and_show(fig)



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
        FYs = np.transpose(self.data.FY)

        xheat = np.linspace(0, 1, FYs.shape[0])
        yheat = np.linspace(0, 1, FYs.shape[1])

        heatmap = go.Heatmap(
            x = xheat,
            y = yheat,
            z = FYs,
            colorscale = grey_colorscale(FYs),
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
        FYs = np.transpose(self.data.FY)

        x = np.linspace(0, 1, FYs.shape[0])
        
        minEqDoseELVec = np.asarray([float(FYs[i, i]) for i in range(FYs.shape[0])])

        ind = np.where(minEqDoseELVec>1)

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
        fig = self._update_axes(fig)
        return fig
    


    @staticmethod
    def _update_axes(fig):
        fig.update_xaxes(title=r"$\Delta_{RFB}$")
        fig.update_yaxes(title="Effective life")
        return fig






class YieldAndRfPlot(BasicFig):
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

        yy = self.data.yield_vec
        
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
        
        y1 = self.data.res_vec_dict['f1']
        y2 = self.data.res_vec_dict['f2']
        
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







class ParamScanPlotMeVsHobb(BasicFig):
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






class ParamScanPlotHighLowDose(BasicFig):
    def __init__(self, data, conf_str) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 800

        self.data = self._process_data(data)

        fig = self._generate_figure()

        self.filename = f"../outputs/figures/paper_figs/param_scan_high_low_{conf_str}"

        self._save_and_show(fig)
    



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











class SR_grid(BasicFig):
    def __init__(self, data) -> None:

        self.width = FULL_PAGE_WIDTH

        self.height = 780

        self.data = data

        fig = self._generate_figure()

        self.filename = "../outputs/figures/paper_figs/sr_grid.png"

        self._save_and_show(fig)



    def _generate_figure(self):
        
        trc = self.get_trace()
        
        fig = go.Figure(trc)
        
        self._update_layout(fig)
        
        return fig



    def get_trace(self):
        data = self.data
        print(self.data)

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
            # colorscale = grey_colorscale(FYs),
            # colorbar = my_colorbar("$EL_l$ / ($EL_l$ + $EL_f$)")
            colorbar = my_colorbar("EL<sub>full</sub> / (EL<sub>low</sub> + EL<sub>full</sub>)")
            )
        

        contour = contour_at_single_level(x, y, z_transpose, 0.5, 'red', 'solid')
        
        traces = [heatmap, contour]

        return traces



    def _update_layout(self, fig):
        fig.update_layout(standard_layout(False, self.width, self.height))
        fig.update_xaxes(title="Between season sexual reproduction proportion", range=[0,1])
        fig.update_yaxes(title="Within season sexual reproduction proportion", range=[0,1])

        return fig




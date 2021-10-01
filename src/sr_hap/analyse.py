from math import log10
from model.simulator import RunGrid
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


from model.config_classes import GridConfig
from sr_hap.scan import get_sr_scan_params



def check_run(n_doses, bs, dfp):

    data = dfp.iloc[int(0),:]

    conf = GridConfig(30, None, None, n_doses)
    
    conf.primary_inoculum = dict(
        RR = data["RR"],
        RS = data["RS"],
        SR = data["SR"],
        SS = data["SS"],
        )

    fcide_pars = dict(
        omega_1 = data["omega_1"],
        omega_2 = data["omega_2"],
        theta_1 = data["theta_1"],
        theta_2 = data["theta_2"],
        delta_1 = data["delta_1"],
        delta_2 = data["delta_2"],
        )

    conf.load_saved = True

    conf.bs_sex_prop = bs
    conf.add_string()

    output = RunGrid(fcide_pars).run(conf)
    return output





def get_traces(output, bs, col):
    fRR = output.start_freqs_DA["RR"][-1,-1,:]
    fRS = output.start_freqs_DA["RS"][-1,-1,:]
    fSR = output.start_freqs_DA["SR"][-1,-1,:]
    fSS = output.start_freqs_DA["SS"][-1,-1,:]

    x = list(range(len(fRR)))

    traces = []

    op = 0.4
    # op = 1

    for ff, dash, opac, name in zip([fRR, fRS, fSR, fSS], 
                        ["solid", "dash", "dot", "solid"],
                        [1, op, op, op],
                        # ["RR", "RS", "SR", "SS"],
                        [f"<i>rr</i>, sexual proportion={bs}", f"<i>rs</i>, sexual proportion={bs}", f"<i>sr</i>, sexual proportion={bs}"],
                        ):



        ff = [None if ee==0 else log10(ee/(1-ee)) for ee in ff]
        trc = dict(x=x, y=ff,
                    line=dict(dash=dash, color=col),
                    opacity=opac, name=name,
                    legendgroup=name[:5],
                    # mode="markers+lines"
                    )
        traces.append(trc)
    
    return traces






# def plot_it(traces, index):
#     fig = go.Figure(data=traces, 
#             layout=dict(template="simple_white", 
#                         # title=f"sexual proportion={bs}"
#                         width = 800,
#                         height= 600,
#                         ))
    
#     fig.update_yaxes(range=[-20, 1], title="Logit frequency")
    
#     # fig.update_xaxes(range=[-0.5, 16.5], title="Time (years)")
#     fig.update_xaxes(range=[-0.5, 8.5], title="Time (years)")
#     # fig.update_xaxes(range=[-0.5, 13.5], title="Time (years)")

#     fig.update_layout(legend=dict(x=1, y=0, xanchor="right", yanchor="bottom"))
#     fig.show()

#     filename = f"./sr_hap/outputs/index={index}.png"
#     print(f"saving to {filename}")
#     fig.write_image(filename)



def get_haploid_outputs(n_its, n_doses, double_freq_factors, index, bss):
    outputs = []

    for bs in bss:

        dfp = get_sr_scan_params(n_its, double_freq_factors, index)

        output = check_run(n_doses, bs, dfp)

        outputs.append(output)
    return outputs




# if __name__=="__main__":

#     # n_its = 10
#     # n_doses = 6
#     # double_freq_factors = [1e-5, 1, 1e5]

#     outputs_l = get_haploid_outputs(n_its, n_doses, double_freq_factors, 6, [0,0.5,1])
#     outputs_m = get_haploid_outputs(n_its, n_doses, double_freq_factors, 10, [0,1])
#     outputs_h = get_haploid_outputs(n_its, n_doses, double_freq_factors, 24, [0,1])
    
    
    # traces = []
        
    # traces += get_traces(output, bs, col)
    
    # plot_it(traces, index)



    # cs = list(px.colors.sequential.Viridis)
    # cs = list(px.colors.sequential.Plotly3)

    # if len(bss)==2:
    #     cols = cs[:9:8]
    # elif len(bss)==3:
    #     cols = cs[:9:4]

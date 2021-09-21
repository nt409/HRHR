import plotly.colors as pltly_clrs

E_cols = pltly_clrs.n_colors("rgb(0,0,255)", "rgb(255,255,255)", 5, colortype="rgb")[:4]
        
I_cols = pltly_clrs.n_colors("rgb(255,0,0)", "rgb(255,255,255)", 5, colortype="rgb")[:4]


ATTRS_DICT = {
    'S': dict(color='limegreen', dash='solid', name='Susceptible', legendgrp="S"),
    'R': dict(color='rgb(100,100,100)', dash='solid', name='Removed', legendgrp="S"),
    
    'ERR': dict(color=E_cols[0], dash='dot', name='Exposed (rr)', legendgrp="E"),
    'ERS': dict(color=E_cols[1], dash='dash', name='Exposed (rs)', legendgrp="E"),
    'ESR': dict(color=E_cols[2], dash='dashdot', name='Exposed (sr)', legendgrp="E"),
    'ESS': dict(color=E_cols[3], dash='solid', name='Exposed (ss)', legendgrp="E"),

    'IRR': dict(color=I_cols[0], dash='dot', name='Infectious (rr)', legendgrp="I"),
    'IRS': dict(color=I_cols[1], dash='dash', name='Infectious (rs)', legendgrp="I"),
    'ISR': dict(color=I_cols[2], dash='dashdot', name='Infectious (sr)', legendgrp="I"),
    'ISS': dict(color=I_cols[3], dash='solid', name='Infectious (ss)', legendgrp="I"),

    'fung_1': dict(color='turquoise', dash='solid', name='Fungicide <i>A</i>', legendgrp="F"),
    'fung_2': dict(color='magenta', dash='dot', name='Fungicide <i>B</i>', legendgrp="F"),
}


LABEL_COLOR = "rgb(110,110,110)"

NULL_HEATMAP_COLOUR = "rgb(100, 100, 100)"

STRAIN_ATTRS = dict(SS=dict(color="rgb(34,140,34)", dash="dot", longname='Double sensitive', abbrv="SS"),
                RS=dict(color="rgb(20,20,200)", dash="dash", longname='Single resistant (RS)', abbrv="RS"),
                SR=dict(color="rgb(200,20,20)", dash="dashdot", longname='Single resistant (SR)', abbrv="SR"),
                RR=dict(color="rgb(50,50,50)", dash="solid", longname='Double resistant', abbrv="RR"),
                )

TITLE_MAP = dict(
    LTY = "Lifetime yield",
    TY = "Total yield",
    FY = "Effective life",
    econ = "Economic life",
    )


# default is half page
PLOT_WIDTH = 600
PLOT_HEIGHT = 400

FULL_PAGE_WIDTH = 800
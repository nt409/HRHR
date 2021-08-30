import plotly.colors as pltly_clrs

E_cols = pltly_clrs.n_colors("rgb(0,0,255)", "rgb(255,255,255)", 5, colortype="rgb")[:4]
        
I_cols = pltly_clrs.n_colors("rgb(255,0,0)", "rgb(255,255,255)", 5, colortype="rgb")[:4]


ATTRS_DICT = {
    'S': dict(color='limegreen', dash='solid', name='Susceptible'),
    'R': dict(color='rgb(100,100,100)', dash='solid', name='Removed'),

    'ERR': dict(color=E_cols[0], dash='dot', name='E (rr)'),
    'ERS': dict(color=E_cols[1], dash='dash', name='E (rs)'),
    'ESR': dict(color=E_cols[2], dash='dashdot', name='E (sr)'),
    'ESS': dict(color=E_cols[3], dash='solid', name='E (ss)'),

    'IRR': dict(color=I_cols[0], dash='dot', name='I (rr)'),
    'IRS': dict(color=I_cols[1], dash='dash', name='I (rs)'),
    'ISR': dict(color=I_cols[2], dash='dashdot', name='I (sr)'),
    'ISS': dict(color=I_cols[3], dash='solid', name='I (ss)'),

    'fung_1': dict(color='turquoise', dash='solid', name='Fungicide A'),
    'fung_2': dict(color='magenta', dash='dot', name='Fungicide B'),
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
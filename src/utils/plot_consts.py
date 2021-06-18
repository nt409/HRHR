from .params import PARAMS

ATTRS_DICT = {
    str(PARAMS.S_ind): dict(name='Susceptible', colour='green', dash="solid"),
    str(PARAMS.ER_ind): dict(name='Latent (RR)', colour='rgb(0,0,150)', dash="dot"),
    str(PARAMS.ERS_ind): dict(name='Latent (RS)', colour='rgb(0,0,180)', dash="dash"),
    str(PARAMS.ESR_ind): dict(name='Latent (SR)', colour='rgb(0,0,210)', dash="dashdot"),
    str(PARAMS.ES_ind): dict(name='Latent (SS)', colour='rgb(0,0,240)', dash="solid"),
    str(PARAMS.IR_ind): dict(name='Infected (RR)', colour='rgb(150,0,0)', dash="dot"),
    str(PARAMS.IRS_ind): dict(name='Infected (RS)', colour='rgb(180,0,0)', dash="dash"),
    str(PARAMS.ISR_ind): dict(name='Infected (SR)', colour='rgb(240,0,0)', dash="dashdot"),
    str(PARAMS.IS_ind): dict(name='Infected (SS)', colour='rgb(240,0,0)', dash="solid"),
    str(PARAMS.R_ind): dict(name='Removed', colour='black', dash="solid"),
    str(PARAMS.PR_ind): dict(name='Primary inoc. (RR)', colour='rgb(150,0,150)', dash="dot"),
    str(PARAMS.PRS_ind): dict(name='Primary inoc. (RS)', colour='rgb(180,0,180)', dash="dash"),
    str(PARAMS.PSR_ind): dict(name='Primary inoc. (SR)', colour='rgb(210,0,210)', dash="dashdot"),
    str(PARAMS.PS_ind): dict(name='Primary inoc. (SS)', colour='rgb(240,0,240)', dash="solid"),
    str(PARAMS.Fung1_ind): dict(name='Fung. 1 conc.', colour='rgb(120,240,0)', dash="solid"),
    str(PARAMS.Fung2_ind): dict(name='Fung. 2 conc.', colour='rgb(240,120,0)', dash="solid"),
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
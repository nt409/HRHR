from math import floor
import plotly.graph_objects as go
import plotly.express as px
import numpy as np

from .consts import LIGHT_GREY_TEXT, NULL_HEATMAP_COLOUR, PLOT_WIDTH, PLOT_HEIGHT, LABEL_COLOR


def standard_layout(legend_on, width=PLOT_WIDTH, height=PLOT_HEIGHT):
    return go.Layout(
        font=dict(size=16),
        template="simple_white",
        width=width,
        height=height,
        showlegend=legend_on,
        xaxis=dict(showgrid=False),
    )


def grey_colorscale(z):

    out = []
    pltly_clr_scale = list(px.colors.sequential.Inferno)

    pal = [NULL_HEATMAP_COLOUR, NULL_HEATMAP_COLOUR] + pltly_clr_scale
    vals = [0, 1/np.amax(z)] + list(np.linspace(1/np.amax(z), 1, len(pal)-2))

    for val, col in zip(vals, pal):
        out.append([val, col])

    return out


def grey_colorscale_discrete(z):
    out = []

    n_cols = int(np.amax(z)-1)

    all_cols = get_N_colors(n_cols)

    # last color doesn't matter
    cols = [NULL_HEATMAP_COLOUR, NULL_HEATMAP_COLOUR,
            NULL_HEATMAP_COLOUR] + all_cols + ["white"]

    vals = [0, 0.5/np.amax(z)] + list(np.linspace(0.5 /
                                                  np.amax(z), 1 - 0.5/np.amax(z), n_cols)) + [1]

    opt_col = "limegreen"

    cols[-2:] = [opt_col]*2

    for ii in range(len(vals)):
        out.append([vals[ii], cols[ii]])
        out.append([vals[ii], cols[ii+1]])

    return out


def divergent_color_scale(z, x):

    out = []

    # clrs = list(px.colors.sequential.Inferno)

    low_col = "indigo"
    hi_col = "seagreen"

    cols = [low_col, "white", hi_col]

    midpoint = (x - np.nanmin(z))/(np.nanmax(z) - np.nanmin(z))

    vals = [0, midpoint, 1]

    for ii in range(len(vals)):
        out.append([vals[ii], cols[ii]])
        # out.append([vals[ii], cols[ii+1]])

    return out


def grey_colorscale_discrete_N(max_N):
    out = []

    n_cols = int(max_N-1)

    all_cols = get_N_colors(n_cols)

    # last color doesn't matter
    cols = [NULL_HEATMAP_COLOUR, NULL_HEATMAP_COLOUR,
            NULL_HEATMAP_COLOUR] + all_cols + ["white"]

    vals = [0, 0.5/max_N] + \
        list(np.linspace(0.5/max_N, 1 - 0.5/max_N, n_cols)) + [1]

    # cols = [NULL_HEATMAP_COLOUR, NULL_HEATMAP_COLOUR] + pltly_clr_scale + ["red"]*5
    # len_vals = int(max_N) - 1
    # vals = [0, 1/max_N] + list(np.linspace(1/max_N, 1, n_cols))

    for ii in range(len(vals)):
        out.append([vals[ii], cols[ii]])
        out.append([vals[ii], cols[ii+1]])

    return out


def get_N_colors(n_cols):

    clrs = list(px.colors.sequential.Inferno)

    spacing = np.linspace(0, 9, n_cols)

    intermed = [
        px.colors.find_intermediate_color(
            px.colors.hex_to_rgb(clrs[floor(xx)]),
            px.colors.hex_to_rgb(clrs[floor(xx)+1]),
            xx - floor(xx)) for xx in spacing[:-1]]

    all_cols = [px.colors.label_rgb(cc) for cc in intermed] + [clrs[-1]]

    return all_cols


def my_colorbar(title):
    return dict(
        title=title,
        titleside='right',
    )


def my_colorbar_subplot(title, x, y, length):
    return dict(
        title=title,
        x=x,
        y=y,
        len=length,
        thickness=25,
        # titleside = 'right',
    )


def invisible_colorbar(x, y):
    """
    hacky way to remove second colorbar - set x position so not visible
    """
    return dict(
        x=x,
        y=y,
        len=0.1,
        thickness=1,
        tickfont=dict(size=1,
                      color="rgba(0,0,0,0)"
                      )
    )


def get_text_annotation(x, y, text, xanchor=None, yanchor=None):

    xanchor = "center" if xanchor is None else xanchor
    yanchor = "top" if yanchor is None else yanchor

    return dict(
        x=x,
        y=y,
        text=text,

        showarrow=False,

        xref='paper',
        yref='paper',

        xanchor=xanchor,
        yanchor=yanchor,

        font=dict(
            size=14,
            color=LABEL_COLOR,
        ),
    )


def get_big_text_annotation(x, y, text, xanchor=None):
    if xanchor is None:
        xanchor = "center"

    return dict(
        x=x,
        y=y,
        text=text,

        showarrow=False,

        xref='paper',
        yref='paper',

        xanchor=xanchor,
        yanchor="top",

        font=dict(
            size=30,
            color=LIGHT_GREY_TEXT,
        ),
    )


def dose_space_annotation(text, x, y, ax, ay, subplot, color):
    return dict(text=text,
                x=x,
                y=y,
                ax=ax,
                ay=ay,
                xref=f"x{int(subplot)}",
                yref=f"y{int(subplot)}",
                # showarrow=False,
                # arrowwidth=4,
                arrowsize=1,
                arrowhead=1,
                arrowcolor=color,
                font=dict(size=16, color=color)
                )


def get_shape_annotation(x, y, fontsize=None, xanchor=None):

    xanchor = "center" if xanchor is None else xanchor
    fontsize = 50 if fontsize is None else fontsize

    return dict(
        x=x,
        y=y,
        text=" ",

        showarrow=False,

        bgcolor="white",

        opacity=1,

        xref='paper',
        yref='paper',

        xanchor=xanchor,
        yanchor="top",

        font=dict(
            size=fontsize,
        ),
    )


def get_arrow_annotation(x, y, dx, dy):
    return dict(
        x=x,
        y=y,
        # text=text,

        showarrow=True,
        arrowcolor=LABEL_COLOR,
        arrowsize=2,
        arrowwidth=1,
        arrowhead=2,

        ax=dx,
        ay=dy,

        xref='paper',
        yref='paper',

        xanchor="center",
        yanchor="top",

        font=dict(
            size=14,
            color=LABEL_COLOR,
        ),
    )


# End of utility functions

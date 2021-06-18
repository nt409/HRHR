import plotly.graph_objects as go
import plotly.express as px
import numpy as np

from .plot_consts import NULL_HEATMAP_COLOUR, PLOT_WIDTH, PLOT_HEIGHT, LABEL_COLOR


def standard_layout(legend_on, width=PLOT_WIDTH, height=PLOT_HEIGHT):
    return go.Layout(
            font = dict(size=16),
            template="plotly_white",
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


def my_colorbar(title):
    return dict(
        title = title,
        titleside = 'right',
        )


def invisible_colorbar(x):
    """
    hacky way to remove second colorbar - set x position so not visible
    """
    return dict(x=x, len=0.1, 
            tickfont=dict(size=1,
                color="rgba(0,0,0,0)"
                ))


def get_text_annotation(x, y, text):
    return dict(
            x=x,
            y=y,
            text=text,
            
            showarrow=False,
                
            xref='paper',
            yref='paper',

            xanchor="center",
            yanchor="top",

            font=dict(
                    size=14,
                    color=LABEL_COLOR,
                ),
        )


def get_big_text_annotation(x, y, text):
    return dict(
            x=x,
            y=y,
            text=text,
            
            showarrow=False,
                
            xref='paper',
            yref='paper',

            xanchor="center",
            yanchor="top",

            font=dict(
                    size=30,
                    color="rgb(150,150,150)",
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

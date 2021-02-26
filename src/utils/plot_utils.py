import plotly.graph_objects as go
import numpy as np
from .plot_consts import NULL_HEATMAP_COLOUR


def standard_layout(legend_on):
    return go.Layout(
            font = dict(size=22),
            template="plotly_white",
            width=1400,
            height=700,
            showlegend=legend_on,
            xaxis=dict(showgrid=False),
            )

def grey_colorscale(z):
    return [
        [0, NULL_HEATMAP_COLOUR],
        [1/np.amax(z), NULL_HEATMAP_COLOUR],
        [1/np.amax(z), "rgb(0, 0, 100)"],
        [1, "rgb(255, 255, 0)"],
    ]

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

# End of utility functions

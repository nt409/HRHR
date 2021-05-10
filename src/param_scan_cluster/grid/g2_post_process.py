"""
Combine the dataframes.
"""
from ..config import config_grid
from ..functions import combine_PS_grid_outputs



def main(config):
    combine_PS_grid_outputs(config)


if __name__=="__main__":
    main(config_grid)
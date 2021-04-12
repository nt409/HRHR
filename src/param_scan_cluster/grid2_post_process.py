"""
Combine the dataframes.
"""
from .config import config
from .functions import combine_PS_grid_outputs



def main():
    combine_PS_grid_outputs(config)


if __name__=="__main__":
    main()
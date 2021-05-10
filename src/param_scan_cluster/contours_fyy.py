import matplotlib.pyplot as plt
import numpy as np
from math import log, exp

from runHRHR.config_classes import GridConfig
from utils.functions import RunGrid
from .functions import get_contours




def get_contour_levels(z, num=5):
    """
    takes an array z and returns contour levels from 
    95 to maximum value, logarithmically spaced
    """
    max_val = z[-1,-1]

    exp_levels = np.linspace(exp(95), exp(max_val), num=num, endpoint=False)

    levels = [log(ii) for ii in exp_levels]

    return levels


def plot_contours(contours):
    _, ax = plt.subplots(1)

    for data in contours:
        print(data)
        print("\n")
        ax.plot(data['x'], data['y'])

    plt.show()



def main():
    Conf = GridConfig(1, 0.001, 0.001, 11)

    Conf.load_saved = True

    output = RunGrid().grid_of_tactics(Conf)

    fyy_array = output['yield_array'][:,:,0]

    levels = get_contour_levels(fyy_array)

    out = get_contours(fyy_array, levels)

    print(out)




if __name__=="__main__":
    main()


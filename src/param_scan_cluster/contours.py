import matplotlib.pyplot as plt
import numpy as np
from math import log, exp

from runHRHR.config_classes import GridConfig
from utils.functions import RunGrid



def get_contours(z, levels=None):
    """
    Takes an array z and returns a list of dictionaries containing:
    - x values
    - y values
    - the contour level
    """
    
    x, y = np.mgrid[0:1:z.shape[0]*1j, 0:1:z.shape[1]*1j]
    
    cs = plt.contour(x, y, z, levels=levels)

    output = []

    for ind in range(len(cs.allsegs)):
        cont = cs.allsegs[ind][0]
        
        vals = dict(x=cont[:,0], y=cont[:,1])
            
        int_ = 1
        
        output.append(
            dict(
                x = [vals['x'][0]] + list(vals['x'][1:-2:int_]) + [vals['x'][-1]],
                y = [vals['y'][0]] + list(vals['y'][1:-2:int_]) + [vals['y'][-1]],
                level = levels[ind]
                )
            )
    return output


def get_contour_levels(z, num=5):
    """
    takes an array z and returns contour levels from 
    95 to maximum value, logarithmically spaced
    """
        
    max_val = z[-1,-1]

    exp_levels = np.linspace(exp(95), exp(max_val), num=num, endpoint=False)

    levels = [log(ii) for ii in exp_levels]

    return levels


Conf = GridConfig(1, 0.001, 0.001, 11)

Conf.load_saved = True

output = RunGrid().grid_of_tactics(Conf)

fyy_array = output['yield_array'][:,:,0]

levels = get_contour_levels(fyy_array)

out = get_contours(fyy_array, levels)

print(out)

exit()

fig, ax = plt.subplots(1)

for data in out:
    print(data)
    print("\n")
    ax.plot(data['x'], data['y'])

plt.show()
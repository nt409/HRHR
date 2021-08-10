import pandas as pd
import numpy as np
import itertools

from utils.functions import EqualResFreqBreakdownArray, \
     EqualSelectionArray

# TOC
# CheckStrategyUsingIVT_DF
# ContourPassesThroughChecker



class CheckStrategyUsingIVT_DF:
    """
    Use "intermediate value theorem" method:
    
    - find optimal region and see if contains points above and below contour
    - if so, then contour passes through (provided connected)

    """
    def __init__(self, grid_output, strat_name) -> None:
        
        print(f"Running IVT method: {strat_name}")

        self.FYs = grid_output['FY']
        
        self.name = strat_name

        if strat_name=="EqSel":
            self.level = 0.5
            strategy_class = EqualSelectionArray

        elif strat_name=="RFB":
            self.level = 0
            strategy_class = EqualResFreqBreakdownArray

        else:
            raise Exception(f"invalid strat_name: {strat_name}")

        self.strategy_obj = strategy_class(grid_output)
        
        self.check_if_gives_optimum()

        self.find_best_value_this_strat()

        self.df = self.get_df(strat_name)




    def check_if_gives_optimum(self):
        FYs = self.FYs
        
        opt_region = FYs!=np.amax(FYs)
        
        if not self.strategy_obj.is_valid:
            self.strat_works = f"Strategy {self.name} is invalid"
            return None
        
        strat_array = self.strategy_obj.array

        opt_strat = np.ma.masked_array(strat_array, mask=opt_region)
        
        self.strat_works = ContourPassesThroughChecker(opt_strat, self.level).passes_through



    def find_best_value_this_strat(self):
        FYs = self.FYs

        if not self.strategy_obj.is_valid:
            self.best_value = f"Strategy {self.name} is invalid"
            return None

        strat_array = self.strategy_obj.array

        for k in range(int(np.amax(FYs))):
            EL = np.amax(FYs)-k
            
            opt_region = FYs<EL

            opt_strat = np.ma.masked_array(strat_array, mask=opt_region)

            worked = ContourPassesThroughChecker(opt_strat, self.level).passes_through
        
            if worked:
                self.best_value = EL
                return None

        self.best_value = EL


    def get_df(self, strat_name):
        data = {f"I_{strat_name[0]}_best_region": self.strat_works,
                f"I_{strat_name[0]}_best_value": self.best_value}

        return pd.DataFrame([data])












class ContourPassesThroughChecker:
    def __init__(self, array, level) -> None:
        
        self.passes_through = self.check_if_passes_through_region(array, level)




    def check_if_passes_through_region(self, array, level) -> bool:
        """
        Check if contour passes through optimal region.
        
        If False, might be because:
        - optimal region is not connected, or
        - all contour values are positive/negative in the region
        - grid sufficiently dense

        If True, contour passes through optimal region :)

        """
        
        includes_level = (np.amin(array)<level and np.amax(array)>level)

        if not includes_level:
            return False
        
        return self.check_is_connected(array, level)



    
    def check_is_connected(self, array, level) -> bool:

        for i, j in itertools.product(*(range(array.shape[ii]) for ii in [0,1])):

            if not array[i,j]:
                continue

            elif array[i,j]==level:
                return True

            else:
                valid = self.check_neighbours(array, i, j, level)
            
            if valid:
                return True

        return False
    


    def check_neighbours(self, array, i, j, level) -> bool:
        """
        For array only including optimal region:

        Check if a cell has any neighbours which have a differing sign for contour.

        If True, there is some intermediate value for doses between the two cells
        which is optimal and along the contour.

        If False, t
        """
        
        self_is_positive = array[i,j]>level

        cell_iterator = self.neighbours((i,j), array.shape[0])

        for ii, jj in cell_iterator:
                
            if not array[ii,jj]:
                continue

            neighbour_is_positive = array[ii,jj]>level

            if neighbour_is_positive!=self_is_positive:
                return True

        return False

    
    @staticmethod
    def neighbours(cell, shape):
        """
        Find neighbouring cells, excluding self and excluding any cells with:
        -: index < 0, or
        -: index > shape.

        Returns a generator

        """
        for new_cell in itertools.product(*(range(n-1, n+2) for n in cell)):
            if new_cell != cell and all(0 <= n < shape for n in new_cell):
                yield new_cell



import pandas as pd
import numpy as np


class OtherStratsDF:
    def __init__(self, grid_output) -> None:
        FYs = grid_output.FY

        self.df = self._get_strategy_outcomes(FYs)

    def _get_strategy_outcomes(self, FYs):

        minEqDoseELVec = [float(FYs[i, i]) for i in range(FYs.shape[0])]

        data = dict(
            max_grid_EL=np.amax(FYs),
            O_minEqDoseEL=self._get_first_geq_2_element(minEqDoseELVec),
            O_corner_00=FYs[0, 0],
            O_corner_10=FYs[-1, 0],
            O_corner_01=FYs[0, -1],
            O_fullDoseEL=FYs[-1, -1])

        out = pd.DataFrame([data])

        return out

    @staticmethod
    def _get_first_geq_2_element(vec):
        filtered = list(filter(lambda x: x >= 2, vec))

        if not filtered:
            return "NA"

        return filtered[0]

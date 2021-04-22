"""
Analyse the parameter scan output
"""
import pandas as pd

from .config import config_rand
from .functions import get_PS_rand_str
from utils.functions import non_decreasing, non_increasing


# def main(config):
    
par_str = get_PS_rand_str(config_rand)

df = pd.read_csv(f"param_scan_cluster/outputs/PS_combined_rand_output_{par_str}.csv")

# print(df)
grouped_df = df.groupby(["run", "MS"])




def monotone_RFB(data):
    sorted_data = data.sort_values(["delta_RFB"])

    dfPos = sorted_data[sorted_data["delta_RFB"]>=0]
    dfNeg = sorted_data[sorted_data["delta_RFB"]<0]

    ELPos = dfPos["EL"]
    ELNeg = dfNeg["EL"]
    
    out = non_increasing(ELPos) and non_decreasing(ELNeg)

    return out



gdf = grouped_df.apply(monotone_RFB)

print(gdf)
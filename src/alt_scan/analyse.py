import itertools
import pandas as pd

n_its = 5
n_doses = 11




df = pd.DataFrame()

for alt_strat, within_season in itertools.product(["alt_12", "alt_21"], [True, False]):

    filename = f"./alt_scan/output_{n_its}_{n_doses}_{within_season}_{alt_strat}.csv"
    
    new_df = pd.read_csv(filename)
    new_df["within_season"] = within_season
    new_df["alt_strat"] = alt_strat

    df = pd.concat([df, new_df])



# print(df.head(-10))
# print(df.describe())

print(df[df["asex_alt"]>df["asex_mix"]])
print(df[df["sex_alt"]>df["sex_mix"]])
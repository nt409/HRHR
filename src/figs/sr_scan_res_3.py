import pandas as pd

from sr_scan.configs import config_res
from plotting.paper_figs import SREffectResults3Panel


if __name__ == "__main__":

    n_its = config_res["n_its"]
    n_sex_props = config_res["n_sex_props"]
    n_doses = config_res["n_doses"]
    double_freq_factors = config_res["double_freq_factors"]

    dff_str = ",".join([str(ee) for ee in double_freq_factors])
    filename = f"./sr_scan/outputs/combined/df_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.csv"
    df = pd.read_csv(filename)

    df = df.loc[df["run"] < 5]

    # change order via middle step
    df.run.replace(0, 50, inplace=True)
    df.run.replace(1, 51, inplace=True)
    df.run.replace(2, 52, inplace=True)
    df.run.replace(3, 53, inplace=True)

    df.run.replace(52, 0, inplace=True)
    df.run.replace(53, 1, inplace=True)
    df.run.replace(50, 2, inplace=True)
    # df.run.replace(51, 3, inplace=True)

    df_use = (
        df
        .loc[df["run"] <= 3]
        .sort_values(["run", "bs_sex_prop"])
    )

    df.sort_values(by=["run", "bs_sex_prop"], inplace=True)

    filename = f"../outputs/figures/paper_figs/sr_effect_res_{n_its}_{n_sex_props}_{n_doses}_{dff_str}.png"
    p = SREffectResults3Panel(df_use, double_freq_factors, filename)
    p.fig.write_image('../outputs/figures/paper_figs/hires/fig5.jpg', scale=4)

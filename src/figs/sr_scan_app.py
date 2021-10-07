from sr_hap.utils import get_haploid_outputs_app
from sr_hap.configs import config_app
from plotting.paper_figs import SREffectAppendix



if __name__=="__main__":

    n_variants = config_app["n_variants"]
    n_sex_props = config_app["n_sex_props"]
    n_doses = config_app["n_doses"]
    double_freq_factor_lowest = config_app["double_freq_factor_lowest"]


    indices = [4, 13, 22]

    outputs_l = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[0], [0,0.2,1])
    outputs_m = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[1], [0,0.2,1])
    outputs_h = get_haploid_outputs_app(n_variants, n_doses, double_freq_factor_lowest, indices[2], [0,0.2,1])

    filename = f"../outputs/figures/paper_figs/sr_effect_app_{n_variants}_{n_sex_props}_{n_doses}_{double_freq_factor_lowest}.png"
    SREffectAppendix(outputs_l, outputs_m, outputs_h, filename)
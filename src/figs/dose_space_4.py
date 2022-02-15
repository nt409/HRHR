"""Now created by notebooks/007_dose_space_4.ipynb"""


from model.params import PARAMS
from model.simulator import RunGrid
from model.config_classes import GridConfig
from model.utils import primary_inoc_from_sing_res_strains


def get_data(load_saved, n_doses):

    primary_inoc_same = primary_inoc_from_sing_res_strains(1e-7, 1e-7)

    primary_inoc_diff = primary_inoc_from_sing_res_strains(1e-7, 1e-3)

    fcide_p_B = dict(
        omega_1=1,
        omega_2=0.95,
        theta_1=PARAMS.theta_1,
        theta_2=7.5,
        delta_1=PARAMS.delta_1,
        delta_2=PARAMS.delta_2
    )

    fcide_p_D = dict(
        omega_1=1,
        omega_2=0.8,
        theta_1=PARAMS.theta_2,
        theta_2=7,
        delta_1=PARAMS.delta_1,
        delta_2=PARAMS.delta_2
    )

    # same RFs, default fung pars
    output_SS, _ = _data_one_panel(
        primary_inoc_same,
        None,
        n_doses,
        load_saved,
    )

    # diff RFs, default fung pars
    output_DS, _ = _data_one_panel(
        primary_inoc_diff,
        None,
        n_doses,
        load_saved,
    )

    # same RFs, diff fung pars
    output_SD, _ = _data_one_panel(
        primary_inoc_same,
        fcide_p_B,
        n_doses,
        load_saved,
    )

    # diff RFs, diff fung pars
    output_DD, conf_grid = _data_one_panel(
        primary_inoc_diff,
        fcide_p_D,
        n_doses,
        load_saved,
    )

    print(f"{primary_inoc_same=}")
    print(f"{primary_inoc_diff=}")
    print(f"{fcide_p_B=}")
    print(f"{fcide_p_D=}")

    return (
        output_SS,
        output_SD,
        output_DS,
        output_DD,
        conf_grid.config_string_img,
    )


def _data_one_panel(primary_inoculum, fcide_pars, n_doses, load_saved):
    conf_grid = GridConfig(30, None, None, n_doses)
    conf_grid.load_saved = load_saved
    conf_grid.primary_inoculum = primary_inoculum
    conf_grid.add_string()

    if fcide_pars is not None:
        fung_par_str = (
            f"_fung_pars="
            f"{fcide_pars['omega_1']},{fcide_pars['omega_2']},"
            f"{fcide_pars['theta_1']},{fcide_pars['theta_2']}"
        )

        conf_grid.config_string = (
            conf_grid.config_string.split(".pickle")[0] +
            f"{fung_par_str.replace('.', ',')}.pickle"
        )

    grid_out = RunGrid(fcide_pars).run(conf_grid)
    return grid_out, conf_grid

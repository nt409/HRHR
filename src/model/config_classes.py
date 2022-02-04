
class BaselineConfig:
    def __init__(self,
                 n_years,
                 rp1,
                 rp2,
                 primary_inoculum,
                 bs_sex_prop=None,
                 ws_sex_prop=None,
                 ):

        self.load_saved = True
        self.save = True

        self.folder_runs = '../outputs/saved_runs/'
        self.folder_figs = '../outputs/figures/'

        self.bs_sex_prop = 0 if bs_sex_prop is None else bs_sex_prop
        self.ws_sex_prop = 0 if ws_sex_prop is None else ws_sex_prop

        self.n_years = n_years

        self.res_props = dict(f1=rp1, f2=rp2)

        self.primary_inoculum = primary_inoculum

    def add_baseline_str(self):

        inoc_str = "N" if self.primary_inoculum is None else (
            f"{round(self.primary_inoculum['RS'],10)},_"
            f"{round(self.primary_inoculum['SR'],10)},_"
            f"{round(self.primary_inoculum['RR'],15)}")

        save_str = (f"Ny={self.n_years}_"
                    f"RPs={self.res_props['f1']},_{self.res_props['f2']}_"
                    f"PI={inoc_str}_"
                    f"BSS={round(self.bs_sex_prop,2)}_"
                    f"WSS={round(self.ws_sex_prop,2)}"
                    )

        save_str = save_str.replace("0.", "")
        save_str = save_str.replace(".", ",")
        save_str = save_str.replace("__", "_")
        save_str = save_str.replace("None", "NA")

        self.save_string = save_str

    def get_conf_strings(self, filename):
        comma_free_name = filename.replace('.', ',')
        self.config_string = f"{self.folder_runs}{comma_free_name}.pickle"
        self.config_string_img = f"{self.folder_figs}{comma_free_name}.png"


class SingleConfig(BaselineConfig):
    def __init__(self,
                 n_years,
                 total_res_f1,
                 total_res_f2,
                 dose_11,
                 dose_12,
                 dose_21,
                 dose_22,
                 primary_inoculum=None,
                 bs_sex_prop=None,
                 ws_sex_prop=None,
                 ):
        """Config for a single model run

        Supply either primary inoclum or total res.

        Parameters
        ----------
        n_years : int
            n years to run model for
        total_res_f1 : float
            Total amount of resistance to f1
        total_res_f2 : float
            Total amount of resistance to f2
        dose_11 : list, np.array
            First dose f1 in year 1, 2, ... len(dose_11)
        dose_12 : list, np.array
            Second dose f1 in year 1, 2, ... len(dose_12)
        dose_21 : list, np.array
            First dose f2 in year 1, 2, ... len(dose_21)
        dose_22 : list, np.array
            Second dose f2 in year 1, 2, ... len(dose_22)
        primary_inoculum : dict, optional
            dict, keys RR, RS, SR, SS, by default None
        bs_sex_prop : float, optional
            between-season sex proportion
        ws_sex_prop : float, optional
            between-season sex proportion

        Example
        -------
        >>>primary_inoc = dict(
        ... RR=1e-5,
        ... RS=1e-3,
        ... SR=1e-5,
        ... SS=1-1e-5-1e-3-1e-5
        ... )
        >>>config = SingleConfig(
        ... 35,
        ... None, None,
        ... dose1, dose1, dose2, dose2,
        ... primary_inoculum=primary_inoc
        ... )
        """

        super().__init__(
            n_years,
            total_res_f1,
            total_res_f2,
            primary_inoculum,
            bs_sex_prop,
            ws_sex_prop,
        )

        self.fung1_doses = dict(spray_1=[dose_11]*n_years,
                                spray_2=[dose_12]*n_years)

        self.fung2_doses = dict(spray_1=[dose_21]*n_years,
                                spray_2=[dose_22]*n_years)

        self.add_string()

    def add_string(self):
        d11 = round(self.fung1_doses['spray_1'][0], 6)
        d12 = round(self.fung1_doses['spray_2'][0], 2)
        d21 = round(self.fung2_doses['spray_1'][0], 6)
        d22 = round(self.fung2_doses['spray_2'][0], 2)

        self.add_baseline_str()

        filename = f"single/{self.save_string}_doses={d11},{d12},{d21},{d22}"

        self.get_conf_strings(filename)


class GridConfig(BaselineConfig):
    def __init__(self,
                 n_years,
                 rp1,
                 rp2,
                 n_doses,
                 primary_inoculum=None
                 ):

        super().__init__(n_years, rp1, rp2, primary_inoculum)

        self.strategy = 'mix'

        self.n_doses = n_doses

        self.add_string()

    def add_string(self):
        self.add_baseline_str()

        filename = f"grid/{self.save_string}_Nd={self.n_doses}_S={self.strategy}"

        self.get_conf_strings(filename)

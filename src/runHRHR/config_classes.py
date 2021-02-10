def get_conf_string(folder, filename):
    conf_str = folder + filename.replace(".", ",") + ".pickle"
    config_string = conf_str

    conf_str2 = conf_str.replace("saved_runs", "figures")
    conf_str2 = conf_str2.replace("pickle", "png")
    config_string_img = conf_str2

    return config_string, config_string_img


class BaselineConfig:
    def __init__(self, n_years, rp1, rp2, primary_inoculum, zeroth_season_reproduction):

        self.load_saved = False

        self.folder_save_run = '../outputs/saved_runs/'

        self.sex_prop = 0
        self.is_mixed_sex = True
        
        self.n_years = n_years

        self.res_props = dict(
            f1 = rp1,
            f2 = rp2
            )

        self.primary_inoculum = primary_inoculum

        self.zeroth_season_reproduction = zeroth_season_reproduction

        self.save_string = f"Ny={n_years}" + \
            f"Rps={rp1},_{rp2}_" + \
            f"PI={str(primary_inoculum)[0]}_" + \
            f"0th={str(zeroth_season_reproduction)[0]}_" + \
            f"Sex={self.sex_prop}"



class SingleConfig(BaselineConfig):
    def __init__(self,
                n_years,
                rp1,
                rp2,
                d11,
                d12,
                d21,
                d22,
                primary_inoculum=None,
                zeroth_season_reproduction=True
                ):
        
        super().__init__(n_years, rp1, rp2, primary_inoculum, zeroth_season_reproduction)
        
        self.fung1_doses = dict(
            spray_1 = [d11]*n_years,
            spray_2 = [d12]*n_years
            )

        self.fung2_doses = dict(
            spray_1 = [d21]*n_years,
            spray_2 = [d22]*n_years
            )
        
        self.add_string()

    def add_string(self):
        d11 = self.fung1_doses['spray_1'][0]
        d12 = self.fung1_doses['spray_2'][0]
        d21 = self.fung2_doses['spray_1'][0]
        d22 = self.fung2_doses['spray_2'][0]

        filename = f"single/{self.save_string}_doses={d11},{d12},{d21},{d22}"
        self.config_string, self.config_string_img = get_conf_string(self.folder_save_run, filename)




class GridConfig(BaselineConfig):
    def __init__(self,
            n_years,
            rp1,
            rp2,
            n_doses,
            primary_inoculum=None,
            zeroth_season_reproduction=True
            ):

        super().__init__(n_years, rp1, rp2, primary_inoculum, zeroth_season_reproduction)
        
        self.strategy = 'mix'

        self.n_doses = n_doses

        self.add_string()

    def add_string(self):
        filename = f"grid/{self.save_string}_Nd={self.n_doses}_S={self.strategy}"
        self.config_string, self.config_string_img = get_conf_string(self.folder_save_run, filename)
        

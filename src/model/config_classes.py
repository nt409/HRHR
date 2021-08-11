def get_conf_string(folder, filename):
    conf_str = folder + filename.replace(".", ",") + ".pickle"
    config_string = conf_str

    conf_str2 = conf_str.replace("saved_runs", "figures")
    config_string_img = conf_str2.replace("pickle", "png")

    return config_string, config_string_img


class BaselineConfig:
    def __init__(self, n_years, rp1, rp2, primary_inoculum):

        self.load_saved = True

        self.folder_save_run = '../outputs/saved_runs/'

        self.sex_prop = 0
        
        self.n_years = n_years

        self.res_props = dict(
            f1 = rp1,
            f2 = rp2
            )

        self.primary_inoculum = primary_inoculum



    def add_baseline_str(self):

        rp1 = self.res_props['f1']
        rp2 = self.res_props['f2']

        if self.primary_inoculum is None:
            inoc_str = "N"
        else:
            inoc_str = "Y"

        self.save_string = f"Ny={self.n_years}_" + \
            f"Rps={rp1},_{rp2}_" + \
            f"PI={inoc_str}_" + \
            f"Sex={self.sex_prop}"
            # f"0th={str(self.zeroth_season_reproduction)[0]}_" + \



class SingleConfig(BaselineConfig):
    def __init__(self,
                n_years,
                rp1,
                rp2,
                d11,
                d12,
                d21,
                d22,
                primary_inoculum=None
                # zeroth_season_reproduction=True
                ):
        
        super().__init__(n_years, rp1, rp2, primary_inoculum)
        
        self.fung1_doses = dict(
            spray_1 = [d11]*n_years,
            spray_2 = [d12]*n_years
            )

        self.fung2_doses = dict(
            spray_1 = [d21]*n_years,
            spray_2 = [d22]*n_years
            )
        
        self.add_string()



    def add_string(self, extra_detail=None):
        d11 = round(self.fung1_doses['spray_1'][0],6)
        d12 = round(self.fung1_doses['spray_2'][0],2)
        d21 = round(self.fung2_doses['spray_1'][0],6)
        d22 = round(self.fung2_doses['spray_2'][0],2)

        self.add_baseline_str()
        filename = f"single/{self.save_string}_doses={d11},{d12},{d21},{d22}"
        if extra_detail is not None:
            filename = filename + extra_detail

        self.config_string, self.config_string_img = get_conf_string(self.folder_save_run, filename)




class GridConfig(BaselineConfig):
    def __init__(self,
            n_years,
            rp1,
            rp2,
            n_doses,
            primary_inoculum=None
            # zeroth_season_reproduction=True
            ):

        super().__init__(n_years, rp1, rp2, primary_inoculum)
        
        self.strategy = 'mix'

        self.n_doses = n_doses

        self.add_string()


    def add_string(self, extra_detail=None):
        self.add_baseline_str()
        filename = f"grid/{self.save_string}_Nd={self.n_doses}_S={self.strategy}"
        
        if extra_detail is not None:
            filename = filename + extra_detail
        
        self.config_string, self.config_string_img = get_conf_string(self.folder_save_run, filename)
        

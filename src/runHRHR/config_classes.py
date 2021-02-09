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

        self.save_string = f"n_y={n_years}_rps={rp1},{rp2}_{str(primary_inoculum)[0]}_{str(zeroth_season_reproduction)[0]}"



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
        self.config_string = self.folder_save_run + filename.replace(".", ",") + ".pickle"



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
        filename = f"grid/{self.save_string}_n_d={self.n_doses},{self.strategy}"
        self.config_string = self.folder_save_run + filename.replace(".", ",") + ".pickle"
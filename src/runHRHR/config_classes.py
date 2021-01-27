class BaselineConfig:
    def __init__(self, n_years, rp1, rp2, primary_inoculum, within_season_before):

        self.load_saved = True

        self.folder_save = '../Outputs'
        
        self.n_years = n_years

        self.res_props = dict(
            f1 = rp1,
            f2 = rp2
            )

        self.primary_inoculum = primary_inoculum

        self.within_season_before = within_season_before

        self.save_string = f"{rp1}_{rp2}_plus_some_other_stuff"



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
                within_season_before=True
                ):
        
        super().__init__(n_years, rp1, rp2, primary_inoculum, within_season_before)
        
        self.fung1_doses = dict(
            spray_1 = [d11]*n_years,
            spray_2 = [d12]*n_years
            )

        self.fung2_doses = dict(
            spray_1 = [d21]*n_years,
            spray_2 = [d22]*n_years
            )

        self.config_string = f"{self.folder_save}/{self.save_string}.pickle"



class GridConfig(BaselineConfig):
    def __init__(self,
            n_years,
            rp1,
            rp2,
            n_doses,
            primary_inoculum=None,
            within_season_before=True
            ):

        super().__init__(n_years, rp1, rp2, primary_inoculum, within_season_before)
        
        self.strategy = 'mix'

        self.n_doses = n_doses

        self.config_string = f"{self.folder_save}/{self.save_string}.pickle"
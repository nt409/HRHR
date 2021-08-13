from math import ceil

# alpha: scale factor for omega (max effect)
# of fungicides after resistance (partial res type 1),
# or for curvature (partial res type 2)


class Parameters:
    def __init__(self):
        self.omega_1 = 1
        self.omega_2 = 1
        
        # partial resistance on 
        # - asymptote (PR type 1)
        # - curvature (PR type 2)
        
        # self.alpha_1 = 0
        # self.alpha_1_C = 1
        # self.alpha_2 = 0
        # self.alpha_2_C = 1
        
        self.theta_1 = 9.6
        self.theta_2 = 9.6
        
        # effect on latent period different?
        # self.omega_1_L = 1
        # self.omega_2_L = 1
        
        # self.theta_1_L = legal_dose*t1
        # self.theta_2_L = legal_dose*t2

        self.r = 1.26*10**(-2)
        self.k = 1
        
        self.S_0 = 0.05/4.2
        self.beta = 1.56*10**(-2)
        self.gamma = 1/266
        
        self.mu = 1/456
        self.nu = 8.5*10**(-3)
        
        self.delta_1 = 1.11*10**(-2)
        self.delta_2 = 1.11*10**(-2)
        
        self.init_den = 1.09*10**(-2)/4.2
        
        # inoculum
        self.inoc_frac = 8
        self.inoc_frac_integral = 0.06
        # proportion of inoculum density explicitly from last year's epidemic
        self.last_year_prop = 0
        
        self.T_emerge = 1212
        self.T_GS32 = 1456
        self.T_GS39 = 1700
        self.T_GS61 = 2066
        self.T_GS87 = 2900
        
        self.nstepz = 10**3
        self.dt = 20

        self.t_points = ceil((self.T_GS87 - self.T_emerge)/self.dt)
        
        self.yield_threshold = 95
        
        self.JSON_path = 'HR_HR/Asexual_config/JSON/'
        self.pickle_path = 'HR_HR/Asexual_output/Saved_pickles/Cluster_version/'
        
        self.no_variables = 16
        
        self.S_ind = 0
        self.ERR_ind = 1
        self.ERS_ind = 2
        self.ESR_ind = 3
        self.ESS_ind = 4
        
        # self.IRR_ind = 5
        # self.IRS_ind = 6
        # self.ISR_ind = 7
        # self.ISS_ind = 8
        
        # self.R_ind = 9
        
        # self.PRR_ind = 10
        # self.PRS_ind = 11
        # self.PSR_ind = 12
        # self.PSS_ind = 13

        self.fung_1_ind = 14
        self.fung_2_ind = 15




PARAMS = Parameters()

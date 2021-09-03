from math import ceil

# alpha: scale factor for omega (max effect)
# of fungicides after resistance (partial res type 1),
# or for curvature (partial res type 2)


class Parameters:
    def __init__(self):
        
        # default fung params
        self.omega_1 = 1
        self.omega_2 = 1
        
        self.theta_1 = 9.6
        self.theta_2 = 9.6
        

        # pathosystem params
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
        
        
        self.T_emerge = 1212
        self.T_GS32 = 1456
        self.T_GS39 = 1700
        self.T_GS61 = 2066
        self.T_GS87 = 2900
        
        self.nstepz = 10**3
        self.dt = 20

        self.t_points = ceil((self.T_GS87 - self.T_emerge)/self.dt)
        
        self.yield_threshold = 95
        
        
        self.no_variables = 16
        
        self.S_ind = 0
        self.ERR_ind = 1
        self.ERS_ind = 2
        self.ESR_ind = 3
        self.ESS_ind = 4

        self.fung_1_ind = 14
        self.fung_2_ind = 15



PARAMS = Parameters()

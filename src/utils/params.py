omega_1, omega_1_L, omega_2, omega_2_L = (1,1,1,1)
alpha_1, alpha_1_C, alpha_2, alpha_2_C = (0,1,0,1)

# Scale factor for omega (max effect) of fungicides after resistance (partial res type 1),
# or for curvature (partial res type 2)

t1 = 9.6
t2 = 9.6

legal_dose = 1
theta_1 = legal_dose*t1
theta_1_L = legal_dose*t1
theta_2 = legal_dose*t2
theta_2_L = legal_dose*t2

# Smaller omega, theta speeds transitions up (by reducing the effect of the fungicide)

#----------------------------------------------------------------------------------------------
r = 1.26*10**(-2)
k = 4.2/4.2
S_0 = 0.05/4.2

#----------------------------------------------------------------------------------------------
beta = 1.56*10**(-2)
gamma = 1/266
mu = 1/456
nu = 8.5*10**(-3)

delta_1 = 1.11*10**(-2)
delta_2 = 1.11*10**(-2)

init_den= 1.1*10**(-2)/4.2
## 1.09*10**(-2) but James' mistake, so 1.1*10**(-2)

innoc_frac = 8
innoc_frac_integral = 0.06
den_frac = 1

#----------------------------------------------------------------------------------------------
T_emerge = 1212
T_GS32 = 1456
T_GS39 = 1700
T_GS61 = 2066
T_GS87 = 2900

#----------------------------------------------------------------------------------------------
nstepz = 10**3
dt = 10
res_prop_calc_method = 'final_value'
strategy = 'mix'
sex_prop = 1
mixed_sex = False
Yield_threshold = 95




class Parameters:
    def __init__(self):
        self.omega_1 = omega_1
        self.omega_1_L = omega_1_L
        self.omega_2 = omega_2
        self.omega_2_L = omega_2_L

        self.alpha_1 = alpha_1
        self.alpha_1_C = alpha_1_C
        self.alpha_2 = alpha_2
        self.alpha_2_C = alpha_2_C

        self.theta_1 = theta_1
        self.theta_1_L = theta_1_L
        self.theta_2 = theta_2
        self.theta_2_L = theta_2_L

        self.r = r
        self.k = k
        
        self.S_0 = S_0
        self.beta = beta
        self.gamma = gamma
        
        self.mu = mu
        self.nu = nu
        
        self.delta_1 = delta_1
        self.delta_2 = delta_2
        
        self.init_den = init_den
        self.innoc_frac = innoc_frac
        self.innoc_frac_integral = innoc_frac_integral
        self.den_frac = den_frac
        
        self.T_emerge = T_emerge
        self.T_GS32 = T_GS32
        self.T_GS39 = T_GS39
        self.T_GS61 = T_GS61
        self.T_GS87 = T_GS87
        
        self.nstepz = nstepz
        self.dt = dt
        
        self.res_prop_calc_method = res_prop_calc_method
        self.strategy = strategy
        self.sex_prop = sex_prop
        self.mixed_sex = mixed_sex
        self.yield_threshold = 95
        
        self.JSON_path = 'HR_HR/Asexual_config/JSON/'
        self.pickle_path = 'HR_HR/Asexual_output/Saved_pickles/Cluster_version/'
        
        self.no_variables = 16
        self.S_ind = 0
        self.ER_ind = 1
        self.ERS_ind = 2
        self.ESR_ind = 3
        self.ES_ind = 4
        self.IR_ind = 5
        self.IRS_ind = 6
        self.ISR_ind = 7
        self.IS_ind = 8
        self.R_ind = 9
        self.PR_ind = 10
        self.PRS_ind = 11
        self.PSR_ind = 12
        self.PS_ind = 13
        self.Fung1_ind = 14
        self.Fung2_ind = 15


PARAMS = Parameters()

params_dict = {parm: getattr(PARAMS, parm, None) for parm in vars(PARAMS)}
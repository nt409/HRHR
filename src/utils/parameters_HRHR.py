omega_1, omega_1_L, omega_2, omega_2_L  = (1,1,1,1) #(0.6,0.6,0.6,0.6) # (0.7,0.7,0.7,0.7)  #originally (0.48,0,1,1)
# Make sure effect is between 0 and 1!
# So omega in (0, 1/(( 1-exp(-theta*1) ) ) ) which basically means (0,1), since with theta 9.6 get 1.0000677
alpha_1, alpha_1_C, alpha_2, alpha_2_C  = (0,1,0,1) #originally (1,1,0,1). Think first two don't actually matter-no resistance in first fungicide in James' model
#Scale factor for omega (max effect) of fungicides after resistance (partial res type 1), or for curvature (partial res type 2)
t1 = 9.6 #7 # 9.6
t2 = 9.6 #8 # 9.6
legal_dose = 1
theta_1, theta_1_L, theta_2, theta_2_L  = (legal_dose*t1, legal_dose*t1, legal_dose*t2, legal_dose*t2) # (20,20,20,20) #(9.9, 9.9, 9.6, 9.6)
# Smaller omega, theta speeds transitions up (by reducing the effect of the fungicide)
#----------------------------------------------------------------------------------------------
r                 = 1.26*10**(-2) #1.26*10**(-2)
k                 = 4.2/4.2##      after non-dimensionalisation     #4.2
S_0               = 0.05/4.2##     after non-dimensionalisation     #0.05         
#----------------------------------------------------------------------------------------------
beta              = 1.56*10**(-2)  # 1.56 #1.1#1.3          #1.56*10**(-2) = 0.0156
gamma             = 1/266                          #1/266        
mu                = 1/456                          #1/456 = 0.00219
nu                = 8.5*10**(-3)                   #8.5*10**(-3) 
delta_1, delta_2  = (1.11*10**(-2),1.11*10**(-2))#(6.91*10**(-3),6.91*10**(-3))#(1.11*10**(-2), 1.11*10**(-2))           #(6.91*10**(-3), 1.11*10**(-2))
init_den          = 1.1*10**(-2)/4.2##  after non-dimensionalisation                 #1.09*10**(-2) but James' mistake, so 1.1*10**(-2)
innoc_frac        = 8       #6.70017018514 #=1/0.14924994028036986, which means it's constant with no spray
innoc_frac_integral = 0.06
den_frac          = 1#0.2     # 0.8 #1 for james
#----------------------------------------------------------------------------------------------
T_emerge, T_GS32, T_GS39, T_GS61, T_GS87  = (1212,1456,1700,2066,2900) #(1212,1456,1700,2066,2900)
#----------------------------------------------------------------------------------------------
nstepz=10**3
dt =10
res_prop_calc_method = 'final_value'
strategy = 'mix'
sex_prop = 1
Yield_threshold = 95

JSON_path    = 'HR_HR/Asexual_config/JSON/'
pickle_path  = 'HR_HR/Asexual_output/Saved_pickles/Cluster_version/'
no_variables = 16 # S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2
S_ind   =       0
ER_ind  =       1
ERS_ind =       2
ESR_ind =       3
ES_ind  =       4
IR_ind  =       5
IRS_ind =       6
ISR_ind =       7
IS_ind  =       8
R_ind   =       9
PR_ind  =      10
PRS_ind =      11
PSR_ind =      12
PS_ind  =      13
Fung1_ind =    14
Fung2_ind =    15


# class example
# (omega_1,omega_1_L,omega_2,omega_2_L,alpha_1,alpha_1_C,alpha_2,alpha_2_C,theta_1,theta_1_L,theta_2,theta_2_L,r,k,S_0,beta,gamma,mu,nu,delta_1,delta_2,init_den,innoc_frac,innoc_frac_integral,den_frac,T_emerge,T_GS32,T_GS39,T_GS61,T_GS87,nstepz,dt,sex_prop,Yield_threshold,JSON_path,pickle_path,no_variables,S_ind,ER_ind,ERS_ind,ESR_ind,ES_ind,IR_ind,IRS_ind,ISR_ind,IS_ind,R_ind,PR_ind,PRS_ind,PSR_ind,PS_ind,Fung1_ind,Fung2_ind):
class params_class:
    __slots__ = ['omega_1','omega_1_L','omega_2','omega_2_L','alpha_1','alpha_1_C','alpha_2','alpha_2_C','theta_1','theta_1_L','theta_2','theta_2_L','r','k','S_0','beta','gamma','mu','nu','delta_1','delta_2','init_den','innoc_frac','innoc_frac_integral','den_frac','T_emerge','T_GS32','T_GS39','T_GS61','T_GS87','nstepz','dt','res_prop_calc_method','strategy','sex_prop','Yield_threshold','JSON_path','pickle_path','no_variables','S_ind','ER_ind','ERS_ind','ESR_ind','ES_ind','IR_ind','IRS_ind','ISR_ind','IS_ind','R_ind','PR_ind','PRS_ind','PSR_ind','PS_ind','Fung1_ind','Fung2_ind']
    def __init__(self,omega_1 = omega_1,omega_1_L = omega_1_L,omega_2 = omega_2,omega_2_L = omega_2_L,alpha_1 = alpha_1,alpha_1_C = alpha_1_C,alpha_2 = alpha_2,alpha_2_C = alpha_2_C,theta_1 = theta_1,theta_1_L = theta_1_L,theta_2 = theta_2,theta_2_L = theta_2_L,r = r,k = k,S_0 = S_0,beta = beta,gamma = gamma,mu = mu,nu = nu,delta_1 = delta_1,delta_2 = delta_2,init_den = init_den,innoc_frac = innoc_frac,innoc_frac_integral = innoc_frac_integral,den_frac = den_frac,T_emerge = T_emerge,T_GS32 = T_GS32,T_GS39 = T_GS39,T_GS61 = T_GS61,T_GS87 = T_GS87,nstepz = nstepz,dt = dt,res_prop_calc_method = res_prop_calc_method,strategy = strategy,sex_prop = sex_prop,Yield_threshold = Yield_threshold,JSON_path = JSON_path,pickle_path = pickle_path,no_variables = no_variables,S_ind = S_ind,ER_ind = ER_ind,ERS_ind = ERS_ind,ESR_ind = ESR_ind,ES_ind = ES_ind,IR_ind = IR_ind,IRS_ind = IRS_ind,ISR_ind = ISR_ind,IS_ind = IS_ind,R_ind = R_ind,PR_ind = PR_ind,PRS_ind = PRS_ind,PSR_ind = PSR_ind,PS_ind = PS_ind,Fung1_ind = Fung1_ind,Fung2_ind = Fung2_ind):
        self.omega_1             = omega_1
        self.omega_1_L           = omega_1_L
        self.omega_2             = omega_2
        self.omega_2_L           = omega_2_L
        self.alpha_1             = alpha_1
        self.alpha_1_C           = alpha_1_C
        self.alpha_2             = alpha_2
        self.alpha_2_C           = alpha_2_C
        self.theta_1             = theta_1
        self.theta_1_L           = theta_1_L
        self.theta_2             = theta_2
        self.theta_2_L           = theta_2_L
        self.r                   = r
        self.k                   = k
        self.S_0                 = S_0
        self.beta                = beta
        self.gamma               = gamma
        self.mu                  = mu
        self.nu                  = nu
        self.delta_1             = delta_1
        self.delta_2             = delta_2
        self.init_den            = init_den
        self.innoc_frac          = innoc_frac
        self.innoc_frac_integral = innoc_frac_integral
        self.den_frac            = den_frac
        self.T_emerge            = T_emerge
        self.T_GS32              = T_GS32
        self.T_GS39              = T_GS39
        self.T_GS61              = T_GS61
        self.T_GS87              = T_GS87
        self.nstepz              = nstepz
        self.dt                  = dt
        self.res_prop_calc_method= res_prop_calc_method
        self.strategy            = strategy
        self.sex_prop            = sex_prop
        self.Yield_threshold     = Yield_threshold
        self.JSON_path           = JSON_path
        self.pickle_path         = pickle_path
        self.no_variables        = no_variables # S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2
        self.S_ind               =      S_ind
        self.ER_ind              =      ER_ind
        self.ERS_ind             =      ERS_ind
        self.ESR_ind             =      ESR_ind
        self.ES_ind              =      ES_ind
        self.IR_ind              =      IR_ind
        self.IRS_ind             =      IRS_ind
        self.ISR_ind             =      ISR_ind
        self.IS_ind              =      IS_ind
        self.R_ind               =      R_ind
        self.PR_ind              =      PR_ind
        self.PRS_ind             =      PRS_ind
        self.PSR_ind             =      PSR_ind
        self.PS_ind              =      PS_ind
        self.Fung1_ind           =    Fung1_ind
        self.Fung2_ind           =    Fung2_ind


params = params_class()

params_dict = {s: getattr(params, s, None) for s in params.__slots__}


### other options ###

# params = {'omega_1': omega_1,
# 'omega_1_L': omega_1_L,
# 'omega_2': omega_2,
# 'omega_2_L': omega_2_L,
# 'alpha_1': alpha_1,
# 'alpha_1_C': alpha_1_C,
# 'alpha_2': alpha_2,
# 'alpha_2_C': alpha_2_C,
# 'theta_1': theta_1,
# 'theta_1_L': theta_1_L,
# 'theta_2': theta_2,
# 'theta_2_L': theta_2_L,
# 'r': r,
# 'k': k,
# 'S_0': S_0,
# 'beta': beta,
# 'gamma': gamma,
# 'mu': mu,
# 'nu': nu,
# 'delta_1': delta_1,
# 'delta_2': delta_2,
# 'init_den': init_den,
# 'innoc_frac': innoc_frac,
# 'innoc_frac_integral': innoc_frac_integral,
# 'den_frac': den_frac,
# 'T_emerge': T_emerge,
# 'T_GS32': T_GS32,
# 'T_GS39': T_GS39,
# 'T_GS61': T_GS61,
# 'T_GS87': T_GS87,
# 'nstepz': nstepz,
# 'dt': dt,
# 'sex_prop': sex_prop,
# 'Yield_threshold': Yield_threshold,
# 'JSON_path': 'HR_HR/Asexual_config/JSON/',
# 'pickle_path': 'HR_HR/Asexual_output/Saved_pickles/Cluster_version/',
# 'no_variables': 16, # S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2
# 'S_ind':         0,
# 'ER_ind':        1,
# 'ERS_ind':       2,
# 'ESR_ind':       3,
# 'ES_ind':        4,
# 'IR_ind':        5,
# 'IRS_ind':       6,
# 'ISR_ind':       7,
# 'IS_ind':        8,
# 'R_ind':         9,
# 'PR_ind':       10,
# 'PRS_ind':      11,
# 'PSR_ind':      12,
# 'PS_ind':       13,
# 'Fung1_ind':    14,
# 'Fung2_ind':    15
# }

# from typing import NamedTuple

# class params_NT(NamedTuple):
#     omega_1: float = omega_1
#     omega_1_L: float = omega_1_L
#     omega_2: float = omega_2
#     omega_2_L: float = omega_2_L
#     alpha_1: float = alpha_1
#     alpha_1_C: float = alpha_1_C
#     alpha_2: float = alpha_2
#     alpha_2_C: float = alpha_2_C
#     theta_1: float = theta_1
#     theta_1_L: float = theta_1_L
#     theta_2: float = theta_2
#     theta_2_L: float = theta_2_L
#     r: float = r
#     k: float = k
#     S_0: float = S_0
#     beta: float = beta
#     gamma: float = gamma
#     mu: float = mu
#     nu: float = nu
#     delta_1: float = delta_1
#     delta_2: float = delta_2
#     init_den: float = init_den
#     innoc_frac: float = innoc_frac
#     innoc_frac_integral: float = innoc_frac_integral
#     den_frac: float = den_frac
#     T_emerge: float = T_emerge
#     T_GS32: float = T_GS32
#     T_GS39: float = T_GS39
#     T_GS61: float = T_GS61
#     T_GS87: float = T_GS87
#     nstepz: int = nstepz
#     dt: float = dt
#     sex_prop: float = sex_prop
#     Yield_threshold: float = Yield_threshold
#     JSON_path: str = 'HR_HR/Asexual_config/JSON/'
#     pickle_path: str = 'HR_HR/Asexual_output/Saved_pickles/Cluster_version/'
#     no_variables: int = 16 # S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,Fung1,Fung2
#     S_ind: int =         0
#     ER_ind: int =        1
#     ERS_ind: int =       2
#     ESR_ind: int =       3
#     ES_ind: int =        4
#     IR_ind: int =        5
#     IRS_ind: int =       6
#     ISR_ind: int =       7
#     IS_ind: int =        8
#     R_ind: int =         9
#     PR_ind: int =       10
#     PRS_ind: int =      11
#     PSR_ind: int =      12
#     PS_ind: int =       13
#     Fung1_ind: int =    14
#     Fung2_ind: int =    15

    # def _asdict(self):
    #     d = super()._asdict()
    #     return d

# d = _asdict(params_NT)
#----------------------------------------------------------------------------------------------





# (omega_1,omega_1_L,omega_2,omega_2_L,alpha_1,alpha_1_C,alpha_2,alpha_2_C,theta_1,theta_1_L,theta_2,theta_2_L,r,k,S_0,beta,gamma,mu,nu,delta_1,delta_2,init_den,innoc_frac,innoc_frac_integral,den_frac,T_emerge,T_GS32,T_GS39,T_GS61,T_GS87,nstepz,dt,sex_prop,Yield_threshold,JSON_path,pickle_path,no_variables,S_ind,ER_ind,ERS_ind,ESR_ind,ES_ind,IR_ind,IRS_ind,ISR_ind,IS_ind,R_ind,PR_ind,PRS_ind,PSR_ind,PS_ind,Fung1_ind,Fung2_ind)
# (omega_1 = omega_1,omega_1_L = omega_1_L,omega_2 = omega_2,omega_2_L = omega_2_L,alpha_1 = alpha_1,alpha_1_C = alpha_1_C,alpha_2 = alpha_2,alpha_2_C = alpha_2_C,theta_1 = theta_1,theta_1_L = theta_1_L,theta_2 = theta_2,theta_2_L = theta_2_L,r = r,k = k,S_0 = S_0,beta = beta,gamma = gamma,mu = mu,nu = nu,delta_1 = delta_1,delta_2 = delta_2,init_den = init_den,innoc_frac = innoc_frac,innoc_frac_integral = innoc_frac_integral,den_frac = den_frac,T_emerge = T_emerge,T_GS32 = T_GS32,T_GS39 = T_GS39,T_GS61 = T_GS61,T_GS87 = T_GS87,nstepz = nstepz,dt = dt,sex_prop = sex_prop,Yield_threshold = Yield_threshold,JSON_path = JSON_path,pickle_path = pickle_path,no_variables = no_variables,S_ind = S_ind,ER_ind = ER_ind,ERS_ind = ERS_ind,ESR_ind = ESR_ind,ES_ind = ES_ind,IR_ind = IR_ind,IRS_ind = IRS_ind,ISR_ind = ISR_ind,IS_ind = IS_ind,R_ind = R_ind,PR_ind = PR_ind,PRS_ind = PRS_ind,PSR_ind = PSR_ind,PS_ind = PS_ind,Fung1_ind = Fung1_ind,Fung2_ind = Fung2_ind)
# params2 = params_class2()

#####JAMES######
# params_James = params
# params_James["omega_1"]   = 0.48
# params_James["omega_2"]   = 1
# params_James["omega_2_L"] = 1
# params_James["theta_1"]   = 9.9
# params_James["theta_1_L"] = 9.9
# params_James["r"]         = 1.26*10**(-2)
# params_James["k"]         = 4.2
# params_James["S_0"]       = 0.05

# omega_1, omega_1_L, omega_2, omega_2_L  = (0.48,0,1,1)
# alpha_1, alpha_1_C, alpha_2, alpha_2_C  = (0,1,0,1) #originally (1,1,0,1). Think first two don't actually matter-no resistance in first fungicide in James' model
# #Scale factor for omega (max effect) of fungicides after resistance (partial res type 1), or for curvature (partial res type 2)
# theta_1, theta_1_L, theta_2, theta_2_L  = (9.9, 9.9, 9.6, 9.6)
# # Smaller omega, theta speeds transitions up (by reducing the effect of the fungicide)
# #----------------------------------------------------------------------------------------------
# r                 = 1.26*10**(-2)
# k                 = 4.2               #/4.2##      after non-dimensionalisation     #4.2
# S_0               = 0.05              #/4.2##     after non-dimensionalisation     #0.05         
# #----------------------------------------------------------------------------------------------
# beta              = 1.56*10**(-2)
# gamma             = 1/266
# mu                = 1/456
# nu                = 8.5*10**(-3)
# delta_1, delta_2  = (6.91*10**(-3), 1.11*10**(-2))
# init_den          = 1.1*10**(-2)      #/4.2##  after non-dimensionalisation                 #1.09*10**(-2) but James' mistake, so 1.1*10**(-2)
# innoc_frac        = 8                               #whatever
# den_frac          = 1                               #0.8     # 0.8 #1 for james
# #----------------------------------------------------------------------------------------------
# T_emerge, T_GS32, T_GS39, T_GS61, T_GS87  = (1212, 1456,1700,2066,2900)
from model.utils import Fungicide, growth, senescence
from model.params import PARAMS
from model.ode_system import BaseSystem


class ODESystemWithinSeasonSex(BaseSystem):
    def __init__(self, fungicide_params) -> None:
        super().__init__(fungicide_params)

    



    def system(self, t, y):

        S,ER,ERS,ESR,ES,IR,IRS,ISR,IS,R,PR,PRS,PSR,PS,conc_1,conc_2 = y

        A = S + ER + ERS + ESR + ES + IR + IRS + ISR + IS + R

        dydt = [self.growth_fn(A,t)
             - (self.senescence_fn(t))*S
             -  S * (PARAMS.beta/A) * (
                  (IR + PR)
                + (IRS + PRS) * (self.fcide2.effect(conc_2))
                + (ISR + PSR) * (self.fcide1.effect(conc_1))
                + (IS  +  PS) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2))),
            
            S*(PARAMS.beta/A) * (IR + PR) - (self.senescence_fn(t)) * ER  - PARAMS.gamma * ER,
            S*(PARAMS.beta/A) * (IRS + PRS) * (self.fcide2.effect(conc_2)) - (self.senescence_fn(t)) * ERS - PARAMS.gamma * (self.fcide2.effect(conc_2)) * ERS,
            S*(PARAMS.beta/A) * (ISR + PSR) * (self.fcide1.effect(conc_1)) - (self.senescence_fn(t)) * ESR - PARAMS.gamma * (self.fcide1.effect(conc_1)) * ESR,
            S*(PARAMS.beta/A) * (IS  +  PS) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2)) - (self.senescence_fn(t)) * ES  - PARAMS.gamma * (self.fcide1.effect(conc_1))*(self.fcide2.effect(conc_2)) * ES,
            
            PARAMS.gamma * ER   -  PARAMS.mu * IR,
            PARAMS.gamma * (self.fcide2.effect(conc_2)) * ERS  -  PARAMS.mu * IRS,
            PARAMS.gamma * (self.fcide1.effect(conc_1)) * ESR  -  PARAMS.mu * ISR,
            PARAMS.gamma * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2)) * ES   -  PARAMS.mu * IS,
            
            PARAMS.mu * (IR + IRS + ISR + IS)   +  (self.senescence_fn(t)) * (S + ER + ERS + ESR + ES),
            
            - PARAMS.nu * PR,
            - PARAMS.nu * PRS,
            - PARAMS.nu * PSR,
            - PARAMS.nu * PS,
            
            - self.fcide1.delta * conc_1,
            - self.fcide2.delta * conc_2
            ]

        return dydt

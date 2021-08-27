from model.ode_system import BaseSystem


class ODESystemWithinSeasonSex(BaseSystem):
    def __init__(self, fungicide_params, sex_prop) -> None:
        super().__init__(fungicide_params)
        self.sex_prop = sex_prop


    def system(self, t, y):
        pars = self.pars
        eps = self.sex_prop

        S,ERR,ERS,ESR,ESS,IRR,IRS,ISR,ISS,R,PRR,PRS,PSR,PSS,conc_1,conc_2 = y

        A = sum([S, ERR, ERS, ESR, ESS, IRR, IRS, ISR, ISS, R])

        sum_E = sum([ERR, ERS, ESR, ESS])
        sum_I = sum([IRR, IRS, ISR, ISS])
        sum_P = sum([PRR, PRS, PSR, PSS])


        I_div = 0 if sum_I==0 else 1/sum_I
        P_div = 0 if sum_P==0 else 1/sum_P


        dydt = [self.growth_fn(A,t)
             - (self.senescence_fn(t))*S
             -  S * (pars.beta/A) * (
                  (IRR + PRR)
                + (IRS + PRS) * (self.fcide2.effect(conc_2))
                + (ISR + PSR) * (self.fcide1.effect(conc_1))
                + (ISS + PSS) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2))),
            

            S*(pars.beta/A) * ((1-eps)*(IRR + PRR)
                                + eps*(I_div*(IRS + IRR)*(ISR + IRR)
                                     + P_div*(PRS + PRR)*(PSR + PRR)))
                 - (self.senescence_fn(t)) * ERR  - pars.gamma * ERR,
            
            S*(pars.beta/A) * ((1-eps)*(IRS + PRS) 
                                + eps*(I_div*(IRS + IRR)*(IRS + ISS)
                                     + P_div*(PRS + PRR)*(PRS + PSS))
                                    ) * (self.fcide2.effect(conc_2))
                 - (self.senescence_fn(t)) * ERS - pars.gamma * (self.fcide2.effect(conc_2)) * ERS,
            
            S*(pars.beta/A) * ((1-eps)*(ISR + PSR) 
                                + eps*(I_div*(ISS + ISR)*(ISR + IRR)
                                     + P_div*(PSS + PSR)*(PSR + PRR))
                                    ) * (self.fcide1.effect(conc_1))
                 - (self.senescence_fn(t)) * ESR - pars.gamma * (self.fcide1.effect(conc_1)) * ESR,
            
            S*(pars.beta/A) * ((1-eps)*(ISS + PSS) 
                                + eps*(I_div*(ISS + ISR)*(ISS + IRS)
                                     + P_div*(PSS + PSR)*(PSS + PRS))
                                    ) * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2))
                 - (self.senescence_fn(t)) * ESS  - pars.gamma * (self.fcide1.effect(conc_1))*(self.fcide2.effect(conc_2)) * ESS,
            
            pars.gamma * ERR  - pars.mu * IRR,
            pars.gamma * (self.fcide2.effect(conc_2)) * ERS  - pars.mu * IRS,
            pars.gamma * (self.fcide1.effect(conc_1)) * ESR  - pars.mu * ISR,
            pars.gamma * (self.fcide1.effect(conc_1)) * (self.fcide2.effect(conc_2)) * ESS  - pars.mu * ISS,
            
            pars.mu * sum_I  + (self.senescence_fn(t)) * (S + sum_E),
            
            - pars.nu * PRR,
            - pars.nu * PRS,
            - pars.nu * PSR,
            - pars.nu * PSS,
            
            - self.fcide1.delta * conc_1,
            - self.fcide2.delta * conc_2
            ]

        return dydt

import matplotlib.pyplot as plt
import numpy as np
from Functions_and_plotting.functions_HRHR import master_loop_one_tactic, master_loop_grid_of_tactics
from Functions_and_plotting.plotter_HRHR    import Overlay_plotter
from Functions_and_plotting.parameters_HRHR import params
from math import ceil, log, exp
# #----------------------------------------------------------------------------------------------
run_SR_prod          = 0
run_what_to_do       = 0 #1
run_yield_full_dose2 = 0 # keeps phi constant, vary p
run_yield_full_dose3 = 0 # doesn't work, but not useful. Keeps phi constant
run_yield_full_dose4 = 0 # keeps C constant, vary phi
run_delta_int        = 0
# #----------------------------------------------------------------------------------------------
phi = 10**(-4)
doses = 0.5
n_d = 8
n_seas = 1
# #----------------------------------------------------------------------------------------------
if run_SR_prod == 1: # keeping phi constant
    n = 4
    b = np.linspace(0,0.5,n)
    r_1 = np.zeros(n)
    r_2 = np.zeros(n)
    Yield = np.zeros((n_d,n_d,n_seas,n))
    Selection_array_1 = np.zeros((n_d,n_d,n_seas+1,n))
    Selection_array_2 = np.zeros((n_d,n_d,n_seas+1,n))
    for i in range(n):
        r_1[i] = b[i]*phi
        r_2[i] = (1-b[i])*phi
        output = master_loop_grid_of_tactics(n_d,n_seas,r_1[i],r_2[i])
        Yield[:,:,:,i] = output['Yield']
        Selection_array_1[:,:,:,i] = output['Selection_array_1']
        Selection_array_2[:,:,:,i] = output['Selection_array_2']
    
    G = np.zeros((n_d,n_d))
    k = 2
    for i in range(n_d):
            for j in range(n_d):
                    G[i,j] = Selection_array_1[i,j,1,k]*Selection_array_2[i,j,1,k]

    Overlay_plotter(Col_attr={'Col_mat':G,'Col_label':'Yield','Col_bds_vec':(1,)},figure_attributes={'title':'Product SRs'})
    Overlay_plotter(Col_attr={'Col_mat': G,'Col_label': 'SR','Col_bds_vec':(1,),'Con_mats':(Yield[:,:,0,2]),'Con_levels':([params.Yield_threshold]),'Con_inline':(None)})
    Overlay_plotter(Con_attr={'Con_mats':(G,Yield[:,:,0,2]),'Con_levels':([1,2,3,3.25,3.5,4,4.5],[params.Yield_threshold]),'Con_inline':('inline','inline')},figure_attributes={'title':'G'})
    Overlay_plotter(Con_attr={'Con_mats':(Yield[:,:,0,0],Yield[:,:,0,1],Yield[:,:,0,2],Yield[:,:,0,3]),'Con_levels':([params.Yield_threshold],[params.Yield_threshold],[params.Yield_threshold],[params.Yield_threshold]),'Con_inline':('inline','inline','inline','inline')},figure_attributes={'title':'r1'})
###--------------------------------------------------------------------------------------------------------------------------
if run_yield_full_dose2 == 1: # keeping phi constant
    n_p = 5
    phi_vec = np.linspace(0.01,0.5,n_p) # [10**(-2),0.2,0.6,0.8,1]
    Yield_vec = np.zeros((n_p,len(phi_vec)))
    p = np.linspace(0.01,0.5,n_p) #[10**(-4),10**(-3),10**(-2),0.2,0.5] #np.linspace(0.01,0.5,n_p)
    for i in range(n_p):
        for j in range(len(phi_vec)):
            output = master_loop_one_tactic([doses],[doses],[doses],[doses],p[i]*phi_vec[j],(1-p[i])*phi_vec[j]) # but full dose might be 0.5
            Yield_vec[i,j] = output['Yield_vec']
    

    # print(Yield_vec)

    # print(Yield_vec[1,:]) # vary p
    # print(Yield_vec[:,1]) # vary phi
    # looks like varying phi makes less difference

    Overlay_plotter(Con_attr={'Con_mats':(Yield_vec),'Con_inline':('inline')},figure_attributes={'x_lab':'p','y_lab':'phi','xtick': p,'ytick': phi_vec,'title':'Yield, p, phi'})
    Overlay_plotter(Col_attr={'Col_mat':Yield_vec,'Col_label':'Yield','Col_bds_vec':(params.Yield_threshold,)},figure_attributes={'title':'Yield, p, phi'})

    X, Y = np.meshgrid(np.linspace(0,1,n_p),np.linspace(0,1,n_p))
    fig5 = plt.figure()
    ax5 = fig5.gca(projection='3d')
    surf5 = ax5.plot_wireframe(X,Y,Yield_vec)
###--------------------------------------------------------------------------------------------------------------------------
if run_yield_full_dose3 == 1: # keeping phi constant
    n_p = 5
    phi_vec = np.linspace(0.01,0.5,n_p)
    C_vec = [10**(-2),0.2,0.6,0.8,1]
    d1 = 0.1
    d2 = 0.1
    d12 = 0.01
    rho = 10
    A = exp(rho*d12)
    B = A + 0.5*(exp(rho*d1) + exp(rho*d2))
    C = exp(rho*d1) - exp(rho*d2)
    D = A - (exp(rho*d1) + exp(rho*d2)) + exp(rho)
    Yield_vec = np.zeros((n_p,len(phi_vec)))
    p = np.linspace(0.01,0.5,n_p) #[10**(-4),10**(-3),10**(-2),0.2,0.5] #np.linspace(0.01,0.5,n_p)
    for i in range(n_p):
        for j in range(len(phi_vec)):
            pq = p[i]*(1-p[i])
            r = p[i] - 0.5
            phi_vec = np.roots([pq*D,B+C*r,A-C_vec[j]])[1]
            output = master_loop_one_tactic([doses],[doses],[doses],[doses],p[i]*phi_vec[j],(1-p[i])*phi_vec[j]) # but full dose might be 0.5
            Yield_vec[i,j] = output['Yield_vec']
###--------------------------------------------------------------------------------------------------------------------------
if run_yield_full_dose4 == 1: # keeping pq phi^2 constant, plotting for phi on x axis
    n_p = 15
    ##----------------------------------------------------------
    C_vec = np.concatenate(([10**(-8),10**(-6)],np.linspace(10**(-4),0.12,n_p-2)))
    a = [10**(-6),10**(-3),4*10**(-3),8*10**(-3),1.2*10**(-2),1.6*10**(-2)]
    b = np.linspace(2*10**(-2),0.99,n_p-6)
    phi = np.concatenate((a,b))
    ##----------------------------------------------------------
    Yield_vec = np.zeros((n_p,len(C_vec)))
    # Yield = np.zeros((n_d,n_d,n_seas,n_p,len(C_vec)))
    p_here = np.zeros((n_p,len(C_vec)))
    ##----------------------------------------------------------
    for i in range(n_p):
        for j in range(len(C_vec)):
            p_here[i,j] = np.roots([1,-1,C_vec[j]/((phi[i])**2)])[0]
            if C_vec[j]/((phi[i])**2)<0.25 and p_here[i,j] < 1 and p_here[i,j] > 0:
                output = master_loop_one_tactic([doses],[doses],[doses],[doses],p_here[i,j]*phi[i],(1-p_here[i,j])*phi[i]) # but full dose might be 0.5
                Yield_vec[i,j] = output['Yield_vec']
            else:
                Yield_vec[i,j] =93

    Overlay_plotter(Con_attr={'Con_mats':(Yield_vec),'Con_inline':('inline')},figure_attributes={'x_lab':'phi','y_lab':'C','xtick': phi,'ytick': C_vec,'title':'Yield, phi, C'})
    Overlay_plotter(Col_attr={'Col_mat':Yield_vec,'Col_label':'Yield','Col_bds_vec':(params.Yield_threshold,)},figure_attributes={'title':'Yield, phi, C'})
    

    Overlay_plotter(Con_attr={'Con_mats':(p_here),'Con_levels':([0.5,0.6,0.7,0.8,0.9,1]),'Con_inline':('inline')},figure_attributes={'x_lab':'phi','y_lab':'C','xtick': phi,'ytick': C_vec,'title':'p, phi, C'})

###--------------------------------------------------------------------------------------------------------------------------
if run_delta_int == 1:
    omega = 1
    C_s = 0.1

    beta = 1.56*10**(-2)
    S = 1

    Delta = 10**(-2) 

    C_init = 1 # full dose

    t_2 = params.T_GS61 - params.T_GS32
    D_2 = (omega/Delta)*log((C_s + C_init*exp(-Delta*t_2))/(C_s + C_init))
    D = (t_2 + D_2)/t_2
    rho = beta*S*t_2


    print(exp(rho))
    print(exp(rho*D))

    print(exp(rho*D)/exp(rho))# print(exp(rho*D_2/t_2))

    print(D)
    # print(t_2)
    # print(D_2)
plt.show()
from scipy.integrate import odeint, simps, ode
import numpy as np
import matplotlib.pyplot as plt
from math import exp, ceil, floor

# can't change delta within script??

dt = 0.01
S_10 = 1 #7
phi = 10**(-4)
den = 0.1
IR_0  =  (phi**2)*den
ISR_0 =  phi*(1-phi)*den
IRS_0 =  phi*(1-phi)*den
IS_0  =  ((1-phi)**2)*den
R_0 = 0
no_variables_2 = 8
t0 = 0
t1 = 1000
t_app_start = 0
t_app_end = 1000
delta = 0.01

def solv_expo_model(t,y):
    S,IR,IRS,ISR,IS,R,d1,d2 = y
    if t > t_app_start and t < t_app_end:
        dose_response = delta
        dose_response = delta
    else:
        dose_response = 1
        dose_response = 1
    dydt = [ - 0*S * beta * (  
               (1                          ) * IR
             + (dose_response              ) * IRS
             + (dose_response              ) * ISR
             + (dose_response*dose_response) * IS  ),
        (S*beta * (1                          ) -  mu) * IR,
        (S*beta * (dose_response              ) -  mu) * IRS,
        (S*beta * (dose_response              ) -  mu) * ISR,
        (S*beta * (dose_response*dose_response) -  mu) * IS,
        mu * (IR + IRS +  ISR +  IS),
        0,
        0]
    return dydt

def integrator_expo():
    y0   = [S_10, IR_0, IRS_0,  ISR_0, IS_0, R_0, delta, delta]
    sol  = ode(solv_expo_model,jac=None).set_integrator('dopri5',nsteps= nstepz)
    n1= 1 + (t1-t0)/dt
    c1 = ceil(n1-0.5)
    tim1 = np.linspace(t0,t1,c1)
    yy1  = np.zeros((no_variables_2,len(tim1)))
    i1=0
#----------------------------------------------------------------------------------------------
    sol.set_initial_value(y0,t0)
    for t in tim1[1:]:
        if sol.successful():
            yy1[:,i1] = sol.y
            i1=i1+1
            sol.integrate(t)
        else:
            raise RuntimeError('ode solver unsuccessful')
    #----------------------------------------------------------------------------------------------
    return yy1, tim1


yy1,tim1 = integrator_expo()
yy2,tim2 = integrator_expo()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(tim1,yy1[0,:],label='S')
ax.plot(tim1,yy1[1,:],label='IR')
ax.plot(tim1,yy1[2,:],label='IRS')
# ax.plot(tim1,yy1[3,:],label='ISR')
ax.plot(tim1,yy1[4,:],label='IS')
ax.plot(tim1,yy1[5,:],label='R')
ax.set_title('Disease amounts')
fig.legend()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
# ax2.plot(tim1,yy1[0,:],label='S')
# ax2.plot(tim1,yy1[1,:],label='IR')
ax2.plot(tim1,yy1[2,:],label='IRS')
ax2.plot(tim1,yy1[4,:],label='IS')
# ax2.plot(tim1,yy1[5,:],label='R')
ax2.set_title('Disease amounts')
fig2.legend()

y1 = np.zeros((no_variables_2,len(tim1)))
y2 = np.zeros((no_variables_2,len(tim1)))
y4 = np.zeros((no_variables_2,len(tim2)))
for i in range(len(tim1)):
    y1[:,i] = solv_expo_model(tim1[i],yy1[:,i])
    for j in range(no_variables_2):
        y2[j,i] = solv_expo_model(tim1[i],yy1[:,i])[j]/yy1[j,i]
        y4[j,i] = solv_expo_model(tim2[i],yy2[:,i])[j]/yy2[j,i]

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(tim1,y1[1,:],label='IR')
ax3.plot(tim1,y1[2,:],label='IRS')
ax3.plot(tim1,y1[4,:],label='IS')
ax3.set_title('Scaled growth rate')
fig3.legend()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.plot(tim1,y2[1,:],label='IR, d=0.5')
ax4.plot(tim1,y2[2,:],label='IRS, d=0.5')
ax4.plot(tim1,y2[4,:],label='IS, d=0.5')
ax4.plot(tim2,y4[1,:],label='IR, d=0.99')
ax4.plot(tim2,y4[2,:],label='IRS, d=0.99')
ax4.plot(tim2,y4[4,:],label='IS, d=0.99')
ax4.set_title('Per capita growth rate')
fig4.legend()


fig5 = plt.figure()
ax5 = fig5.add_subplot(111)
# ax5.plot(tim1,y2[1,:],label='IR')
# ax5.plot(tim1,y2[2,:],label='IRS')
ax5.plot(tim1,y2[2,:] - y2[4,:],label='IRS - IS, delta = 0.5')
ax5.plot(tim2,y4[2,:] - y4[4,:],label='IRS - IS, delta = 0.01')
ax5.set_title('Per capita growth rate difference')
fig5.legend()

fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
# ax6.plot(tim1,y2[1,:],label='IR')
# ax6.plot(tim1,y2[2,:],label='IRS')
ax6.plot(tim1,y2[1,:] - y2[2,:],label='IR - IRS, delta = 0.5')
ax6.plot(tim2,y4[1,:] - y4[2,:],label='IR - IRS, delta = 0.01')
ax6.set_title('Per capita growth rate difference')
fig6.legend()

plt.show()
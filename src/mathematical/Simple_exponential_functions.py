from scipy.integrate import odeint, simps, ode
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D
from math import exp, ceil, floor, log

from utils.parameters_HRHR import params
from utils.plotter_HRHR    import Contour_plotter, Colour_plotter, Overlay_plotter, Equal_dose, plot_disease_dynamics, bar_plotter

phi = 4*10**(-1)
p = 0.5
rho = 5

n= 20

a_1  = np.zeros((n+1,n+1))
a_2  = np.zeros((n+1,n+1))
f  = np.zeros((n+1,n+1))
G  = np.zeros((n+1,n+1))

A  = np.zeros((n+1,n+1))
B  = np.zeros((n+1,n+1))
C  = np.zeros((n+1,n+1))  
A_y= np.zeros((n+1,n+1))
B_y= np.zeros((n+1,n+1))
C_y= np.zeros((n+1,n+1))  
H_x= np.zeros((n+1,n+1))
H_y= np.zeros((n+1,n+1))
Y_line= np.zeros((n+1,n+1))
Y_line1= np.zeros((n+1,n+1))
Y_line2= np.zeros((n+1,n+1))
Y_line3= np.zeros((n+1,n+1))
Y_line4= np.zeros((n+1,n+1))

for i in range(n+1):
    for j in range (n+1):
        x = i/n
        y = j/n
        Y_line[i,j] = + exp(rho*x*y)
        Y_line1[i,j] = - phi*exp(rho*x*y) + (1-p)*phi* exp(rho*x) + p*phi* exp(rho*y) + exp(rho*x*y)
        Y_line2[i,j] = p*(1-p)*(phi**2)*(exp(rho) + exp(rho*x*y) - exp(rho*x) -  exp(rho*y)) - phi*exp(rho*x*y) + (1-p)*phi* exp(rho*x) + p*phi* exp(rho*y) + exp(rho*x*y)
        Y_line3[i,j] = 1 + p*(1-p)*(phi**2)*(exp(rho) - 1)
        Y_line4[i,j] = (1-p*phi)*(1-(1-p)*phi)*exp(rho*x*y) + p*phi*(1-(1-p)*phi)*exp(rho*x) + (1-p)*phi*(1-p*phi)*exp(rho*y) + p*(1-p)*(phi**2)*exp(rho)
        a_1[i,j] = log(p*phi*exp(rho) + (1-p*phi)*exp(rho*x))
        a_2[i,j] = log((1-p)*phi*exp(rho) + (1-(1-p)*phi)*exp(rho*y))
        f[i,j] =  log( (1-(1-p)*phi)*(1-p*phi)*exp(rho*x*y) +  (1-p)*phi*(1-p*phi)*exp(rho*x) +  p*phi*(1-(1-p)*phi)*exp(rho*y) + p*(1-p)*(phi**2)*exp(rho)  )
        A[i,j] = p*(1-p)*(1- exp(-rho*x*y)*(exp(rho) + exp(rho*y) - exp(rho*x)) + 2*y*( exp(rho*(1-x)) - 1 )   )
        B[i,j] = p * (exp(- rho*x*y) *(exp(rho*y) + exp(rho*x)) - 2 * y * exp(rho*(1-x))) - exp(rho*x*(1-y)) + 2*y - 1
        C[i,j] = 1 - 2*y
        A_y[i,j] = p*(1-p)*(1- exp(-rho*x*y)*(exp(rho) + exp(rho*x) - exp(rho*y)) + 2*x*( exp(rho*(1-y)) - 1 )   )
        B_y[i,j] = (1-p)* (exp(- rho*x*y) *(exp(rho*x) + exp(rho*y)) - 2 * x * exp(rho*(1-y))) - exp(rho*y*(1-x)) + 2*x - 1
        C_y[i,j] = 1 - 2*x
        H_x[i,j] = (phi**2)*A[i,j] + phi*B[i,j] + C[i,j]
        H_y[i,j] = (phi**2)*A_y[i,j] + phi*B_y[i,j] + C_y[i,j]

G = a_1 + a_2 - 2*f

# print(A)
# print(B)
# print(C)
# print(H_x)
# print(H_x[:,0])
print(Y_line3[:,0])
print(Y_line4[:,0])

X, Y = np.meshgrid(np.linspace(0,1,n+1),np.linspace(0,1,n+1))

fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
surf1 = ax1.plot_wireframe(X,Y,A)

fig2 = plt.figure()
ax2 = fig2.gca(projection='3d')
surf2 = ax2.plot_wireframe(X,Y,B)

fig3 = plt.figure()
ax3 = fig3.gca(projection='3d')
surf3 = ax3.plot_wireframe(X,Y,C)

fig4 = plt.figure()
ax4 = fig4.gca(projection='3d')
surf4 = ax4.plot_wireframe(X,Y,H_x)

Contour_plotter((H_x,H_y),n+1,([-0.4,-0.3,-0.2,-0.1,0,0.1],[-0.1,0]),('inline','inline'))
Contour_plotter(G,n+1,(),('inline'))
Contour_plotter(Y_line,n+1, ([1,2,3,4,5,10,15,20,40,60,80]),('inline'))
Contour_plotter(Y_line1,n+1,([1,2,3,4,5,10,15,20,40,60,80]),('inline'))
Contour_plotter(Y_line2,n+1,([1,2,3,4,5,10,15,20,40,60,80]),('inline'))
Overlay_plotter(G,n+1,'G',(0,),(G,Y_line,Y_line1,Y_line2),n+1,(),('inline','N','N','N'),0.5)
# Overlay_plotter(LTY,n_doses,'Yield',(0,),(Z,FYY),n_doses,(Z_levels_sparse,Y_levels),('N','inline'),0.5,title='LTY, Mixture')
Colour_plotter(G,n+1,'G',(1,),'Dose')
Contour_plotter((Y_line4),n+1,([7,8,9]),('inline'))

plt.show()
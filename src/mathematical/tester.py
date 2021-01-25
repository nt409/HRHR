from math import exp
import numpy as np
import matplotlib.pyplot as plt
from Functions_and_plotting.functions_HRHR import fngcide


p = 50

bet = 1.56*10**(-2)
S = 1
T = 2066 - 1456 # 366

rho = bet*S*T
dose = 1
delta = fngcide(1,9.6,1)

thresh = (exp(rho/2) -1 )/(exp(rho) -1 ) 
thresh2 = (1 - delta* exp(rho/2))* (exp(rho/2) -1 )/(exp(rho) -1 ) 
print(rho)
print(thresh)
print(delta)
print(thresh2)

x,y,z = [np.zeros(10) for i in range(3)]
print(x,y,z)
x[0] = 10
y[1] = 10
z[2] = 10
print(x,y,z)

def phi(d):

    a = exp(p*(1+d-d**2))
    b = 2*d*exp(p)
    c = (1-2*d) * exp(p*d)

    A = -a + b + c
    B = -b -2*c
    C = c

    phi_plus = (b + 2*c + (b**2 + 4*a*c)**0.5)/(2*(b+c-a))
    phi_minus = (b + 2*c - (b**2 + 4*a*c)**0.5)/(2*(b+c-a))
    return phi_minus, phi_plus, A, B, C

n=50
for i in np.linspace(0,1,n):
    phi_minus, phi_plus, A, B, C = phi(i)
    print(np.roots([A,B,C]))

# n= 11
y_minus = np.zeros(n)
y_plus = np.zeros(n)
A = np.zeros(n)
B = np.zeros(n)
C = np.zeros(n)
quad = np.zeros((n,n))
for i in range(n):
    j=(i+0.001)/(n)
    y_minus[i], y_plus[i], A[i], B[i], C[i] = phi(j)
    for k in range(n):
        l = k/n
        quad[i,k] = A[i]*(l**2) + B[i]*l + C[i]


fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(np.linspace(0.001,(n-0.999)/n,n),y_minus,label='minus')
ax.plot(np.linspace(0.001,(n-0.999)/n,n),y_plus,label='plus')
ax.grid()
fig.legend()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for i in range(n):
    ax2.plot(np.linspace(0,1,n),quad[i,:],label='dose ')
ax2.grid()

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.plot(np.linspace(0.001,(n-0.999)/n,n),A,label='A')
ax3.plot(np.linspace(0.001,(n-0.999)/n,n),B,label='B')
ax3.plot(np.linspace(0.001,(n-0.999)/n,n),C,label='C')
ax3.grid()
fig3.legend()

plt.show()
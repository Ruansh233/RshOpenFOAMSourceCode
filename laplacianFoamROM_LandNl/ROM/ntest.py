"""
Program to solve system of ODEs through multivariate Newton-Raphson method
System of ODEs:
	y1' = -0.04*y1 + 10^4*y2*y3
	y2' = 0.04*y1 - 10^4*y2*y3 - 3*10^7*y2*y2
	y3' = 3*10^7*y2*y2


Program by: Shrey Shah
Date: 30-Sep-19
"""

import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

# Defining the root-solving functions
def f1(y1_n,y2_n,y3_n,y1_o,h):
    return y1_n + 0.04*h*y1_n - 1e4*h*y2_n*y3_n - y1_o

def f2(y1_n,y2_n,y3_n,y2_o,h):
    return y2_n - 0.04*h*y1_n + 1e4*h*y2_n*y3_n + 3e7*h*pow(y2_n,2) - y2_o

def f3(y1_n,y2_n,y3_n,y3_o,h):
    return y3_n - 3e7*h*pow(y2_n,2) - y3_o

# Calculating numerical Jacobian through spatial forward differencing
def jac(Y_n,Y_o,h):
    dy = 1e-8
    y1_n = Y_n[0]
    y2_n = Y_n[1]
    y3_n = Y_n[2]

    y1_o = Y_o[0]
    y2_o = Y_o[1]
    y3_o = Y_o[2]

    # Jacobian matrix for 3 equations and 3 unknowns
    J = np.zeros((3,3))

    J[0,0] = (f1(y1_n+dy,y2_n,y3_n,y1_o,h) - f1(y1_n,y2_n,y3_n,y1_o,h))/dy
    J[0,1] = (f1(y1_n,y2_n+dy,y3_n,y1_o,h) - f1(y1_n,y2_n,y3_n,y1_o,h))/dy
    J[0,2] = (f1(y1_n,y2_n,y3_n+dy,y1_o,h) - f1(y1_n,y2_n,y3_n,y1_o,h))/dy

    J[1,0] = (f2(y1_n+dy,y2_n,y3_n,y2_o,h) - f2(y1_n,y2_n,y3_n,y2_o,h))/dy
    J[1,1] = (f2(y1_n,y2_n+dy,y3_n,y2_o,h) - f2(y1_n,y2_n,y3_n,y2_o,h))/dy
    J[1,2] = (f2(y1_n,y2_n,y3_n+dy,y2_o,h) - f2(y1_n,y2_n,y3_n,y2_o,h))/dy

    J[2,0] = (f3(y1_n+dy,y2_n,y3_n,y3_o,h) - f3(y1_n,y2_n,y3_n,y3_o,h))/dy
    J[2,1] = (f3(y1_n,y2_n+dy,y3_n,y3_o,h) - f3(y1_n,y2_n,y3_n,y3_o,h))/dy
    J[2,2] = (f3(y1_n,y2_n,y3_n+dy,y3_o,h) - f3(y1_n,y2_n,y3_n,y3_o,h))/dy

    return J

# Defining the initial conditions
Y_o = np.zeros((3,1))
Y_o[0] = 1

Y_n = np.zeros((3,1))
F = np.copy(Y_o)

# Defining the Newton-Raphson solver parameters
err = 1e9
alpha = 1
tol = 1e-12
count = 0
t = np.arange(0,0.5,step=0.1)
h = t[1] - t[0]
y1 = [1]*len(t)
y2 = [0]*len(t)
y3 = [0]*len(t)
e = [1e-6]*len(t)

# Outer time loop
for k in range(1,len(t)):
    y1_o = Y_o[0]
    y2_o = Y_o[1]
    y3_o = Y_o[2]
    Y_g = Y_o

	# Inner iterative solver loop
    while err >= tol:
        J = jac(Y_n,Y_o,h)
        # print(J)
        
        y1_g = Y_g[0]
        y2_g = Y_g[1]
        y3_g = Y_g[2]

        F[0] = f1(y1_g,y2_g,y3_g,y1_o,h)
        F[1] = f2(y1_g,y2_g,y3_g,y2_o,h)
        F[2] = f3(y1_g,y2_g,y3_g,y3_o,h)

        Y_n = Y_g - alpha*np.matmul(inv(J),F)
        err = max(abs(Y_n - Y_g))

        # Updating the guess values for a new iteration
        Y_g = Y_n
        count += 1
	
    # Updating the new time-step values
    e[k] = err
    Y_o = Y_n
    y1[k] = Y_n[0]
    y2[k] = Y_n[1]
    y3[k] = Y_n[2]

    # Prompt message for each time-step
    log_message = 'Time = {0} sec, y1 = {1}, y2 = {2}, y3 = {3}'.format\
            (round(t[k],1),round(Y_n[0][0],3),round(Y_n[1][0],3),round(Y_n[2][0],3))
    print(log_message)
    # Resetting the error criteria
    count = 1
    err = 1e9
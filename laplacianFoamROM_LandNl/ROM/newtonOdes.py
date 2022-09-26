"""
Program to solve nonlinear system of ODEs through multivariate Newton-Raphson method
System of ODEs:
	At = k1 * A.T * B * A + k2 * C * A 

Program by: Shenhui Ruan. Modified from code of Shrey Shah
Date: 26-Sep-22
"""

def newtonOdes(Bijk, Cij, A0, dt, tTotal, nEqu):
	import numpy as np
	# from numpy.linalg import inv
	# # import matplotlib.pyplot as plt

	# Defining the root-solving functions, derivate based on Backward Euler method
	def fn(n, An, An_1):
		return An[n] - An_1[n] - dt * (An.dot(Bijk[n]).dot(An) + Cij[n].dot(An))

	# Calculating numerical Jacobian through spatial forward differencing
	def jac(Ak, An_1):
		da = 1e-8
		# Jacobian matrix for nEqu equations and nEqu unknowns
		J = np.zeros((nEqu, nEqu))

		for i in range(0, nEqu):
			for j in range(0, nEqu):
				dA = np.zeros((nEqu))
				dA[j] = da

				J[i, j] = (fn(i, Ak + dA, An_1) - fn(i, Ak, An_1))/da
		return J

	def func(Ak, An_1):
		# Function matrix for nEqu equations
		F = np.zeros((nEqu))

		for i in range(0, nEqu):
			F[i] = fn(i, Ak, An_1)
		return F

	# The following is the main code	
	# Defining the initial conditions
	An = A0

	# Defining the Newton-Raphson solver parameters
	err = 1
	alpha = 1
	tol = 1e-8
	maxIter = 50
	count = 0
	t = np.arange(0, tTotal, step = dt)

	# Outer time loop
	for n in range(1, len(t)):
		Ank_1 = An
		Ank = An
		
		# Inner iterative solver loop
		while err >= tol and count <= maxIter:
			J = jac(Ank, An)
			F = func(Ank, An)

			Ank = Ank_1 - alpha*np.matmul(np.linalg.inv(J),F)
			err = max(abs(Ank - Ank_1))

			# Updating the guess values for a new iteration
			Ank_1 = Ank
			count += 1
		
		# Updating the new time-step values
		An = Ank
		print(An)

		# Resetting the error criteria
		count = 1
		err = 1e9

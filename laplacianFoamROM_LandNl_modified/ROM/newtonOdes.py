"""
Program to solve nonlinear system of ODEs through multivariate Newton-Raphson method
System of ODEs:
	At = k1 * A.T * B * A + k2 * C * A 

Program by: Shenhui Ruan. Modified from code of Shrey Shah
Date: 26-Sep-22
"""
from ast import If
import numpy as np

def newtonOdesFunc(Bijk, Cij, A0, dt, tTotal, nEqu):

	# print(Bijk)
	# print("\n________________\n")
	# print(Cij)
	# print("\n________________\n")
	# print(A0)

	# Defining the root-solving functions, derivate based on Backward Euler method
	def fn(n, An, An_1):
		return An[n] - An_1[n] - dt * (1.0e-7 * An.dot(Bijk[n]).dot(An) + 1.0e-3 * Cij[n].dot(An))

	# Calculating numerical Jacobian through spatial forward differencing
	def jac(Ak, An_1):
		da = 0.01
		# Jacobian matrix for nEqu equations and nEqu unknowns
		J = np.zeros((nEqu, nEqu))

		# set every elements of Jacobian, i is row, j is column
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
	Ank = An
	AnStore = np.empty((0, nEqu))

	# Defining the Newton-Raphson solver parameters
	err = 1
	alpha = 1
	tol = 1e-10
	maxIter = 15
	count = 1
	t = np.arange(0, tTotal, step = dt)

	# Outer time loop
	for n in range(1, len(t) + 2):
		
		Ank_1 = An
		if n%(0.1/dt) == 0:
			AnStore = np.append(AnStore, An.reshape(1, nEqu), 0)
		
		# Inner iterative solver loop
		while err >= tol and count <= maxIter:
			J = jac(Ank, An)
			F = func(Ank, An)

			Ank = Ank_1 - alpha*np.matmul(np.linalg.inv(J),F)
			err = max(abs(Ank - Ank_1))/min(abs(Ank))

			# Updating the guess values for a new iteration
			Ank_1 = Ank
			count += 1
		
		# Updating the new time-step values
		An = Ank
		# print("\n____________________________________\n")

		# print(count)

		# Resetting the error criteria
		count = 1
		err = 1e9

	return AnStore	
	

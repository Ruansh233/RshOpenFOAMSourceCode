from os import PRIO_USER
from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"

coeff = filePath + "/coeffMatrix"
spatialmode = filePath + "/modeMatrix"
laplacianMode = filePath + "/laplacianModesMatrix"
diffucoefficient = filePath + "/diffuTermCoeffMatrix"
nonLinearCoeff = filePath + "/nonLinearCoeffMatrix"

coeff_calculate = filePath + "/coeff_calculate"
snapshots_calculate = filePath + "/snapshots_calculate"

diffuTermCoeffMatrix = np.loadtxt(diffucoefficient)
nonLinearCoeffMatrix = np.loadtxt(nonLinearCoeff)

# mode numbers 
modesNum = nonLinearCoeffMatrix.shape[1]
# nonliear tensor
nonLinearCoeffTensor = nonLinearCoeffMatrix.reshape(modesNum, modesNum, modesNum)

modeMatrix = np.loadtxt(spatialmode)
modeLapMatrix = np.loadtxt(laplacianMode)

coeffMatrix = np.loadtxt(coeff)
initialA = coeffMatrix[0, :]

# print(diffuTermCoeffMatrix)
print(initialA)
# print(nonLinearCoeffTensor)

a = np.array([1, 2, 3])
# print(a)
b = a.reshape(3, 1)

c = a.dot(nonLinearCoeffTensor).dot(a)

print(c)

# d = np.array([[1], [2], [2]])
# print(a.dot(d))

# time = np.linspace(0, 50, 501)

# def odefun(t, a):
#     da = 2.5e-3 * diffuTermCoeffMatrix.dot(a)
#     return da

# sol = solve_ivp(odefun, [0, 50], initialA, dense_output=True)


# snapshots_cal = modeMatrix @ sol.sol(time)

# np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')
# np.savetxt(f"{snapshots_calculate}", snapshots_cal, fmt='%.6e', delimiter=',')

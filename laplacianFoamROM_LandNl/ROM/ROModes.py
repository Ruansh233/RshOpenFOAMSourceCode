from os import PRIO_USER
from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"

coeff = filePath + "/coeffMatrix"
spatialmode = filePath + "/modeMatrix"
laplacianMode = filePath + "/laplacianModesMatrix"
diffucoefficient = filePath + "/diffuTermCoeffMatrix"
nonLinearCoeff = filePath + "/nonLinearCoeffM"

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

print(initialA)
print(diffuTermCoeffMatrix.shape)
print(nonLinearCoeffTensor.shape)

time = np.linspace(0, 100, 1001)

ddt1 = 1.0e-7
ddt2 = 5.0e-3

def odefun(t, a):
    # da = 1.0e-7 * a.dot(nonLinearCoeffTensor).dot(a) + 2.0e-3 * diffuTermCoeffMatrix.dot(a)
    da = ddt1 * a.dot(nonLinearCoeffTensor).dot(a) + ddt2 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 100], initialA, method='Radau', dense_output=True)


snapshots_cal = modeMatrix @ sol.sol(time)

np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')
np.savetxt(f"{snapshots_calculate}", snapshots_cal, fmt='%.6e', delimiter=',')

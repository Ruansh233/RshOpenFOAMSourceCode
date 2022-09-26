from os import PRIO_USER
from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest_2/SVD"

ddt1 = 1.0e-7
ddt2 = 1.0e-3
modesNum = 2

coeff = filePath + "/coeffMatrix"
spatialmode = filePath + "/modeMatrix"
laplacianMode = filePath + "/laplacianModesMatrix"
diffucoefficient = filePath + "/diffuTermCoeffMatrix"
nonLinearCoeff = filePath + "/nonLinearCoeffM"

coeff_calculate = filePath + "/coeff_calculate"
snapshots_calculate = filePath + "/snapshots_calculate"

modeMatrix = np.loadtxt(spatialmode)[: , 0: modesNum]
modeLapMatrix = np.loadtxt(laplacianMode)[: , 0: modesNum]
coeffMatrix = np.loadtxt(coeff)[: , 0: modesNum]

diffuTermCoeffMatrix = np.loadtxt(diffucoefficient)[0: modesNum, 0: modesNum]
nonLinearCoeffMatrix = np.loadtxt(nonLinearCoeff)
nonLinearCoeffLen = nonLinearCoeffMatrix.shape[1]
nonLinearCoeffMatrix = nonLinearCoeffMatrix[: , 0: modesNum]

diffuTermCoeffMatrix = diffuTermCoeffMatrix
nonLinearCoeffMatrix = nonLinearCoeffMatrix
nonLinearCoeffTensor = np.empty((0, modesNum))

for i in range(0, modesNum):
    for j in range(0, modesNum):
        rowN = i * nonLinearCoeffLen + j
        nonLinearCoeffTensor = np.append(nonLinearCoeffTensor, nonLinearCoeffMatrix[rowN, :].reshape(1, modesNum), 0)

# nonliear tensor
nonLinearCoeffTensor = nonLinearCoeffTensor.reshape(modesNum, modesNum, modesNum)

# print(diffuTermCoeffMatrix)
# print(nonLinearCoeffTensor)

initialA = coeffMatrix[0, 0: modesNum]

print(initialA)

time = np.linspace(0, 100, 1001)

def odefun(t, a):
    da = ddt1 * a.dot(nonLinearCoeffTensor).dot(a) + ddt2 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 100], initialA, method='Radau', dense_output=True)


snapshots_cal = modeMatrix @ sol.sol(time)

np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')
np.savetxt(f"{snapshots_calculate}", snapshots_cal, fmt='%.6e', delimiter=',')

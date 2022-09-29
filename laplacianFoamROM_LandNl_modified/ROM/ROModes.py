from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"

ddt1 = 1.0e-6
ddt2 = 1.0e-3
modesNum = 5
fileLen = 20

coeff = filePath + "/coeffMatrix"
spatialmode = filePath + "/modeMatrix"
laplacianMode = filePath + "/laplacianModesMatrix"
diffucoefficient = filePath + "/diffuTermCoeffMatrix"
nonLinearCoeff1 = filePath + "/nonLinearCoeffM1"
nonLinearCoeff2 = filePath + "/nonLinearCoeffM2"

coeff_calculate = filePath + "/coeff_calculate"
snapshots_calculate = filePath + "/snapshots_calculate"

modeMatrix = np.loadtxt(spatialmode)[: , 0: modesNum]
modeLapMatrix = np.loadtxt(laplacianMode)[: , 0: modesNum]
coeffMatrix = np.loadtxt(coeff)[: , 0: modesNum]

diffuTermCoeffMatrix = np.loadtxt(diffucoefficient)[0: modesNum, 0: modesNum]
nonLinearCoeffMatrix1 = np.loadtxt(nonLinearCoeff1)[: , 0: modesNum]
nonLinearCoeffMatrix2 = np.loadtxt(nonLinearCoeff2)[: , 0: modesNum]

nonLinearCoeffTensor1 = np.empty((0, modesNum))
nonLinearCoeffTensor2 = np.empty((0, modesNum))

for i in range(0, modesNum):
    for j in range(0, modesNum):
        rowN = i * fileLen + j
        nonLinearCoeffTensor1 = np.append(nonLinearCoeffTensor1, nonLinearCoeffMatrix1[rowN, :].reshape(1, modesNum), 0)
        nonLinearCoeffTensor2 = np.append(nonLinearCoeffTensor2, nonLinearCoeffMatrix2[rowN, :].reshape(1, modesNum), 0)

# nonliear tensor
nonLinearCoeffTensor1 = nonLinearCoeffTensor1.reshape(modesNum, modesNum, modesNum)
nonLinearCoeffTensor2 = nonLinearCoeffTensor2.reshape(modesNum, modesNum, modesNum)

# print(diffuTermCoeffMatrix)
# print(nonLinearCoeffTensor)

initialA = coeffMatrix[0, 0: modesNum]

print(initialA)

time = np.linspace(0, 100, 1001)

def odefun(t, a):
    da = ddt1 * a.dot(nonLinearCoeffTensor1).dot(a) + ddt1 * a.dot(nonLinearCoeffTensor2).dot(a) + ddt2 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 100], initialA, method='Radau', dense_output=True)


snapshots_cal = modeMatrix @ sol.sol(time)

np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')
np.savetxt(f"{snapshots_calculate}", snapshots_cal, fmt='%.6e', delimiter=',')

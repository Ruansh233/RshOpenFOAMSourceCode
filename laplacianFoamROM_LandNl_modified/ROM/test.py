from scipy.integrate import solve_ivp
import numpy as np
from newtonOdes import newtonOdesFunc

# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"

ddt1 = 1.0e-7
ddt2 = 1.0e-3
modesNum = 15

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

nonLinearCoeffTensor = np.empty((0, modesNum))

for i in range(0, modesNum):
    for j in range(0, modesNum):
        rowN = i * nonLinearCoeffLen + j
        # print(rowN)
        nonLinearCoeffTensor = np.append(nonLinearCoeffTensor, nonLinearCoeffMatrix[rowN, :].reshape(1, modesNum), 0)

# nonliear tensor
nonLinearCoeffTensor = nonLinearCoeffTensor.reshape(modesNum, modesNum, modesNum)

# print(nonLinearCoeffTensor)
# print("\n____________________________________\n")
# print(nonLinearCoeffTensor[0])

# print(diffuTermCoeffMatrix)
# print("\n____________________________________\n")
# print(diffuTermCoeffMatrix[0])

# print(diffuTermCoeffMatrix)
# print(nonLinearCoeffTensor)

initialA = coeffMatrix[0, 0: modesNum]
print(initialA)

testA = coeffMatrix[900, 0: modesNum]

test1 = 1.0e-7*testA.dot(nonLinearCoeffTensor[0]).dot(testA) + 1.0e-3*diffuTermCoeffMatrix[0].dot(testA)
test2 = 1.0e-7*testA.dot(nonLinearCoeffTensor[1]).dot(testA) + 1.0e-3*diffuTermCoeffMatrix[1].dot(testA)
test3 = 1.0e-7*testA.dot(nonLinearCoeffTensor[2]).dot(testA) + 1.0e-3*diffuTermCoeffMatrix[2].dot(testA)

print(test1, test2, test3)

# dt = 0.01
# totalt = 1000
# An = np.zeros((5, modesNum))
# An = newtonOdesFunc(nonLinearCoeffTensor, diffuTermCoeffMatrix, initialA, dt, totalt, modesNum)
# # print(An.shape)

# np.savetxt("An", An, fmt='%.6e', delimiter=',')
# print("\n____________________________________\n")

# dA = np.zeros((modesNum))
# dA[0] = 1
# delta = dA[0]
# iniA_da = initialA + dA

# test1 = - 0.1 * (initialA.dot(nonLinearCoeffTensor[0]).dot(initialA) + diffuTermCoeffMatrix[0].dot(initialA))

# test2 = delta - 0.1 * (iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da))
# # test2 = iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da)

# print(test1, test2)
# print((test2  - test1)/delta)

# dA = np.zeros((modesNum))
# dA[1] = 1
# delta = dA[1]
# iniA_da = initialA + dA

# test1 = - 0.1 * (initialA.dot(nonLinearCoeffTensor[0]).dot(initialA) + diffuTermCoeffMatrix[0].dot(initialA))

# test2 = - 0.1 * (iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da))
# # test2 = iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da)

# print(test1, test2)
# print((test2  - test1)/delta)

# dA = np.zeros((modesNum))
# dA[2] = 1
# delta = dA[2]
# iniA_da = initialA + dA

# test1 = - 0.1 * (initialA.dot(nonLinearCoeffTensor[2]).dot(initialA) + diffuTermCoeffMatrix[2].dot(initialA))

# test2 = - 0.1 * (iniA_da.dot(nonLinearCoeffTensor[2]).dot(iniA_da) + diffuTermCoeffMatrix[2].dot(iniA_da))
# # test2 = iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da)

# print(test1, test2)
# print((test2  - test1)/delta)

# dA = np.zeros((modesNum))
# dA[0] = 1
# delta = dA[0]
# iniA_da = initialA + dA

# test1 = 1.0e-7*initialA.dot(nonLinearCoeffTensor[0]).dot(initialA) + 1.0e-3*diffuTermCoeffMatrix[0].dot(initialA)
# # test2 = iniA_da.dot(nonLinearCoeffTensor[0]).dot(iniA_da) + diffuTermCoeffMatrix[0].dot(iniA_da)

# print(test1)

time = np.linspace(0, 100, 1001)

def odefun(t, a):
    da = ddt1 * a.dot(nonLinearCoeffTensor).dot(a) + ddt2 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 100], testA, method='Radau', dense_output=True)

snapshots_cal = modeMatrix @ sol.sol(time)

coeff_calculate_test = filePath + "/coeff_calculate_test"
np.savetxt(f"{coeff_calculate_test}", sol.sol(time).T, fmt='%.6e', delimiter=',')



				


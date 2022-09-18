from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"

diffucoefficient = filePath + "/diffuTermCoeffMatrix"
coeff = filePath + "/coeffMatrix"
spatialmode = filePath + "/modeMatrix"
laplacianMode = filePath + "/laplacianModesMatrix"
coeff_calculate = filePath + "/coeff_calculate"
snapshots_calculate = filePath + "/snapshots_calculate"

diffuTermCoeffMatrix = np.loadtxt(diffucoefficient)

modeMatrix = np.loadtxt(spatialmode)
modeLapMatrix = np.loadtxt(laplacianMode)

coeffMatrix = np.loadtxt(coeff)
initialA = coeffMatrix[0, :]

print(diffuTermCoeffMatrix)
print(initialA)

time = np.linspace(0, 50, 501)


def odefun(t, a):
    da = 2.5e-3 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 50], initialA, dense_output=True)

# print(sol.y)
# print(sol.sol(time))

snapshots_cal = projmodeMatrix @ sol.sol(time)

np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')
np.savetxt(f"{snapshots_calculate}", snapshots_cal, fmt='%.6e', delimiter=',')

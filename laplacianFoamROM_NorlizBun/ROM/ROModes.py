from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"
diffucoefficient = filePath + "/diffuTermCoeffMatrix"
coeff = filePath + "/coeffMatrix"
coeff_calculate = filePath + "/coeff_calculate"

diffuTermCoeffMatrix = np.loadtxt(diffucoefficient)
coeffMatrix = np.loadtxt(coeff)
initialA = coeffMatrix[0, :np.shape(diffuTermCoeffMatrix)[1]]

print(diffuTermCoeffMatrix)
print(initialA)

time = np.linspace(0, 200, 2001)


def odefun(t, a):
    da = 2.5e-3 * diffuTermCoeffMatrix.dot(a)
    return da

sol = solve_ivp(odefun, [0, 200], initialA, dense_output=True)

# print(sol.y)
# print(sol.sol(time))

np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')

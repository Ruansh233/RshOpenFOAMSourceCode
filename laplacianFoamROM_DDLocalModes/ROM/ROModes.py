from lib2to3.pgen2.token import RPAR
from scipy.integrate import solve_ivp
import numpy as np
# import matplotlib.pyplot as plt

filePath = "../svdtest/SVD"


for domainI in range(0, 4):

    ROMdiffucoefficient = filePath + "/subDROMCoeffMatrix" + str(domainI)
    coeff = filePath + "/temporalCoeff" + str(domainI)
    coeff_calculate = filePath + "/coeff_calculate" + str(domainI)

    diffuTermCoeffMatrix = np.loadtxt(ROMdiffucoefficient)

    coeffMatrix = np.loadtxt(coeff)
    initialA = coeffMatrix[0, :np.shape(diffuTermCoeffMatrix)[1]]

    # print(diffuTermCoeffMatrix)

    print(initialA)
    print(domainI)
    # print("")

    time = np.linspace(0, 20, 101)


    def odefun(t, a):
        da = 2.5e-3 * diffuTermCoeffMatrix.dot(a)
        return da

    sol = solve_ivp(odefun, [0, 20], initialA, method='BDF', dense_output=True)

    np.savetxt(f"{coeff_calculate}", sol.sol(time).T, fmt='%.6e', delimiter=',')

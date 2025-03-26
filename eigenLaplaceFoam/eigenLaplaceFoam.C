/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Group
    grpBasicSolvers

Description
    Laplace equation solver for a scalar quantity.

    \heading Solver details
    The solver is applicable to, e.g. for thermal diffusion in a solid.  The
    equation is given by:

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "EigenLaplace solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "readEigenFunctionDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "Computing Laplacian eigenfunctions\n" << endl;

    // List to store computed eigenfunctions
    PtrList<volScalarField> eigenModes(nEigen);
    List<scalar> eigenValues(nEigen);

    // Compute nEigen smallest eigenfunctions
    for (label m = 0; m < nEigen; ++m)
    {
        Info << "Computing eigenmode " << m + 1 << endl;
        runTime.setTime(scalar(m + 1), m+1);

        // Initialize φ with a random field (or a sine function)

        forAll(U.internalField(), cellI)
        {
            U[cellI] = Foam::sin(2.0 * Foam::constant::mathematical::pi * mesh.C()[cellI][0]) *
                        Foam::sin(2.0 * Foam::constant::mathematical::pi * mesh.C()[cellI][1]);
        }
        U.correctBoundaryConditions();

        scalar lambda_old = 0, lambda_new = 0;

        // Power iteration with deflation
        for (int iter = 0; iter < maxIter; iter++)
        {
            // Solve Laplace equation: ∇²v = U
            fvScalarMatrix laplacianEqn(fvm::laplacian(nu, v) + U);
            laplacianEqn.solve();

            // Deflate previously computed modes, The Gram–Schmidt process
            for (int j = 0; j < m; j++)
            {
                volScalarField &prevU = eigenModes[j];
                scalar coeff = fvc::domainIntegrate(v * prevU).value() / 
                                fvc::domainIntegrate(sqr(prevU)).value();
                v -= coeff * prevU;
            }

            // Normalize φ
            scalar vNorm = Foam::sqrt(fvc::domainIntegrate(sqr(v)).value());
            U = v / vNorm;
            U.correctBoundaryConditions();

            // Compute new eigenvalue using Rayleigh quotient
            lambda_new = fvc::domainIntegrate(U * fvc::laplacian(nu, U)).value() / 
                        fvc::domainIntegrate(sqr(U)).value();

            // Check for convergence
            if (mag(lambda_new - lambda_old) < tol)
            {
                Info<< "Converged! Final eigenvalue: " << lambda_new << nl;
                Info<< "The number of iterations: " << iter << nl;
                break;
            }
            lambda_old = lambda_new;
        }

        // // rename the computed mode
        // U.rename("U"+ Foam::name(m+1));
        // Write the computed mode
        U.write();
        // Store the computed mode
        eigenModes.set(m, U.clone());
        eigenValues[m] = lambda_new;
    }

    IOList<scalar> testListIO
    (
        IOobject
        (
            "eigenValues",
            runTime.rootPath(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        eigenValues
    );
    testListIO.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

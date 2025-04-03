/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    simpleFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, turbulent flows.

    \heading Solver details
    The solver uses the SIMPLE algorithm to solve the continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for incompressible, turbulent flows."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    #include "readEigenFunctionDict.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    dimensionedScalar nu("nu", dimViscosity, 1.0);  // Kinematic viscosity
    dimensionedScalar mu_("mu", dimensionSet(0, 0, -1, 0, 0, 0, 0), 1.0); 

    volScalarField p0 = p;
    volVectorField U0 = U;

    // List to store computed eigenfunctions
    PtrList<volVectorField> eigenFunctions(nEigen);
    List<scalar> eigenValues(nEigen);

    for (label eigenIndex = 0; eigenIndex < nEigen; eigenIndex++)
    {
        Info << "Computing eigenmode " << eigenIndex + 1 << endl;
        runTime.setTime(scalar(eigenIndex + 1), eigenIndex+1);

        scalar lambda_old = 0, lambda_new = 0;

        p = p0;
        U = U0;

        // Inner loop: Inverse power method
        for (label iter = 0; iter < maxIter; iter++)
        {
            // volVectorField source = U;

            // --- Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "pEqn.H"
            }

            laminarTransport.correct();
            turbulence->correct();

            // Deflate the source term
            for (label j = 0; j < eigenIndex; j++)
            {
                scalar coeff = fvc::domainIntegrate(U & eigenFunctions[j]).value() /
                               fvc::domainIntegrate(eigenFunctions[j] & eigenFunctions[j]).value();
                U -= coeff * eigenFunctions[j];
            }

            // Normalize U and update U0
            U /= Foam::sqrt(fvc::domainIntegrate(magSqr(U)).value());

            // Compute Rayleigh quotient for eigenvalue estimate
            scalar lambda_new = fvc::domainIntegrate(U & (fvc::laplacian(nu, U) + fvc::grad(p))).value() /
                            fvc::domainIntegrate(magSqr(U)).value();

            // Check for convergence
            if (mag(lambda_new - lambda_old)/mag(lambda_new) < tol)
            {
                Info<< "Converged! Final eigenvalue: " << lambda_new << nl;
                Info<< "The number of iterations: " << iter << nl;
                break;
            }
            lambda_old = lambda_new;
        }

        // Write the computed mode
        U.write();
        // Store the converged eigenfunction
        eigenFunctions.set(eigenIndex, U.clone());
        eigenValues[eigenIndex] = lambda_new;

        runTime.printExecutionTime(Info);
    }

    return 0;
}


// ************************************************************************* //

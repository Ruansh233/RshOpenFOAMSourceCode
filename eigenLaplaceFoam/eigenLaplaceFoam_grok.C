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

    // List to store eigenfunctions
    PtrList<volScalarField> ef_list;

    // Compute nEigen smallest eigenfunctions
    for (label m = 0; m < nEigen; ++m)
    {
        // Random initial guess (ensure Dirichlet BCs if applicable)
        forAll(u, cellI) 
        {
            u[cellI] = rand() / scalar(RAND_MAX); 
        }
        u.correctBoundaryConditions();

        scalar lambda = 0.0, lambdaPrev = 0.0;

        // Inverse power method iterations
        for (label iter = 0; iter < maxIter; ++iter)
        {
            // b = V * u
            b.ref() = mesh.V() * u;
            b.correctBoundaryConditions();

            // Solve laplacian(v) = b
            fvScalarMatrix laplacianEqn(fvm::laplacian(v) - b);
            laplacianEqn.solve();

            // Orthogonalize against previous eigenfunctions
            forAll(ef_list, i)
            {
                volScalarField& ef = ef_list[i];
                scalar coef = fvc::domainIntegrate(v * ef).value() / 
                                fvc::domainIntegrate(sqr(ef)).value();
                v -= coef * ef;
            }

            // Normalize
            scalar norm_v = Foam::sqrt(fvc::domainIntegrate(sqr(v)).value());
            volScalarField u_new = v / norm_v;

            // Estimate eigenvalue
            lambdaPrev = lambda;
            lambda = fvc::domainIntegrate(u_new * b).value() / 
                        fvc::domainIntegrate(sqr(u_new)).value();

            // Update u
            u = u_new;

            // Check convergence
            if (iter > 0 && mag(lambda - lambdaPrev) < tol)
            {
                Info<< "Eigenfunction " << m + 1 
                    << " converged. Eigenvalue: " << lambda << endl;
                break;
            }
        }

        // Write eigenfunction
        u.write();
        // Store eigenfunction
        ef_list.append(u.clone());
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

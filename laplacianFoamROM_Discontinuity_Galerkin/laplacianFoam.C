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
#include "SVD.H"
#include "fvOptions.H"
#include "simpleControl.H"
#include "cpuTimeCxx.H"
#include "writeMatrix.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Laplace equation solver for a scalar quantity."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    #include "readSVDDict.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;  
        
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T) - fvm::laplacian(DT, T)
             ==
                fvOptions(T)
            );

            fvOptions.constrain(TEqn);
            TEqn.solve();
            fvOptions.correct(T);
        }

        #include "write.H"

        runTime.printExecutionTime(Info);
    }

    Info<< "\nThe simulation of FOM is finishing at " << runTime.elapsedCpuTime() << " s.\n" << nl;

    // file location and OF pointer of the file
    fileName dataFile;
    autoPtr<OFstream> outputFilePtr;

    #include "svdSnapshots.H"

    Info<< "\nThe solver is finishing at " << runTime.elapsedCpuTime() << " s.\n" << nl;

    // // ------------------------------------------------------------------------
    // // Rsh, very important, access field value data from other cases

    // Info<< "\ntest:\n" << endl;

    // // create new time object for other cases
    // Foam::Time runTimeTest
    // (
    //     Foam::Time::controlDictName,
    //     args.rootPath(),
    //     // args.caseName(),
    //     "Block1",
    //     "system",
    //     "constant"
    // );

    // // create new mesh object for other cases
    // fvMesh  refElementMesh
    // (
    //     IOobject
    //     (
    //         polyMesh::defaultRegion,
    //         args.rootPath()/"Block1"/"constant",
    //         runTimeTest,
    //         IOobject::MUST_READ
    //     ),
    //     false
    // );

    // // create test modes
    // volScalarField testFieldValue
    // (
    //     IOobject
    //     (
    //         "T",
    //         runTimeTest.timeName(),
    //         refElementMesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     refElementMesh
    // );

    // // laplacian of test field value
    // volScalarField testlap
    // (
    //     IOobject
    //     (
    //         "testlap",
    //         runTimeTest.timeName(),
    //         refElementMesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::laplacian(testFieldValue)
    // );

    // // write the field to new cases folder
    // testlap.write();
    // // ------------------------------------------------------------------------

    return 0;
}


// ************************************************************************* //

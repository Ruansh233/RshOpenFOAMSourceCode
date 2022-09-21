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

        #include "setTransportProperties.H"

        #include "write.H"

        runTime.printExecutionTime(Info);

        if(ROMorNot)
        {
            // add snapshots in specific time interval
            #include "snapshotsMatrix.H"
        }
    }

    Info<< "End\n" << endl;

    Info<< "The simulation of FOM is finished in " << runTime.elapsedCpuTime() << " s." << nl;
 
    if(ROMorNot)
    {
        Info<< "snapshots number is: " << snapshotsNo << endl
            << "snapshotsTime is: " << snapshotsTime << endl;

        // check whether the number of snapshots is equal to predefined value,
        // -- sometimes the snapshotsNo = snapshotsNum-1 due to the double type of runTime at the last step
        if(snapshotsNo != snapshotsNum)
        {
            snapshotsM.resize(snapshotsRows, snapshotsNo); 
            snapshotsNum = snapshotsNo;  
        }
        svdDict.add("snapshotsNum", snapshotsNum, true);   
        svdDict.regIOobject::write();

        // svd of the snapshots matrix
        SVD fieldValueSVD(snapshotsM);   

        Info<< "The SVD process of snapshots are finished in " << runTime.elapsedCpuTime() << " s." << nl;

        // write data of modes
        #include "writeModesField.H"

        // file location and OF pointer of the file
        fileName dataFile;
        autoPtr<OFstream> outputFilePtr;

        // write diffusion terms coefficient
        #include "calculateDiffuCoeff.H"

        // wirte snapshots matrix, modes, eigenvalues and coefficient
        #include "writeprocessedMatrix.H"   
    }

    Info<< "The solver is finished in " << runTime.elapsedCpuTime() << " s.\n" << nl;

    return 0;
}


// ************************************************************************* //

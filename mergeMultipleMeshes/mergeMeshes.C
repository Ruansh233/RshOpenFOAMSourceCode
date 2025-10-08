/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    mergeMeshes

Group
    grpMeshManipulationUtilities

Description
    Merges two meshes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "mergePolyMesh.H"
#include "topoSet.H"
#include "processorMeshes.H"
#include "fvCFD.H"
#include <list>

using namespace Foam;

void getRootCase(fileName& casePath)
{
    casePath.clean();  // Remove unneeded ".."

    if (casePath.empty() || casePath == ".")
    {
        // handle degenerate form and '.'
        casePath = cwd();
    }
    else if (casePath[0] != '/' && casePath.name() == "..")
    {
        // avoid relative cases ending in '..' - makes for very ugly names
        casePath = cwd()/casePath;
        casePath.clean();  // Remove unneeded ".."
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge multiple meshes"
    );

    #include "addOverwriteOption.H"
    argList::addOption("dict", "file", "Alternative mergeMultipleMeshesDict");

    argList::addOption
    (
        "resultTime",
        "time",
        "Specify a time for the resulting mesh"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    #include "dictRead.H"

    const bool overwrite = args.found("overwrite");

    const word masterRegion =
        args.getOrDefault<word>("masterRegion", polyMesh::defaultRegion);

    getRootCase(masterCase);
    Info<< "Master:      " << masterCase  << nl
        << "mesh to add: " << addCaseList << endl;

    #include "createMasterTimes.H"

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info<< "Create mesh\n" << endl;
    mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );

    word meshInstance = masterMesh.pointsInstance();

    const bool specifiedInstance =
    (
        !overwrite
     && args.readIfPresent("resultTime", meshInstance)
    );

    if (specifiedInstance)
    {
        runTimeMaster.setTime(instant(meshInstance), 0);
    }
    else if (!overwrite)
    {
        runTimeMaster++;
    }

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    forAll(addCaseList, caseID)
    {
        fileName addCase(addCaseList[caseID]);

        const word addRegion = 
                args.getOrDefault<word>("addRegion", polyMesh::defaultRegion);

        getRootCase(addCase);

        #include "createAddTimes.H"

        Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;
        Info<< "Create mesh\n" << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                addRegion,
                runTimeToAdd.timeName(),
                runTimeToAdd
            )
        );

        masterMesh.addMesh(meshToAdd);
    }

    masterMesh.merge();

    if (overwrite || specifiedInstance)
    {
        masterMesh.setInstance(meshInstance);
    }

    masterMesh.write();
    topoSet::removeFiles(masterMesh);
    processorMeshes::removeFiles(masterMesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

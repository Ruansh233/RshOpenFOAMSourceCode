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

using namespace Foam;

void getRootCase(fileName& casePath)
{
    // Rsh
    // Info  << "casePath: " << casePath << endl;  

    casePath.clean();  // Remove unneeded ".."
    // Rsh
    // Info  << "casePath.clean(): " << casePath << endl; 

    // 
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

    // Rsh
    // Info  << "casePath.clean() two times: " << casePath << endl; 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge two meshes"
    );

    #include "addOverwriteOption.H"
    argList::addOption("dict", "file", "Alternative mergeMultipleMeshesDict");

    // argList::addArgument("masterCase");
    // argList::addOption
    // (
    //     "masterRegion",
    //     "name",
    //     "Specify alternative mesh region for the master mesh"
    // );

    // argList::addArgument("addCase");
    // argList::addOption
    // (
    //     "addRegion",
    //     "name",
    //     "Specify alternative mesh region for the additional mesh"
    // );

    // argList::addArgument("addCase2");
    // argList::addOption
    // (
    //     "addRegion",
    //     "name",
    //     "Specify alternative mesh region for the additional mesh"
    // );

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

    // argList args(argc, argv);
    // if (!args.check())
    // {
    //      FatalError.exit();
    // }

    const bool overwrite = args.found("overwrite");

    // auto masterCase = args.get<fileName>(1);
    // auto addCase = args.get<fileName>(2);

    // auto addCase2 = args.get<fileName>(3);

    // Rsh, masterCase: ".", addCase: "../testCase_2"
    // Info << "masterCase: " << masterCase << endl
    //      << "addCase: " << addCase << endl; 

    const word masterRegion =
        args.getOrDefault<word>("masterRegion", polyMesh::defaultRegion);

    // const word addRegion =
    //     args.getOrDefault<word>("addRegion", polyMesh::defaultRegion);

    // const word addRegion2 =
    //     args.getOrDefault<word>("addRegion2", polyMesh::defaultRegion);

    // Since we don't use argList processor directory detection, add it to
    // the casename ourselves so it triggers the logic inside TimePath.
    // Rsh, caseName Return case name (parallel run) or global case (serial run)
    // const fileName& cName = args.caseName();
    
    // Rsh, output, caseName: "testCase"
    // Info << "caseName: " << cName << endl;

    // const auto pos = cName.find("processor");

    // // Rsh
    // // Info << "pos: " << pos << endl;
    
    // if (pos != string::npos && pos != 0)
    // {
    //     fileName processorName = cName.substr(pos);
    //     masterCase += '/' + processorName;
    //     // addCase += '/' + processorName;
    //     // addCase2 += '/' + processorName;
    // }

    // Rsh, find is a function of c++ string 
    // const word testword_("hello Ruan Shenhui");
    // const auto testpos1 = testword_.find("hello");
    // const auto testpos2 = testword_.find("Ruan");

    // Info << "testpos1: " << testpos1 << endl
    //      << "testpos2: " << testpos2 << endl;

    getRootCase(masterCase);
    // getRootCase(addCase);
    // getRootCase(addCase2);

    // Info<< "Master:      " << masterCase << "  region " << masterRegion << nl
    //     << "mesh to add: " << addCase    << "  region " << addRegion << endl;
        // << "mesh2 to add: " << addCase2    << "  region " << addRegion2 << endl;
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

    // Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;
    // Info<< "Create mesh\n" << endl;
    // polyMesh meshToAdd
    // (
    //     IOobject
    //     (
    //         addRegion,
    //         runTimeToAdd.timeName(),
    //         runTimeToAdd
    //     )
    // );

    // Info<< "Reading mesh to add for time = " << runTimeToAdd2.timeName() << nl;
    // Info<< "Create mesh\n" << endl;
    // polyMesh meshToAdd2
    // (
    //     IOobject
    //     (
    //         addRegion2,
    //         runTimeToAdd2.timeName(),
    //         runTimeToAdd2
    //     )
    // );

    // Rsh, pointsInstance() Return the current instance directory for points(mesh).
    word meshInstance = masterMesh.pointsInstance();

    // Rsh, meshInstance: constant
    // Info << "meshInstance: " << meshInstance << endl;

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

    // masterMesh.addMesh(meshToAdd);
    // masterMesh.addMesh(meshToAdd2);
    // masterMesh.merge();

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

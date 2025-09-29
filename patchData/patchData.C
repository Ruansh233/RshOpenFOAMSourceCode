/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "IFstream.H"
#include "IOobjectList.H"

template <class Type>
void printList(
    const List<List<Type>> &local_tmp,
    fvMesh &mesh,
    word listName,
    fileName dirName)
{
    if (Pstream::master())
    {
        // make mesh.time().globalPath() + "/" + dirName if it does not exist
        fileName dir(mesh.time().globalPath() + "/" + dirName);
        if (!isDir(dir))
            mkDir(dir);

        List<Type> tmpList;
        forAll(local_tmp, procI)
        {
            tmpList.append(local_tmp[procI]);
        }

        IOList<Type> tmpListIO(
            IOobject(
                listName,
                dir,
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE),
            tmpList);

        OFstream os(dir + "/" + listName);
        os << tmpList << endl;
    }
}

template <class Type>
void interfaceDataIO(
    IOobjectList &objects,
    fvMesh &mesh,
    label interfaceID,
    wordRes &selectedFields,
    word interfaceName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    for (
        const IOobject &io :
        (
            selectedFields.empty()
                ? objects.csorted<fieldType>()
                : objects.csorted<fieldType>(selectedFields)))
    {
        fieldType data(io, mesh, false);
        List<List<Type>> tmp_data(Pstream::nProcs());
        tmp_data[Pstream::myProcNo()] = data.boundaryField()[interfaceID];
        Pstream::gatherList(tmp_data);
        printList(tmp_data, mesh, io.name() + "_" + mesh.time().timeName(), interfaceName);
    }
}

template <>
void interfaceDataIO<tensor>(
    IOobjectList &objects,
    fvMesh &mesh,
    label interfaceID,
    wordRes &selectedFields,
    word interfaceName)
{
    typedef GeometricField<tensor, fvPatchField, volMesh> fieldType;

    for (
        const IOobject &io :
        (
            selectedFields.empty()
                ? objects.csorted<fieldType>()
                : objects.csorted<fieldType>(selectedFields)))
    {
        fieldType data(io, mesh, false);
        List<List<vector>> tmp_data(Pstream::nProcs());
        tmp_data[Pstream::myProcNo()] = data.boundaryField()[interfaceID] & mesh.boundary()[interfaceID].nf();
        Pstream::gatherList(tmp_data);
        printList(tmp_data, mesh, io.name() + "_" + mesh.time().timeName(), interfaceName);
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::addOption(
        "interfaceName",
        "word",
        "name of the interface, e.g., interface_Dir.");

    argList::addOption(
        "fields",
        "wordRes",
        "fields to be processed, default is all fields.");

    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    label interfaceID(-1);
    word interfaceName;
    if (args.found("interfaceName"))
    {
        interfaceName = args.get<word>("interfaceName");

        interfaceID = mesh.boundary().findPatchID(interfaceName);

        if (interfaceID == -1)
        {
            FatalErrorIn("main")
                << "Cannot find the patch " << interfaceName << exit(FatalError);
        }

        Info << "interfaceName: " << mesh.boundary()[interfaceID].name() << endl;
    }
    else
    {
        FatalErrorIn("main")
            << "Please specify the interfaceName." << exit(FatalError);
    }

    List<List<vector>> local_tmp(Pstream::nProcs());
    local_tmp[Pstream::myProcNo()] = mesh.boundary()[interfaceID].Cf();
    Pstream::gatherList(local_tmp);
    printList(local_tmp, mesh, "Cf", interfaceName);

    instantList timeDirs = timeSelector::select0(runTime, args);

    wordRes selectedFields(args.get<wordRes>("fields"));

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        IOobjectList objects(mesh, runTime.timeName());

        interfaceDataIO<scalar>(objects, mesh, interfaceID, selectedFields, interfaceName);
        interfaceDataIO<vector>(objects, mesh, interfaceID, selectedFields, interfaceName);
        interfaceDataIO<tensor>(objects, mesh, interfaceID, selectedFields, interfaceName);
    }

    Info << "End\n"
         << endl;
    return 0;
}

// ************************************************************************* //

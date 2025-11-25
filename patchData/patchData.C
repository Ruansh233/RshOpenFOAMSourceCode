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
        tmp_data[Pstream::myProcNo()] = mesh.boundary()[interfaceID].nf() & data.boundaryField()[interfaceID];
        Pstream::gatherList(tmp_data);
        printList(tmp_data, mesh, io.name() + "_" + mesh.time().timeName(), interfaceName);
    }
}

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    argList::addOption(
        "interfaceName",
        "wordList",
        "name of the interface, e.g., '(interface_Dir)' or '(interface_Dir interface_Neu)'.");

    argList::addOption(
        "fields",
        "wordRes",
        "fields to be processed, default is all fields.");

    argList::addBoolOption(
        "noFields",
        "if found, no field data is processed.");

    argList::addBoolOption(
        "patchInfo",
        "if found, patch information is processed.");

    argList::addBoolOption(
        "centre",
        "if found, output the field of patch centres.");

    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    labelList interfaceID;

    if (args.found("noFields") && args.found("fields"))
    {
        FatalErrorIn("main")
            << "Cannot specify both noFields and fields options." << exit(FatalError);
    }

    if (!args.found("interfaceName"))
    {
        FatalErrorIn("main")
            << "Please specify the interfaceName." << exit(FatalError);
    }
    wordList interfaceName(args.get<wordList>("interfaceName"));
    forAll(interfaceName, i)
    {
        interfaceID.append(mesh.boundary().findPatchID(interfaceName[i]));
    }
    forAll(interfaceID, i)
    {
        if (interfaceID[i] == -1)
        {
            FatalErrorIn("main")
                << "Cannot find patch " << interfaceName[i] << exit(FatalError);
        }

        Info<< "Found patch " << interfaceName[i] << endl;
    }

    if (args.found("centre"))
    {
        forAll(interfaceID, i)
        {
            List<List<vector>> local_tmp(Pstream::nProcs());
            local_tmp[Pstream::myProcNo()] = mesh.boundary()[interfaceID[i]].Cf();
            Pstream::gatherList(local_tmp);
            printList(local_tmp, mesh, "Cf", interfaceName[i]);
        }
    }

    if (args.found("patchInfo"))
    {
        forAll(interfaceID, i)
        {   
            List<scalar> local_tmp(Pstream::nProcs());
            local_tmp[Pstream::myProcNo()] = mesh.boundary()[interfaceID[i]].size();
            Pstream::gatherList(local_tmp);
            Info<< "Patch info:" << endl
                << "\t ------------------------" << endl
                << "\t Patch name: " << interfaceName[i] << endl
                << "\t Patch type: " << mesh.boundary()[interfaceID[i]].type() << endl
                << "\t Total number of faces: " << sum(local_tmp) << endl;

            forAll(local_tmp, proc)
            {
                if (local_tmp[proc] >= 1)
                {
                    Info<< "\t Processor " << proc
                        << " has " << local_tmp[proc] << " faces." << endl;
                }
            }
            Info<< "\t ------------------------" << endl;
        }
    }

    if (args.found("fields"))
    {
        instantList timeDirs = timeSelector::select0(runTime, args);

        wordRes selectedFields(args.get<wordRes>("fields"));

        forAll(timeDirs, timei)
        {
            runTime.setTime(timeDirs[timei], timei);

            IOobjectList objects(mesh, runTime.timeName());

            forAll(interfaceName, i)
            {
                interfaceDataIO<scalar>(objects, mesh, interfaceID[i], selectedFields, interfaceName[i]);
                interfaceDataIO<vector>(objects, mesh, interfaceID[i], selectedFields, interfaceName[i]);
                interfaceDataIO<tensor>(objects, mesh, interfaceID[i], selectedFields, interfaceName[i]);
            }
        }
    }

    Info<< "End\n"
        << endl;
    return 0;
}

// ************************************************************************* //

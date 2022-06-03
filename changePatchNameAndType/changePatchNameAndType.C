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
#include <iomanip>

word replaceStringPart(word& str, const word& sub, const word& mod) 
{
    word tmp(str);
    while (tmp.find(sub) != -1)
    {
        tmp.insert(tmp.find(sub), mod);
        tmp.erase(tmp.find(sub), sub.length());
    }
    return tmp;
}

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
   
    // Rsh, write file head
    IOdictionary createPatchDict
    (
        IOobject
        (
            "createPatchDict_interfaces",
            runTime.system(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    createPatchDict.regIOobject::write();
    
    fileName outputFile("createPatchDict_interfaces");
    OFstream writeFile(runTime.path()/runTime.system()/outputFile, IOstreamOption(), true);

    writeFile << "\npointSync false;\npatches\n(";

    // Rsh, interface pair
    const word emptyDict = R"deli(
    {
        name                interface_pair1;
        patchInfo
        {
            type            cyclic;
            neighbourPatch  interface_pair2;
            transform       coincidentFullMatch;
        }
        constructFrom patches;
        patches (pair1);
    }

    {
        name                interface_pair2;
        patchInfo
        {
            type            cyclic;
            neighbourPatch  interface_pair1;
            transform       coincidentFullMatch;
        }
        constructFrom patches;
        patches (pair2);
    })deli";

    forAll(interfacePatchListA, patchListI)
    {
        word interfacePairDict(emptyDict);
        interfacePairDict = replaceStringPart(interfacePairDict, "pair1", interfacePatchListA[patchListI]);
        interfacePairDict = replaceStringPart(interfacePairDict, "pair2", interfacePatchListB[patchListI]);

        writeFile << interfacePairDict << nl;
    }

    const word endOfFile = R"deli(
);

// ************************************************************************* //
    )deli";

    writeFile << "\n);"
              << "\n// ************************************************************************* //"; 
    // writeFile.close();

//     // system("cat system/createPatchDict");
//     // system("createPatch");

}

// ************************************************************************* //

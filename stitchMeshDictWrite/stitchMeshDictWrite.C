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

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
    
    // Rsh, read interfaces patch pair
    List<word> interfacePatchListA(mesh.boundaryMesh().size());
    List<word> interfacePatchListB(mesh.boundaryMesh().size());
    scalar junk(0.0);
    label interfacePairNumber(0);

    fileName interfacedataFile (runTime.path()/runTime.system()/"patchPair");
    IFstream dataStream(interfacedataFile);
    while (! dataStream.eof())
    {
        dataStream >> interfacePatchListA[interfacePairNumber];
        dataStream >> interfacePatchListB[interfacePairNumber];
        dataStream.read(junk);
        dataStream.read(junk);
        dataStream.read(junk);

        interfacePairNumber += 1;
    }

    interfacePatchListA.resize(interfacePairNumber);
    interfacePatchListB.resize(interfacePairNumber);

    // Rsh, write file head
    IOdictionary stitchMeshDict
    (
        IOobject
        (
            "stitchMeshDict",
            runTime.system(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    forAll(interfacePatchListA, I)
    {
        word stitchPairName("interface" + name(I));
        stitchMeshDict.subDictOrAdd(stitchPairName);

        stitchMeshDict.subDict(stitchPairName).add("match", "perfect");
        stitchMeshDict.subDict(stitchPairName).add("master", interfacePatchListA[I]);
        stitchMeshDict.subDict(stitchPairName).add("slave", interfacePatchListB[I]);

    }

    stitchMeshDict.regIOobject::write();


    // stitchPair
    // {
    //     match   perfect;
    //     master  patchMaster;
    //     slave   patchSlave;
    // }
    // })deli";


}

// ************************************************************************* //

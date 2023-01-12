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
    argList::addOption // string variable
    (
        "interfaceFile",
        "word",
        "the interface file name(directory)"
    );

    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
    
    // Rsh, read interfaces patch pair
    List<word> interfacePatchListA(mesh.boundaryMesh().size());
    List<word> interfacePatchListB(mesh.boundaryMesh().size());
    label interfacePairNumber(0);

    word interfaceFile;
    if(args.readIfPresent("interfaceFile", interfaceFile, word ("patchPair")))
        Info<< "no interfaceFile inputed, default file is patchPair.\n";
    fileName interfacedataFile (runTime.path()/runTime.system()/interfaceFile);

    IFstream dataStream(interfacedataFile);
    word dataLine;
    token readToken;

    while(dataStream.getLine(dataLine) && !isNull(dataLine))
    {
        IStringStream dataString (dataLine);
        dataString.read(readToken);
        interfacePatchListA[interfacePairNumber] = readToken.wordToken();
        dataString.read(readToken);
        interfacePatchListB[interfacePairNumber] = readToken.wordToken();

        // dataString >> interfacePatchListA[interfacePairNumber];
        // dataString >> interfacePatchListB[interfacePairNumber];
        interfacePairNumber += 1;
    }     

    // while (! dataStream.eof())
    // {
    //     dataStream >> interfacePatchListA[interfacePairNumber];
    //     dataStream >> interfacePatchListB[interfacePairNumber];
    //     dataStream.read(junk);
    //     dataStream.read(junk);
    //     dataStream.read(junk);

    //     interfacePairNumber += 1;
    // }

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
}

// ************************************************************************* //

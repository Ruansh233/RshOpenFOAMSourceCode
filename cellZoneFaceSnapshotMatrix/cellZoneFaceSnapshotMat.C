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

// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

#include "fvCFD.H"
#include "wordRe.H"

int main(int argc, char *argv[])
{
    // Rsh, add dict select function
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative cellZoneDict"
    );

    argList::addOption
    (
        "cellZonePatch",
        "word",
        "Specify the fileName of cellZonePatch"
    );

    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // add dict Selector option
    const word dictName("cellZoneDict");
    #include "setSystemMeshDictionaryIO.H"
    Info<< "Reading " << dictIO.instance()/dictIO.name() << nl << endl;

    // create dictionary object and read the dictionary
    IOdictionary customDict(dictIO);
    #include "readDict.H"


    List<word> allBoundaryNames(customDict.lookup("allBoundaryNames"));


    word outputDir(customDict.getWord("outputDir"));
    // create the SVD folder if it is not exist.
    if(!isDir(runTime.globalPath()/outputDir))
        mkDir(runTime.globalPath()/outputDir);


    fileName dataFile;


    // read the fields
    List<word> scalarFields(customDict.lookup("scalarFields"));
    List<word> vectorFields(customDict.lookup("vectorFields"));


    fileName snapshotFileName;
    autoPtr<OFstream> outputFilePtr;


    word cellZonePatch(args.getOrDefault("cellZonePatch", word("cellZonePatchID")));
    // IOList for list input and output
    IOList<List<label>> zonePatchID
    (
        IOobject
        (
            cellZonePatch,
            runTime.caseConstant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    
    // // patchNameList, used for avoid parallel matching of mesh
    // List<word> patchNameList;
    // forAll(mesh.boundary(), patchI)
    // {
    //     patchNameList.append(mesh.boundary()[patchI].name());
    // }

    // // Info<< "patchNameList: " << patchNameList << endl;

    Pout<< "patch size: " << mesh.boundary().size() <<endl;


    // output scalarFields boundary value
    forAll(scalarFields, fieldI)
    {
        forAll(typesName, typeI)
        {
            List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
            forAll(boundaryNames, bounI)
            {
                wordRe patchMatch(".*"+boundaryNames[bounI]);
                patchMatch.compile();
                Info<< "write boundary field value: " << scalarFields[fieldI]+name(typeI)+"_"+boundaryNames[bounI] << endl;

                snapshotFileName = scalarFields[fieldI]+name(typeI)+"_"+boundaryNames[bounI];
                if (Pstream::master())
                    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

                for(label caseI=0; caseI<casesNumber; ++caseI)
                {                   
                    volScalarField tmpField
                    (
                        IOobject
                        (
                            scalarFields[fieldI]+"_"+name(CasesID[caseI]),
                            runTime.timeName(),
                            mesh,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh
                    );

                    forAll(cellZonesType[typeI][caseI], zoneI)
                    {
                        List<label> patchIDList (zonePatchID[cellZonesType[typeI][caseI][zoneI]]);
                        // Info<< mesh.boundary()[patchIDList.last()].name() << endl;
                        
                        List<word> patchName;
                        List<label> patchIDs;
                        List<scalar> faceValue;
                        List<scalar> faceValue_;

                        forAll(patchIDList, patchI)
                        {
                            if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name())
                                && ! mesh.boundary()[patchIDList[patchI]].name().starts_with("procBoundary"))
                            {                               
                                Pout<< "1-the matched patch name: " << mesh.boundary()[patchIDList[patchI]].name() << endl
                                    // << "the patchNameList name: " << patchNameList[patchIDList[patchI]] << endl
                                    << "1-the matched patch size: " << mesh.boundary()[patchIDList[patchI]].size() << endl;
                                
                                patchName.append(mesh.boundary()[patchIDList[patchI]].name());

                                // collected face value from individual processors
                                forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                                {
                                    faceValue.append(tmpField.boundaryField()[patchIDList[patchI]][faceI]);
                                }
                            }

                            List< List<scalar> > gatheredFaceValue(Pstream::nProcs());
                            gatheredFaceValue[Pstream::myProcNo()] = faceValue;
                            Pstream::gatherList(gatheredFaceValue);
                            Pstream::scatterList(gatheredFaceValue);
                            faceValue_  = ListListOps::combine<List<scalar> >(gatheredFaceValue, accessOp<List<scalar> >());

                            Info << "1-faceValue_.size: " << faceValue_.size() << endl;

                            if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name())
                                && mesh.boundary()[patchIDList[patchI]].name().starts_with("procBoundary"))
                            {
                                Pout<< "2-the matched patch name: " << mesh.boundary()[patchIDList[patchI]].name() << endl
                                    // << "the patchNameList name: " << patchNameList[patchIDList[patchI]] << endl
                                    << "2-the matched patch size: " << mesh.boundary()[patchIDList[patchI]].size() << endl;
                                
                                // collected face value from individual processors
                                forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                                {
                                    faceValue_.append(tmpField.boundaryField()[patchIDList[patchI]][faceI]);
                                }
                            }

                            Pout<< "2-faceValue_.size: " << faceValue_.size() << endl;

                        }

                        // forAll(patchIDList, patchI)
                        // {
                        //     if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
                        //     {                               
                        //         patchName.append(mesh.boundary()[patchIDList[patchI]].name());
                        //         patchIDs.append(patchIDList[patchI]);

                        //         Pout<< "the matched patch name: " << mesh.boundary()[patchIDList[patchI]].name() << endl
                        //             // << "the patchNameList name: " << patchNameList[patchIDList[patchI]] << endl
                        //             << "the matched patch size: " << mesh.boundary()[patchIDList[patchI]].size() << endl;
                        //     }
                        // }

                        // List<scalar> faceValue;
                        // forAll(patchIDs, patchI)
                        // {
                        //     // collected face value from individual processors
                        //     forAll(mesh.boundary()[patchIDs[patchI]], faceI)
                        //     {
                        //         faceValue.append(tmpField.boundaryField()[patchIDs[patchI]][faceI]);
                        //     }
                        // }

                        // // Info only write the value on the master processor
                        // Info<< "patchName: " << patchName << endl
                        //     << "faceValue, " << faceValue.size() << endl;

                        

                        for(label i=0; i < Pstream::nProcs(); ++i)
                        {
                            // Info<< "gatheredFaceValue_" << i << ": " << gatheredFaceValue[i].size() << endl;
                            Pout<< "3-faceValue_" << i << ": " << faceValue_.size() << endl;
                        }

                        if (Pstream::master())
                        {
                            // Info<< "patchName: " << patchName << endl
                            //     << "faceValue_, " << faceValue_.size() << endl;

                            forAll(faceValue_, faceI)
                            {
                                outputFilePtr() << faceValue_[faceI] << " ";
                            }
                            outputFilePtr() << endl;
                        }
                    }
                }
            }
        }
    }

    // // output vector field
    // forAll(vectorFields, fieldI)
    // {
    //     forAll(typesName, typeI)
    //     {
    //         List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
    //         forAll(boundaryNames, bounI)
    //         {
    //             wordRe patchMatch(".*"+boundaryNames[bounI]);
    //             patchMatch.compile();
    //             Info<< "write boundary field value: " << vectorFields[fieldI]+name(typeI)+"_"+boundaryNames[bounI] << endl;

    //             forAll(cellZonesType[typeI], caseI)
    //             {
    //                 volVectorField tmpField
    //                 (
    //                     IOobject
    //                     (
    //                         vectorFields[fieldI]+"_"+name(CasesID[caseI]),
    //                         runTime.timeName(),
    //                         mesh,
    //                         IOobject::MUST_READ,
    //                         IOobject::NO_WRITE
    //                     ),
    //                     mesh
    //                 );
                    
    //                 snapshotFileName = vectorFields[fieldI]+name(typeI)+"bun"+name(bounI);
    //                 if (Pstream::master())
    //                     outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));


    //                 forAll(cellZonesType[typeI][caseI], zoneI)
    //                 {
    //                     List<label> patchIDList (zonePatchID[cellZonesType[typeI][caseI][zoneI]]);
                        
    //                     forAll(patchIDList, patchI)
    //                     {
    //                         if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
    //                         {
    //                             // collected face value from individual processors
    //                             List<vector> faceValue;
    //                             forAll(mesh.boundary()[patchIDList[patchI]], faceI)
    //                             {
    //                                 faceValue.append(tmpField.boundaryField()[patchIDList[patchI]][faceI]);
    //                             }
    //                             List< List<vector> > gatheredFaceValue(Pstream::nProcs());
    //                             gatheredFaceValue[Pstream::myProcNo()] = faceValue;
    //                             Pstream::gatherList(gatheredFaceValue);
    //                             List<vector> faceValue_  = ListListOps::combine<List<vector> >(gatheredFaceValue, accessOp<List<vector> >());


    //                             if (Pstream::master())
    //                             {
    //                                 forAll(faceValue_, faceI)
    //                                 {
    //                                     outputFilePtr() << faceValue_[faceI].x() << " ";
    //                                 }
    //                                 forAll(faceValue_, faceI)
    //                                 {
    //                                     outputFilePtr() << faceValue_[faceI].y() << " ";
    //                                 }
    //                                 forAll(faceValue_, faceI)
    //                                 {
    //                                     outputFilePtr() << faceValue_[faceI].z() << " ";
    //                                 }
    //                             }
    //                         }
    //                     }
                        
    //                     if (Pstream::master())
    //                         outputFilePtr() << endl;
    //                 }
    //             }
    //         }
    //     }
    // }
    

    Info << "\nEnd\n" << endl;

    return 0;

}


// ************************************************************************* //

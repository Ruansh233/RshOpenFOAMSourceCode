// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone


word cellZonePatch(args.getOrDefault("cellZonePatch", word("cellZonePatch")));
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


// output pressure boundary value
forAll(typesName, typeI)
{
    List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
    forAll(boundaryNames, bounI)
    {
        wordRe patchMatch(".*"+boundaryNames[bounI]);
        patchMatch.compile();

        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, p" << name(typeI)+"_"+boundaryNames[bounI] << ", " 
            << zoneCells << ", " << cellZonesType[typeI].size() << endl;

        snapshotFileName = "p"+name(typeI)+"_"+boundaryNames[bounI];
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {
            List<label> patchIDList (zonePatchID[cellZonesType[typeI][zoneI]]);

            forAll(patchIDList, patchI)
            {
                if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
                {
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << p.boundaryField()[patchIDList[patchI]][faceI] << " ";
                    }
                }
            }
            outputFilePtr() << endl;
        }
    }
}

// output U velocity field
forAll(typesName, typeI)
{
    List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
    forAll(boundaryNames, bounI)
    {
        wordRe patchMatch(".*"+boundaryNames[bounI]);
        patchMatch.compile();

        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, U" << name(typeI)+"_"+boundaryNames[bounI] << ", " 
            << 3*zoneCells << ", " << cellZonesType[typeI].size() << endl;
            
        snapshotFileName = "U"+name(typeI)+"_"+boundaryNames[bounI];
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {
            List<label> patchIDList (zonePatchID[cellZonesType[typeI][zoneI]]);
            
            forAll(patchIDList, patchI)
            {
                if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
                {
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << U.boundaryField()[patchIDList[patchI]][faceI].x() << " ";
                    }
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << U.boundaryField()[patchIDList[patchI]][faceI].y() << " ";
                    }
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << U.boundaryField()[patchIDList[patchI]][faceI].z() << " ";
                    }
                }
            }
            outputFilePtr() << endl;
        }
    }
}


// output kinetic energy
forAll(typesName, typeI)
{
    List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
    forAll(boundaryNames, bounI)
    {
        wordRe patchMatch(".*"+boundaryNames[bounI]);
        patchMatch.compile();       
        
        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, k" << name(typeI)+"_"+boundaryNames[bounI] << ", " 
            << zoneCells << ", " << cellZonesType[typeI].size() << endl;

        snapshotFileName = "k"+name(typeI)+"_"+boundaryNames[bounI];
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {
            List<label> patchIDList (zonePatchID[cellZonesType[typeI][zoneI]]);

            forAll(patchIDList, patchI)
            {
                if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
                {
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << turbKinEnergy.boundaryField()[patchIDList[patchI]][faceI] << " ";
                    }
                }
            }
            outputFilePtr() << endl;
        }
    }
}


// output epsilon
forAll(typesName, typeI)
{
    List<word> boundaryNames(customDict.subDict(typesName[typeI]).lookup("boundaryNames"));
    forAll(boundaryNames, bounI)
    {
        wordRe patchMatch(".*"+boundaryNames[bounI]);
        patchMatch.compile();

        label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
        Info<< "Size of the snapshots, epsilon" << name(typeI)+"_"+boundaryNames[bounI] << ", " 
            << zoneCells << ", " << cellZonesType[typeI].size() << endl;

        snapshotFileName = "epsilon"+name(typeI)+"_"+boundaryNames[bounI];
        outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

        forAll(cellZonesType[typeI], zoneI)
        {
            List<label> patchIDList (zonePatchID[cellZonesType[typeI][zoneI]]);

            forAll(patchIDList, patchI)
            {
                if(patchMatch.match(mesh.boundary()[patchIDList[patchI]].name()))
                {
                    forAll(mesh.boundary()[patchIDList[patchI]], faceI)
                    {
                        outputFilePtr() << turbEpsilon.boundaryField()[patchIDList[patchI]][faceI] << " ";
                    }
                }
            }
            outputFilePtr() << endl;
        }
    }
}



// ************************************************************************* //

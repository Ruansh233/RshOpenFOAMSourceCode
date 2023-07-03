// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone


// p, pressure cell value
const label zoneCells(globalcellIndexList[0].size());
Info<< "Size of the snapshots, p: " << zoneCells << endl;

fileName snapshotFileName = "pSnapshots";
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(globalcellIndexList, domainI)
{                   
    forAll(globalcellIndexList[domainI], cellI)
    {
        label cellID (globalcellIndexList[domainI][cellI]);
        outputFilePtr() << p[cellID] << " ";
    }
    outputFilePtr() << endl;  
}


// vectorFields cell value
snapshotFileName = "USnapshots";
Info<< "Size of the snapshots, u: " << 3*zoneCells << endl;
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(globalcellIndexList, domainI)
{                   
    forAll(globalcellIndexList[domainI], cellI)
    {
        label cellID (globalcellIndexList[domainI][cellI]);
        outputFilePtr() << U[cellID].x() << " ";
    }

    forAll(globalcellIndexList[domainI], cellI)
    {
        label cellID (globalcellIndexList[domainI][cellI]);
        outputFilePtr() << U[cellID].y() << " ";
    }

    forAll(globalcellIndexList[domainI], cellI)
    {
        label cellID (globalcellIndexList[domainI][cellI]);
        outputFilePtr() << U[cellID].z() << " ";
    }

    outputFilePtr() << endl;  
} 


// ************************************************************************* //

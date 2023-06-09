// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone


// p, pressure cell value
label subdomainCells (0);
forAll(subdomains[0], zoneI)
{
    subdomainCells += mesh.cellZones()[subdomains[0][zoneI]].size();
}

Info<< "Size of the Snapshots, p: " 
    << subdomainCells << ", " << subdomains.size() << endl;

snapshotFileName = "pSnapshots";
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(subdomains, domainI)
{  
    forAll(subdomains[domainI], zoneI)
    {                   
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << p[cellN] << " ";
        }
    }

    outputFilePtr() << endl;  
}


// vectorFields cell value
Info<< "Size of the Snapshots, U: " 
    << subdomainCells << ", " << subdomains.size() << endl;

snapshotFileName = "USnapshots";
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(subdomains, domainI)
{
    forAll(subdomains[domainI], zoneI)
    {
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].x() << " "; 
        }        
    } 

    forAll(subdomains[domainI], zoneI)
    {
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].y() << " "; 
        }
    } 

    forAll(subdomains[domainI], zoneI)
    {
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].z() << " "; 
        }
    } 

    outputFilePtr() << endl;   
}


// output kinetic energy
Info<< "Size of the Snapshots, k: " 
    << subdomainCells << ", " << subdomains.size() << endl;

snapshotFileName = "kSnapshots";
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(subdomains, domainI)
{
    forAll(subdomains[domainI], zoneI)
    {
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << turbKinEnergy[cellN] << " ";
        }
    } 
    outputFilePtr() << endl;  
}


// output epsilon
Info<< "Size of the Snapshots, epsilon: " 
    << subdomainCells << ", " << subdomains.size() << endl;

snapshotFileName = "epsilonSnapshots";
outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

forAll(subdomains, domainI)
{
    forAll(subdomains[domainI], zoneI)
    {
        forAll(mesh.cellZones()[subdomains[domainI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[subdomains[domainI][zoneI]][cellI]);
            outputFilePtr() << turbEpsilon[cellN] << " ";
        }
    }
    outputFilePtr() << endl;  
}

// ************************************************************************* //

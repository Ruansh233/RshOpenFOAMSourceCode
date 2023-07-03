// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone

word outputDir(customDict.getWord("outputDir"));
// create the SVD folder if it is not exist.
if(!isDir(runTime.globalPath()/outputDir))
    mkDir(runTime.globalPath()/outputDir);


// // read the fields
// List<word> scalarFields(customDict.lookup("scalarFields"));
// List<word> vectorFields(customDict.lookup("vectorFields"));
// List<word> turbulenceFields(customDict.lookup("turbulenceFields"));


fileName snapshotFileName;
autoPtr<OFstream> outputFilePtr;

// p, pressure cell value
forAll(typesName, typeI)
{
    label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
    Info<< "Size of the snapshots, p" << name(typeI) << ", " 
        << zoneCells << ", " << cellZonesType[typeI].size() << endl;

    snapshotFileName = "p"+name(typeI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(cellZonesType[typeI], zoneI)
    {                   
        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << p[cellN] << " ";
        }
        outputFilePtr() << endl;  
    } 
}


// vectorFields cell value
forAll(typesName, typeI)
{           
    label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
    Info<< "Size of the snapshots, U" << name(typeI) << ", "
        << zoneCells << ", " << cellZonesType[typeI].size() << endl;

    snapshotFileName = "U"+name(typeI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(cellZonesType[typeI], zoneI)
    {
        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].x() << " "; 
        }

        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].y() << " "; 
        }

        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << U[cellN].z() << " "; 
        }
        
        outputFilePtr() << endl;   
    } 
}


// output kinetic energy
volScalarField turbKinEnergy = turbulence().k();
forAll(typesName, typeI)
{
    label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
    Info<< "Size of the snapshots, k" << name(typeI) << ", " 
        << zoneCells << ", " << cellZonesType[typeI].size() << endl;

    snapshotFileName = "k"+name(typeI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(cellZonesType[typeI], zoneI)
    {                   
        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << turbKinEnergy[cellN] << " ";
        }
        outputFilePtr() << endl;  
    } 
}


// output epsilon
volScalarField turbEpsilon = turbulence().epsilon();
forAll(typesName, typeI)
{
    label zoneCells(mesh.cellZones()[cellZonesType[typeI][0]].size());
    Info<< "Size of the snapshots, epsilon" << name(typeI) << ", " 
        << zoneCells << ", " << cellZonesType[typeI].size() << endl;

    snapshotFileName = "epsilon"+name(typeI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(cellZonesType[typeI], zoneI)
    {                   
        forAll(mesh.cellZones()[cellZonesType[typeI][zoneI]], cellI)
        {
            label cellN (mesh.cellZones()[cellZonesType[typeI][zoneI]][cellI]);
            outputFilePtr() << turbEpsilon[cellN] << " ";
        }
        outputFilePtr() << endl;  
    } 
}

// ************************************************************************* //

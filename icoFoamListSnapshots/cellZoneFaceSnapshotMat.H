// Rsh. get postprocessed field value (i.e. grad(U), laplacian(U) and Ux, ..., 
// --- P, grad(P)) in each cell zone
// -- and write them into a file each column of represent a cell zone


// create face field
surfaceScalarField pFace = fvc::interpolate(p);

// p, pressure face value
for(label patchI=0; patchI<patchNum; ++patchI)
{
    const label zoneFaces(globalfaceIndexListPtr[patchI][0].size());
    Info<< "Size of the face snapshots, p_" << patchI << ": " << zoneFaces << endl;

    fileName snapshotFileName = "pSnapshots_"+name(patchI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(globalfaceIndexListPtr[patchI], domainI)
    {                   
        forAll(globalfaceIndexListPtr[patchI][domainI], faceI)
        {
            label faceID (globalfaceIndexListPtr[patchI][domainI][faceI]);
            if(faceID < pFace.size())
                outputFilePtr() << pFace[faceID] << " ";
            else
            {
                forAll(mesh.boundary(), boundaryI)
                {
                    if(mesh.boundary()[boundaryI].start() <= faceID && 
                        faceID < mesh.boundary()[boundaryI].start()+mesh.boundary()[boundaryI].size())
                    {
                        label bounFaceID = faceID - mesh.boundary()[boundaryI].start();
                        outputFilePtr() << p.boundaryField()[boundaryI][bounFaceID] << " ";
                        break;
                    }
                }
            }
        }
        outputFilePtr() << endl;  
    }
}


// create face field
surfaceVectorField UFace = fvc::interpolate(U);

// U, veclocity face value
for(label patchI=0; patchI<patchNum; ++patchI)
{
    const label zoneFaces(globalfaceIndexListPtr[patchI][0].size());
    Info<< "Size of the face snapshots, U_" << patchI << ": " << 3*zoneFaces << endl;

    fileName snapshotFileName = "USnapshots_"+name(patchI);
    outputFilePtr.reset(new OFstream(runTime.globalPath()/outputDir/snapshotFileName));

    forAll(globalfaceIndexListPtr[patchI], domainI)
    {                   
        forAll(globalfaceIndexListPtr[patchI][domainI], faceI)
        {
            label faceID (globalfaceIndexListPtr[patchI][domainI][faceI]);
            if(faceID < UFace.size())
                outputFilePtr() << UFace[faceID].x() << " ";
            else
            {
                forAll(mesh.boundary(), boundaryI)
                {
                    if(mesh.boundary()[boundaryI].start() <= faceID && 
                        faceID < mesh.boundary()[boundaryI].start()+mesh.boundary()[boundaryI].size())
                    {
                        label bounFaceID = faceID - mesh.boundary()[boundaryI].start();
                        outputFilePtr() << U.boundaryField()[boundaryI][bounFaceID].x() << " ";
                        break;
                    }
                }
            }
        }

        forAll(globalfaceIndexListPtr[patchI][domainI], faceI)
        {
            label faceID (globalfaceIndexListPtr[patchI][domainI][faceI]);
            if(faceID < UFace.size())
                outputFilePtr() << UFace[faceID].y() << " ";
            else
            {
                forAll(mesh.boundary(), boundaryI)
                {
                    if(mesh.boundary()[boundaryI].start() <= faceID && 
                        faceID < mesh.boundary()[boundaryI].start()+mesh.boundary()[boundaryI].size())
                    {
                        label bounFaceID = faceID - mesh.boundary()[boundaryI].start();
                        outputFilePtr() << U.boundaryField()[boundaryI][bounFaceID].y() << " ";
                        break;
                    }
                }
            }
        }

        forAll(globalfaceIndexListPtr[patchI][domainI], faceI)
        {
            label faceID (globalfaceIndexListPtr[patchI][domainI][faceI]);
            if(faceID < UFace.size())
                outputFilePtr() << UFace[faceID].z() << " ";
            else
            {
                forAll(mesh.boundary(), boundaryI)
                {
                    if(mesh.boundary()[boundaryI].start() <= faceID && 
                        faceID < mesh.boundary()[boundaryI].start()+mesh.boundary()[boundaryI].size())
                    {
                        label bounFaceID = faceID - mesh.boundary()[boundaryI].start();
                        outputFilePtr() << U.boundaryField()[boundaryI][bounFaceID].z() << " ";
                        break;
                    }
                }
            }
        }


        outputFilePtr() << endl;  
    }
}

// ************************************************************************* //

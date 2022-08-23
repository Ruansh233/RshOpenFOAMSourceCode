
    List<fileName> modeNames (modesNum);
    List<fileName> laplacianModeNames (modesNum);
    List<fileName> gradModeNames (modesNum);
    List<fileName> gaussGradModeNames (modesNum);
    List<fileName> sfVectorNames (modesNum);


    for (label iname = 0; iname < modesNum; ++iname)
    {
        modeNames[iname] = "Tmode" + name(iname);
        laplacianModeNames[iname] = "Tmode_lapla" + name(iname);
        gradModeNames[iname] = "Tmode_grad" + name(iname);
        gaussGradModeNames[iname] = "Tmode_gaussGrad" + name(iname);
        sfVectorNames[iname] = "mesh_sf" + name(iname);
    }

    // matrix of gradient of mode 
    RectangularMatrix<scalar> laplacianModesMatrix(mesh.C().size(), modesNum);
    // wordList patchTypes_ ({"zeroGradient", "zeroGradient", "zeroGradient", "zeroGradient", 
    //                         "fixedValue", "fixedValue"});

    // calculate the laplacian of mode
    forAll(modeNames, No_)
    {
        // Field<scalar> zeroScalarField (mesh.C().size(), Foam::zero());
        volScalarField fieldValueMode
        (
            IOobject
            (
                modeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimless
            // patchTypes_
        );

        forAll(mesh.C(), cellI)
        {
            fieldValueMode[cellI] = fieldValue.U()[cellI][No_];
        }

        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     Info << mesh.boundary()[patchI].name() << ": " << fieldValueMode.boundaryField()[patchI][5] << endl;
        // }

        // // use to set boundary value of mdoes
        // label rows(mesh.C().size());

        // // set the fixedValue boundary value 
        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     if(fieldValueMode.boundaryField()[patchI].type() == "fixedValue")
        //     {
        //         forAll(fieldValueMode.boundaryField()[patchI], faceI)
        //         {
        //             fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValue.U()[rows][No_];
        //             ++rows;
        //         }
        //     }
        //     else
        //     {
        //         rows += fieldValueMode.boundaryField()[patchI].size();
        //     }
        // }

        // // set the face value in boundary to the value in snapshots
        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     forAll(fieldValueMode.boundaryField()[patchI], faceI)
        //     {
        //         fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValue.U()[rows][No_];
        //         ++rows;
        //     }
        // }

        // // fieldValueMode.boundaryFieldRef().set() is a function of Ptrlist 
        // // because boundaryFieldRef() is a pointer
        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     fieldValueMode.boundaryFieldRef().set(patchI, 
        //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
        // }

        // set the face value in boundary to the value in the face neighbouring cell
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            forAll(fieldValueMode.boundaryField()[patchI], faceI)
            {
                label gFaceI = mesh.boundary()[patchI].start() + faceI;
                label bunFaceOwnerI = mesh.faceOwner()[gFaceI];
                fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValueMode[bunFaceOwnerI];
            }
        }

        fieldValueMode.write();

        // laplacian of mode
        volScalarField lapFieldValueMode
        (
            IOobject
            (
                laplacianModeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            // mesh,
            // dimless,
            fvc::laplacian(fieldValueMode)
            // "zeroGradient"
        );
        lapFieldValueMode.write();   

        forAll(mesh.C(), cellI)
        {
            laplacianModesMatrix[cellI][No_] = lapFieldValueMode[cellI];
        }

        // forAll(lapFieldValueMode.boundaryField(), patchI)
        // {
        //     lapFieldValueMode.boundaryFieldRef().set(patchI, 
        //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], lapFieldValueMode));
        // }

        // Rsh, I don't know what the difference betwween geoField<IO, Field> and 
        // geoField<IO, mesh, dim,Field>? confused
        // gradient of mode
        volVectorField gradFieldValueMode
        (
            IOobject
            (
                gradModeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            // mesh,
            // dimless,
            fvc::grad(fieldValueMode)
        );

        gradFieldValueMode.write();   

        // // surfaceScalarField snGradMode (fvc::interpolate(fieldValueMode));
        // // snGradMode.write();

        // // surfaceScalarField snPhiMode (fvc::interpolate(fieldValueMode) * mesh.Sf().component(0)
        // //                             + fvc::interpolate(fieldValueMode) * mesh.Sf().component(1)
        // //                             + fvc::interpolate(fieldValueMode) * mesh.Sf().component(2));
        // // snPhiMode.write();

        // surfaceScalarField facesInterMode (fvc::interpolate(fieldValueMode));
        // facesInterMode.write();

        // // Info << "mesh.Sf(): \n" << mesh.Sf() << endl;      
        // // transfer the data as surfaceVectorField will result in wrong boundary value 
        // // surfaceVectorField facesVectorMode
        // // (
        // //     IOobject
        // //     (
        // //         sfVectorNames[No_],
        // //         mesh.time().timeName(),
        // //         mesh,
        // //         IOobject::NO_READ,
        // //         IOobject::AUTO_WRITE
        // //     ),
        // //     mesh,
        // //     dimless,
        // //     mesh.Sf()
        // // );
        // // facesVectorMode.write();
        // // Info << mesh.Sf() << endl;

        // // gauss gradient of mode with self defined method
        // volVectorField GuassGradFieldValueMode
        // (
        //     IOobject
        //     (
        //         gaussGradModeNames[No_],
        //         mesh.time().timeName(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimless
        // );

        // forAll(mesh.cells(), cellI)
        // {
        //     GuassGradFieldValueMode[cellI] = vector::zero;
        //     List<label> adjFaces = mesh.cells()[cellI];
        //     forAll(adjFaces, faceI)
        //     {
        //         scalar innerProductCf = (mesh.C()[cellI] - mesh.Cf()[adjFaces[faceI]]) & 
        //                                 mesh.Sf()[adjFaces[faceI]];

        //         if(innerProductCf < scalar(0))
        //         {
        //             GuassGradFieldValueMode[cellI] += 
        //                 facesInterMode[adjFaces[faceI]] * mesh.Sf()[adjFaces[faceI]];
        //         }
        //         else
        //         {
        //             GuassGradFieldValueMode[cellI] += 
        //                 facesInterMode[adjFaces[faceI]] * - mesh.Sf()[adjFaces[faceI]];
        //         }                
        //     }

        //     GuassGradFieldValueMode[cellI] /= mesh.V()[cellI];
        // }      

        // GuassGradFieldValueMode.write();  




    }


    List<fileName> modeNames (modesNum);
    List<fileName> laplacianModeNames (modesNum);
    List<fileName> gradModeNames (modesNum);


    for (label iname = 0; iname < modesNum; ++iname)
    {
        modeNames[iname] = "Tmode" + name(iname);
        laplacianModeNames[iname] = "Tmode_lapla" + name(iname);
        gradModeNames[iname] = "Tmode_grad" + name(iname);
    }

    // matrix of gradient of mode 
    RectangularMatrix<scalar> laplacianModesMatrix(snapshotsRows, modesNum);
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

        // use to set boundary value of mdoes
        label rows(mesh.C().size());

        // // set the face value in boundary to fixedvalue
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

        // set the face value in boundary to the value in snapshots
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            forAll(fieldValueMode.boundaryField()[patchI], faceI)
            {
                fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValue.U()[rows][No_];
                ++rows;
            }
        }

        // // fieldValueMode.boundaryFieldRef().set() is a function of Ptrlist 
        // // because boundaryFieldRef() is a pointer
        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     fieldValueMode.boundaryFieldRef().set(patchI, 
        //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
        // }

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
            mesh,
            dimless,
            fvc::laplacian(fieldValueMode)
            // "zeroGradient"
        );
        lapFieldValueMode.write();   

        // forAll(mesh.C(), cellI)
        // {
        //     laplacianModesMatrix[cellI][No_] = lapFieldValueMode[cellI];
        // }

        // forAll(lapFieldValueMode.boundaryField(), patchI)
        // {
        //     lapFieldValueMode.boundaryFieldRef().set(patchI, 
        //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], lapFieldValueMode));
        // }

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
            mesh,
            dimless,
            fvc::grad(fieldValueMode)
        );

        gradFieldValueMode.write();   

        // surfaceScalarField snGradMode (fvc::interpolate(fieldValueMode));
        // snGradMode.write();

        // surfaceScalarField snPhiMode (fvc::interpolate(fieldValueMode) * mesh.Sf().component(0)
        //                             + fvc::interpolate(fieldValueMode) * mesh.Sf().component(1)
        //                             + fvc::interpolate(fieldValueMode) * mesh.Sf().component(2));
        // snPhiMode.write();

        surfaceVectorField snPhiMode (fvc::interpolate(fieldValueMode) * mesh.Sf());
        snPhiMode.write();

        const labelListList & cellFaces = mesh.cellFaces();

        
    }

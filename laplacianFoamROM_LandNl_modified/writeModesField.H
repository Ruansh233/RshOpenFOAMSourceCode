
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
    RectangularMatrix<scalar> laplacianModesMatrix(mesh.C().size(), modesNum);
    RectangularMatrix<vector> gradModesMatrix(mesh.C().size(), modesNum);

    // forAll(T.boundaryField(), patchI)
    // {
    //     // for fixedValue, // phibc = Tbc * V * inv(S)
    //     if(T.boundaryField()[patchI].type() == "fixedValue")
    //     {
    //         RectangularMatrix<scalar> boundaryValueM(T.boundaryField()[patchI].size(), 
    //                                                 snapshotsNo, 
    //                                                 T.boundaryField()[patchI][0]);

    //         dataFile = mesh.time().path()/"SVD"/"boundaryValueM";
    //         writeMatrix(boundaryValueM, dataFile);
    //         Info << "write boundaryValueM Matrix: " << boundaryValueM.sizes() << endl;

    //         RectangularMatrix<scalar> boundaryModeM(boundaryValueM * invCoefficientsM);
    //         dataFile = mesh.time().path()/"SVD"/"boundaryModeM";
    //         writeMatrix(boundaryModeM, dataFile);
    //         Info << "write boundaryModeM Matrix: " << boundaryModeM.sizes() << endl;
    //     }
    // }

    PtrList<RectangularMatrix<scalar>> boundaryModesMList;

    forAll(T.boundaryField(), patchI)
    {
        // for fixedValue, // phibc = Tbc * V * inv(S)
        if(T.boundaryField()[patchI].type() == "fixedValue")
        {
            RectangularMatrix<scalar> boundaryValueM(T.boundaryField()[patchI].size(), 
                                                    snapshotsNo, 
                                                    T.boundaryField()[patchI][0]);

            boundaryModesMList.append(RectangularMatrix<scalar> (boundaryValueM * invCoefficientsM).clone());
        }

        // for fixedGradient, the face value is equal to the face owner
        else if(T.boundaryField()[patchI].type() == "fixedGradient")
        {
            // the steps to access the gradient of patch
            // the following line is usable or not
            // fixedGradientFvPatchField<scalar> (T.boundaryField()[patchI]).gradient
            // or T.boundaryField()[patchI].snGrad()
            // Define new patch-field for the variable at the patch (in my case temperature T)
            fixedGradientFvPatchField<scalar>& Tpatch = refCast<fixedGradientFvPatchField<scalar>>(T.boundaryFieldRef()[patchI]);
            // Create scalar Field for the gradient at that patch
            scalarField& gradTpatch = Tpatch.gradient();
            
            // phibc = Tbc * V * inv(S)
            RectangularMatrix<scalar> boundaryValueM(T.boundaryField()[patchI].size(), 
                                                        snapshotsNo, 
                                                        gradTpatch[0]);

            boundaryModesMList.append(RectangularMatrix<scalar> (boundaryValueM * invCoefficientsM).clone());
        }

        // for other type, defined as zeroGradient
        else
        {
            boundaryModesMList.append(RectangularMatrix<scalar> (T.boundaryField()[patchI].size(), 
                                                    snapshotsNo, 0).clone());
        }
    }

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
            fieldValueMode[cellI] = fieldValueSVD.U()[cellI][No_];
        }

        // assign boundary value of all boundaries
        // creating the mode field by copying the T, than the setting of patch type is not needed
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            // for fixedValue, // phibc = Tbc * V * inv(S)
            if(T.boundaryField()[patchI].type() == "fixedValue")
            {
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patchI], fieldValueMode));

                forAll(fieldValueMode.boundaryField()[patchI], faceI)
                {
                    fieldValueMode.boundaryFieldRef()[patchI][faceI] = boundaryModesMList[patchI](faceI, No_);
                }
            }

            // for fixedGradient, the face value is equal to the face owner
            else if(T.boundaryField()[patchI].type() == "fixedGradient")
            {
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<scalar>::New("fixedGradient", mesh.boundary()[patchI], fieldValueMode));

                // Define new patch-field for the variable at the patch (in my case temperature T)
                fixedGradientFvPatchField<scalar>& Mpatch = refCast<fixedGradientFvPatchField<scalar>>(fieldValueMode.boundaryFieldRef()[patchI]);
                // Create scalar Field for the gradient at that patch
                scalarField& gradMpatch = Mpatch.gradient();

                // Loop over every gradient in this patch
                forAll(gradMpatch, faceI)
                {
                    // Perform calculations to determine the new gradient "newGrad"
                    // Store the new value in the scalar field
                    gradMpatch[faceI] = boundaryModesMList[patchI](faceI, No_);
                }
            }

            // for zeroGradient
            else if(T.boundaryField()[patchI].type() == "zeroGradient")
            {
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
            }

            else
            {
                Info<< "The boundary: " << T.boundaryField()[patchI].type() 
                    << " doesn't support, and redefine to zeroGradient" << nl;
                    
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
            }
        }
        fieldValueMode.write();

        Info<< "field of mode " << No_ << " is written" << endl;

        // gradient of mode
        volScalarField fieldValueModeRef
        (
            IOobject
            (
                modeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        
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
            fvc::grad(fieldValueModeRef)
        );
        gradFieldValueMode.write();   

        forAll(mesh.C(), cellI)
        {
            gradModesMatrix[cellI][No_] = gradFieldValueMode[cellI];
        }

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
            fvc::laplacian(fieldValueModeRef)
        );
        lapFieldValueMode.write();   

        forAll(mesh.C(), cellI)
        {
            laplacianModesMatrix[cellI][No_] = lapFieldValueMode[cellI];
        }
    }

    

    // // Rsh, write the file and then open it, the value near the boundary face will be correct
    // // T field test
    // volScalarField testT
    // (
    //     IOobject
    //     (
    //         "T",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // volVectorField testgradT
    // (
    //     IOobject
    //     (
    //         "testgradT",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::grad(T)
    // );
    // testgradT.write();

    // volScalarField testlapT
    // (
    //     IOobject
    //     (
    //         "testlapdT",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::laplacian(T)
    // );
    // testlapT.write();


    // // mode field test 
    // volScalarField testMode
    // (
    //     IOobject
    //     (
    //         "Tmode0",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh
    // );

    // volVectorField testgradMode
    // (
    //     IOobject
    //     (
    //         "testgradMode",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::grad(testMode)
    // );
    // testgradMode.write();

    // volScalarField testlapMode
    // (
    //     IOobject
    //     (
    //         "testlapMode",
    //         runTime.timeName(),
    //         mesh,
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     fvc::laplacian(testMode)
    // );
    // testlapMode.write();
    

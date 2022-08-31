
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

    Info << "test3 " << nl;


    // matrix of gradient of mode 
    RectangularMatrix<scalar> laplacianModesMatrix(mesh.C().size(), modesNum);

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

        scalar eigenValueSqrt (fieldValueSVD.S()[No_]);

        // assign boundary value of all boundaries
        forAll(fieldValueMode.boundaryField(), patchI)
        {
            // for fixedValue, // phibc = Tbc * V * inv(S)
            if(T.boundaryField()[patchI].type() == "fixedValue")
            {
                fieldValueMode.boundaryFieldRef().set(patchI, 
                    fvPatchField<scalar>::New("fixedValue", mesh.boundary()[patchI], fieldValueMode));

                // phibc = Tbc * V * inv(S)
                RectangularMatrix<scalar> boundaryValueM(fieldValueMode.boundaryField()[patchI].size(), 
                                                         snapshotsNo, 
                                                         T.boundaryField()[patchI][0]);
                RectangularMatrix<scalar> boundaryModeM(boundaryValueM * 
                                RectangularMatrix<scalar>(fieldValueSVD.V().subColumn(No_, 0, -1)));
                
                for(label i = 0; i < boundaryModeM.m(); i++)
                {
                    for(label j = 0; j < boundaryModeM.n(); j++)
                    {
                        boundaryModeM(i, j) = boundaryModeM(i, j) / eigenValueSqrt;
                    }
                }

                forAll(fieldValueMode.boundaryField()[patchI], faceI)
                {
                    fieldValueMode.boundaryFieldRef()[patchI][faceI] = boundaryModeM(faceI, 0);
                }
            }
            // for zeroGradient, the face value is equal to the face owner
            else if(T.boundaryField()[patchI].type() == "zeroGradient")
            {
                // fieldValueMode.boundaryFieldRef().set(patchI, 
                //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], fieldValueMode));
                forAll(fieldValueMode.boundaryField()[patchI], faceI)
                {
                    label gFaceI = mesh.boundary()[patchI].start() + faceI;
                    label bunFaceOwnerI = mesh.faceOwner()[gFaceI];
                    fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValueMode[bunFaceOwnerI];
                }
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

        // // this kind of definition also give strange value on boundary cell
        // // -- do not use
        // volVectorField gradFieldValueMode
        // (
        //     fvc::grad(fieldValueMode)
        // );
        
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
            fvc::grad(fieldValueMode)
        );
        gradFieldValueMode.write();   

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
    }

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
    


    List<fileName> modeNames (modesNum);
    List<fileName> laplacianModeNames (modesNum);

    for (label iname = 0; iname < modesNum; ++iname)
    {
        modeNames[iname] = "Tmode" + name(iname);
        laplacianModeNames[iname] = "Tmode_lapla" + name(iname);
    }

    // matrix of gradient of mode 
    RectangularMatrix<scalar> laplacianModesMatrix(mesh.C().size(), modesNum);
    wordList patchTypes_ ({"zeroGradient", "zeroGradient", "zeroGradient", "zeroGradient", 
                            "extrapolatedCalculated", "extrapolatedCalculated"});

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
            dimless,
            patchTypes_
        );

        forAll(mesh.C(), cellI)
        {
            fieldValueMode[cellI] = fieldValue.U()[cellI][No_];
        }

        // // set the face value in boundary equal to the neibouring cell
        // forAll(fieldValueMode.boundaryField(), patchI)
        // {
        //     forAll(fieldValueMode.boundaryField()[patchI], faceI)
        //     {
        //         label adjcellID = mesh.boundary()[patchI].patch().faceCells()[faceI];

        //         fieldValueMode.boundaryFieldRef()[patchI][faceI] = fieldValueMode[adjcellID];
        //     }
        // }

        // // fieldValueMode.boundaryFieldRef().set() is
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

        // forAll(mesh.C(), cellI)
        // {
        //     laplacianModesMatrix[cellI][No_] = lapFieldValueMode[cellI];
        // }

        // forAll(lapFieldValueMode.boundaryField(), patchI)
        // {
        //     lapFieldValueMode.boundaryFieldRef().set(patchI, 
        //     fvPatchField<scalar>::New("zeroGradient", mesh.boundary()[patchI], lapFieldValueMode));
        // }

        // // gradient of mode
        // volVectorField lapFieldValueMode
        // (
        //     IOobject
        //     (
        //         laplacianModeNames[No_],
        //         mesh.time().timeName(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     mesh,
        //     dimless,
        //     fvc::grad(fieldValueMode)
        // );

        lapFieldValueMode.write();   
    }

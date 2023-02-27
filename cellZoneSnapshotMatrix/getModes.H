void getModes(const fvMesh & mesh_,
            PtrList<volScalarField> tmpscalarFields_,
            word fieldName_,
            label modeNumbers_,
            List<List<label>> cellZonesList_,
            List<List<label>> zonePatchID_,
            List<word> boundaryNames_
            )
{
    label fieldMatrixRows(mesh_.cellZones()[cellZonesList_[0][0]].size());
    label fieldMatrixColumns(0);
    forAll(cellZonesList_, caseI)
    {
        fieldMatrixColumns += cellZonesList_[caseI].size(); 
    }

    RectangularMatrix<scalar> fieldMatrix(fieldMatrixRows, fieldMatrixColumns);

    label totalZoneNo(0);
    forAll(cellZonesList_, caseI)
    {
        forAll(cellZonesList_[caseI], zoneI)
        {
            forAll(mesh_.cellZones()[cellZonesList_[caseI][zoneI]], cellI)
            {
                label cellN (mesh_.cellZones()[cellZonesList_[caseI][zoneI]][cellI]);
                fieldMatrix(cellI, totalZoneNo) = tmpscalarFields_[caseI][cellN];
            }
            ++totalZoneNo;
        }
    }

    Info<< "fieldMatrix.sizes(): " << fieldMatrix.sizes() << endl;

    fileName dataFile;
    dataFile = mesh_.time().globalPath()/"SVD"/"fieldMatrix";
    writeMatrix(fieldMatrix, dataFile);

    // SVD fieldValueSVD(fieldMatrix); 

    // // SVDinv to calculate inverse of coefficient, but will have problem if 1/eigenValue is very large
    // RectangularMatrix<scalar> coefficientsM(fieldMatrix.T() * fieldValueSVD.U());
    // RectangularMatrix<scalar> invCoefficientsM(SVDinv(coefficientsM.T()));
    // // write snapshots Matrix
    // fileName dataFile;
    // dataFile = mesh_.time().globalPath()/"SVD"/"coefficientsM";
    // writeMatrix(coefficientsM, dataFile);
    // dataFile = mesh_.time().globalPath()/"SVD"/"invCoefficientsM";
    // writeMatrix(SVDinv(coefficientsM), dataFile);

    // dataFile = mesh_.time().globalPath()/"SVD"/"fieldValueModes";
    // writeMatrix(fieldValueSVD.U(), dataFile);
    // dataFile = mesh_.time().globalPath()/"SVD"/"eigenValue";
    // writeMatrix(fieldValueSVD.S(), dataFile);
    // dataFile = mesh_.time().globalPath()/"SVD"/"fieldValueV";
    // writeMatrix(fieldValueSVD.V(), dataFile);

    // Info << "fieldValueSVD.V(): " << fieldValueSVD.V().sizes() << endl;
    // RectangularMatrix<scalar> coefficientsMTest(fieldValueSVD.V().m(), fieldValueSVD.V().n());
    // for(label i=0; i<coefficientsMTest.m(); ++i)
    // {
    //     for(label j=0; j<coefficientsMTest.n(); ++j)
    //     {
    //         coefficientsMTest(i, j) = fieldValueSVD.V()(i, j) * fieldValueSVD.S()[i];
    //     }
    // }
    // dataFile = mesh_.time().globalPath()/"SVD"/"coefficientsMTest";
    // writeMatrix(coefficientsMTest, dataFile);


    // RectangularMatrix<scalar> coefficientsMTest2(fieldValueSVD.V().m(), fieldValueSVD.V().n());
    // for(label i=0; i<coefficientsMTest2.m(); ++i)
    // {
    //     for(label j=0; j<coefficientsMTest2.n(); ++j)
    //     {
    //         coefficientsMTest2(i, j) = fieldValueSVD.V()(i, j) * 1/fieldValueSVD.S()[i];
    //     }
    // }
    // dataFile = mesh_.time().globalPath()/"SVD"/"coefficientsMTest2";
    // writeMatrix(coefficientsMTest2, dataFile);

    // RectangularMatrix<scalar> mulTest(coefficientsMTest2 * coefficientsMTest.T());
    // dataFile = mesh_.time().globalPath()/"SVD"/"mulTest";
    // writeMatrix(mulTest, dataFile);

    // PtrList<fvPatchField<scalar>> boundarySnapshots;
    // PtrList<RectangularMatrix<scalar>> boundaryModesMatList;

    // forAll(boundaryNames_, bounI)
    // {
    //     wordRe patchMatch(".*"+boundaryNames_[bounI]);
    //     patchMatch.compile();

    //     forAll(cellZonesList_, caseI)
    //     {
    //         forAll(cellZonesList_[caseI], zoneI)
    //         {
    //             List<label> patchIDList (zonePatchID_[cellZonesList_[caseI][zoneI]]);
    //             forAll(patchIDList, patchI)
    //             {
    //                 if(patchMatch.match(mesh_.boundary()[patchIDList[patchI]].name()))
    //                 {
    //                     boundarySnapshots.append(tmpscalarFields_[caseI].boundaryFieldRef()[patchIDList[patchI]].clone());
    //                     break;
    //                 }
    //             }
    //         }
    //     }

    //     RectangularMatrix<scalar> boundarySnapshotsMat(boundarySnapshots[0].size(), boundarySnapshots.size());

    //     forAll(boundarySnapshots, colunmI)
    //     {
    //         forAll(boundarySnapshots[colunmI], rowI)
    //         {
    //             boundarySnapshotsMat(rowI, colunmI) = boundarySnapshots[colunmI][rowI];
    //         }
    //     }

    //     boundaryModesMatList.append(RectangularMatrix<scalar> (boundarySnapshotsMat * invCoefficientsM).clone());
    //     boundarySnapshots.clear();
    // }    


    // List<fileName> modeNames (modeNumbers_);
    // forAll(modeNames, nameI)
    // {
    //     modeNames[nameI] = fieldName_ + name(nameI + 1);
    // }


    // // This can not write to file
    // // List<List<List<scalar>>> scalarModesField(modeNumbers_);
    // // // set mode field value for scalar field
    // // forAll(modeNames, No_)
    // // {
    // //     scalarModesField[No_].resize(boundaryNames_.size()+1);

    // //     scalarModesField[No_][0].resize(fieldValueSVD.U().m()); 
    // //     forAll(scalarModesField[No_][0], cellI)
    // //     {
    // //         scalarModesField[No_][0][cellI] = fieldValueSVD.U()(cellI, No_);
    // //     }

    // //     forAll(boundaryNames_, bounI)
    // //     {
    // //         scalarModesField[No_][bounI+1].resize(boundaryModesMatList[bounI].m());
    // //         forAll(scalarModesField[No_][bounI+1], faceI)
    // //         {
    // //             scalarModesField[No_][bounI+1][faceI] = boundaryModesMatList[bounI](faceI, No_);
    // //         }
    // //     }
    // // }

    // // for the list, the 0 sublist is cell value, the rest sub-list is face value
    // List<List<scalar>> scalarModesField(boundaryNames_.size()+1);
    // // set mode field value for scalar field
    // forAll(modeNames, No_)
    // {
    //     forAll(scalarModesField, listI)
    //     {
    //         scalarModesField[0].resize(fieldValueSVD.U().m()); 
    //         forAll(scalarModesField[0], cellI)
    //         {
    //             scalarModesField[0][cellI] = fieldValueSVD.U()(cellI, No_);
    //         }

    //         forAll(boundaryNames_, bounI)
    //         {
    //             scalarModesField[bounI+1].resize(boundaryModesMatList[bounI].m());
    //             forAll(scalarModesField[bounI+1], faceI)
    //             {
    //                 scalarModesField[bounI+1][faceI] = boundaryModesMatList[bounI](faceI, No_);
    //             }
    //         } 
    //     }

    //     // IOList for list write
    //     IOList<List<scalar>> scalarModesFieldIO
    //     (
    //         IOobject
    //         (
    //             modeNames[No_],
    //             mesh_.time().constant(),
    //             mesh_,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         scalarModesField
    //     );

    //     scalarModesFieldIO.write();  
    // }
}

// void getModes(const fvMesh & mesh_,
//             PtrList<volVectorField> tmpvectorFields_,
//             word fieldName_,
//             label modeNumbers_,
//             List<List<label>> cellZonesList_,
//             List<List<label>> zonePatchID_,
//             List<word> boundaryNames_
//             )

// {
//     label elementCellSize(mesh_.cellZones()[cellZonesList_[0][0]].size());
//     label fieldMatrixRows(3*elementCellSize);
//     label fieldMatrixColumns(0);
//     forAll(cellZonesList_, caseI)
//         fieldMatrixColumns += cellZonesList_[caseI].size();

//     RectangularMatrix<scalar> fieldMatrix(fieldMatrixRows, fieldMatrixColumns);

//     forAll(cellZonesList_, caseI)
//     {
//         forAll(cellZonesList_[caseI], zoneI)
//         {
//             forAll(mesh_.cellZones()[cellZonesList_[caseI][zoneI]], cellI)
//             {
//                 label cellN (mesh_.cellZones()[cellZonesList_[caseI][zoneI]][cellI]);
//                 fieldMatrix(cellI,                   zoneI+caseI*cellZonesList_[caseI].size()) = tmpvectorFields_[caseI][cellN].x();
//                 fieldMatrix(cellI+elementCellSize,   zoneI+caseI*cellZonesList_[caseI].size()) = tmpvectorFields_[caseI][cellN].y();
//                 fieldMatrix(cellI+2*elementCellSize, zoneI+caseI*cellZonesList_[caseI].size()) = tmpvectorFields_[caseI][cellN].z();
//             }
//         }
//     }

//     SVD fieldValueSVD(fieldMatrix);
//     RectangularMatrix<scalar> invCoefficientsM(SVDinv(fieldValueSVD.U().T() * fieldMatrix));


//     PtrList<fvPatchField<vector>> boundarySnapshots;
//     PtrList<RectangularMatrix<scalar>> boundaryModesMatList;

//     forAll(boundaryNames_, bounI)
//     {
//         wordRe patchMatch(".*"+boundaryNames_[bounI]);
//         patchMatch.compile();

//         forAll(cellZonesList_, caseI)
//         {
//             forAll(cellZonesList_[caseI], zoneI)
//             {
//                 List<label> patchIDList (zonePatchID_[cellZonesList_[caseI][zoneI]]);
//                 forAll(patchIDList, patchI)
//                 {
//                     if(patchMatch.match(mesh_.boundary()[patchIDList[patchI]].name()))
//                     {
//                         boundarySnapshots.append(tmpvectorFields_[caseI].boundaryFieldRef()[patchIDList[patchI]].clone());
//                         break;
//                     }
//                 }
//             }
//         }

//         label patchFaceNumber(boundarySnapshots[0].size());
//         RectangularMatrix<scalar> boundarySnapshotsMat(3*patchFaceNumber, boundarySnapshots.size());

//         forAll(boundarySnapshots, colunmI)
//         {
//             forAll(boundarySnapshots[colunmI], rowI)
//             {
//                 boundarySnapshotsMat(rowI,                   colunmI) = boundarySnapshots[colunmI][rowI].x();
//                 boundarySnapshotsMat(rowI+patchFaceNumber,   colunmI) = boundarySnapshots[colunmI][rowI].y();
//                 boundarySnapshotsMat(rowI+2*patchFaceNumber, colunmI) = boundarySnapshots[colunmI][rowI].z();
//             }
//         }

//         boundaryModesMatList.append(RectangularMatrix<scalar> (boundarySnapshotsMat * invCoefficientsM).clone());
//         boundarySnapshots.clear();
//     }    


//     List<fileName> modeNames (modeNumbers_);
//     forAll(modeNames, nameI)
//     {
//         modeNames[nameI] = fieldName_ + name(nameI + 1);
//     }


//     List<List<vector>> vectorModesField(modeNumbers_);
//     // set mode field value for vector field
//     forAll(modeNames, No_)
//     {
//         forAll(vectorModesField, listI)
//         {
//             vectorModesField[0].resize(elementCellSize); 
//             forAll(vectorModesField[0], cellI)
//             {
//                 vectorModesField[0][cellI] = Vector<scalar> (fieldValueSVD.U()(cellI,                     No_),
//                                                                 fieldValueSVD.U()(cellI + elementCellSize,   No_),
//                                                                 fieldValueSVD.U()(cellI + 2*elementCellSize, No_));;
//             }

//             forAll(boundaryNames_, bounI)
//             {
//                 label patchTypeISize(boundaryModesMatList[bounI].m()/3);
//                 vectorModesField[bounI+1].resize(patchTypeISize);
//                 forAll(vectorModesField[bounI+1], faceI)
//                 {
//                     vectorModesField[bounI+1][faceI] = Vector<scalar> (
//                         boundaryModesMatList[bounI](faceI,                      No_),
//                         boundaryModesMatList[bounI](faceI + patchTypeISize,     No_),
//                         boundaryModesMatList[bounI](faceI + 2 * patchTypeISize, No_));
//                 }
//             }
//         }

        
//         // IOList for list write
//         IOList<List<vector>> vectorModesFieldIO
//         (
//             IOobject
//             (
//                 modeNames[No_],
//                 mesh_.time().path()/"postProcessing",
//                 mesh_,
//                 IOobject::NO_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             vectorModesField
//         );

//         vectorModesFieldIO.write();   
//     }
// }

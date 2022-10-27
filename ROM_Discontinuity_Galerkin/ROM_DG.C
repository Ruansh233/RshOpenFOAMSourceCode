/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SVD.H"
#include "cpuTimeCxx.H"
#include "writeMatrix.H"
#include "wordRe.H"
#include "IFstream.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    fileName dataPath (mesh.time().path()/"SVD");
    label modesNum(5);

    // Discontinuity galerkin method parameters, add to dict reading later
    scalar epsilonPara(1.0);
    scalar xigema0(1.0);
    scalar xigema1(1.0);
    scalar beta0(1.0);
    scalar beta1(1.0);

    // element volume
    scalar elementVol(gSum(mesh.V()));


    // read the matrix
    // read cell modes matrix
    fileName dataFile (dataPath/"modesM");
    RectangularMatrix<scalar> modesM(mesh.C().size(), modesNum);
    RectangularMatrix<vector> gradModesM(mesh.C().size(), modesNum);

    if(isFile(dataFile))
    {                
        IFstream dataStream(dataFile);
        word dataLine;
        label row(0);

        while(dataStream.getLine(dataLine) && dataLine != word::null)
        {
            IStringStream dataString (dataLine);
            token singleData;  // token stores the data read from IFstream 

            for(label modesNo = 0; modesNo < modesNum; ++modesNo)
            {
                dataString.read(singleData);    
                modesM(row, modesNo) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }
    // dataFile = mesh.time().path()/"SVD"/"testModes";
    // writeMatrix(modesM, dataFile);

    // read boundary patch modes matrix and matchPatchID
    PtrList<RectangularMatrix<scalar>> boundaryModesMList;
    PtrList<RectangularMatrix<vector>> gradBoundaryModesMList;
    List<fileName> boundaryModesName({dataPath/"boundaryModes0", dataPath/"boundaryModes1"});
    RectangularMatrix<scalar> boundaryModesM(mesh.C().size(), modesNum);            

    forAll(boundaryModesName, patchI)
    {
        dataFile = boundaryModesName[patchI];
        label row(0);

        if(isFile(dataFile))
        {                
            IFstream dataStream(dataFile);
            word dataLine;
            while(dataStream.getLine(dataLine) && dataLine != word::null)
            {
                IStringStream dataString (dataLine);
                token singleData;  // token stores the data read from IFstream 

                for(label modesNo = 0; modesNo < modesNum; ++modesNo)
                {
                    dataString.read(singleData);    
                    boundaryModesM(row, modesNo) = singleData.scalarToken();
                }   
                ++row;
            }                       
        }  
        else
        {
            Info << "file: " << dataFile << " is not exist!" << endl;
            // break;
        }
        boundaryModesM.resize(row, modesNum);
        boundaryModesMList.append(boundaryModesM.clone());
    }
    // forAll(boundaryModesMList, matrixI)
    // {
    //     dataFile = mesh.time().path()/"SVD"/"test" + name(matrixI);
    //     writeMatrix(boundaryModesMList[matrixI], dataFile);
    // }

    // create field value
    List<fileName> modeNames (modesNum);
    for (label iname = 0; iname < modesNum; ++iname)
    {
        modeNames[iname] = "Tmode" + name(iname + 1);
    }

    // set patch value
    List<word> boundaryWordRe ({".*_in", ".*_out"});
    List<label> matchPatchID(mesh.boundaryMesh().size());
    wordRe matchPatch;            
    label countNumber(0);

    forAll(boundaryWordRe, reI)
    {
        matchPatch = boundaryWordRe[reI];
        matchPatch.compile();
        forAll(mesh.boundaryMesh(), patchI)
        {
            if(matchPatch.match(mesh.boundary()[patchI].name()))
            {
                matchPatchID[countNumber] = patchI;
                ++countNumber;
                break;
            }
        }
    }
    matchPatchID.resize(countNumber);

    // ptrlist to field value
    PtrList<volScalarField> fieldModesList;
    PtrList<volVectorField> gradfieldModesList;
    
    forAll(modeNames, No_)
    {
        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        
        // create mode field by copying T
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
            T
        );

        // set cell value
        forAll(mesh.C(), cellI)
        {
            fieldValueMode[cellI] = modesM(cellI, No_);
        }

        // set patch value
        forAll(matchPatchID, patchID)
        {
            forAll(fieldValueMode.boundaryField()[matchPatchID[patchID]], faceI)
            {
                fieldValueMode.boundaryFieldRef()[matchPatchID[patchID]][faceI] = boundaryModesMList[patchID](faceI, No_);
            }
        }

        fieldValueMode.write();
        fieldModesList.append(fieldValueMode.clone());
    }

    // the ptrlist for matchpatch matrix 
    gradBoundaryModesMList.append(gradModesM.clone());
    gradBoundaryModesMList.append(gradModesM.clone());
    gradBoundaryModesMList[0].resize(mesh.boundary()[matchPatchID[0]].size(), modesNum);
    gradBoundaryModesMList[1].resize(mesh.boundary()[matchPatchID[1]].size(), modesNum);

    // calculate gradModesM and gradBoundaryModesM
    forAll(modeNames, No_)
    {
        // create mode field by copying T
        volScalarField fieldValueModeTest
        (
            IOobject
            (
                modeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        // create mode field by copying T
        volVectorField fieldValueModegrad
        (
            IOobject
            (
                "grad" + modeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(fieldValueModeTest)
        );
        fieldValueModegrad.write();
        gradfieldModesList.append(fieldValueModegrad.clone());

        forAll(mesh.C(), cellI)
        {
            gradModesM[cellI][No_] = fieldValueModegrad[cellI];
        }

        forAll(matchPatchID, patchID)
        {
            forAll(fieldValueModegrad.boundaryField()[matchPatchID[patchID]], faceI)
            {
                gradBoundaryModesMList[patchID](faceI, No_) = 
                    fieldValueModegrad.boundaryField()[matchPatchID[patchID]][faceI];
            }
        }

        // // create mode field by copying T
        // volScalarField fieldValueModelap
        // (
        //     IOobject
        //     (
        //         "lap" + modeNames[No_],
        //         mesh.time().timeName(),
        //         mesh,
        //         IOobject::NO_READ,
        //         IOobject::AUTO_WRITE
        //     ),
        //     fvc::laplacian(fieldValueModeTest)
        // );
        // fieldValueModelap.write();
    }
    
    volScalarField fieldValueModeTest
    (
        fieldModesList[0] * fieldModesList[1]
    );

    // create Matrix system
    // initial global matrix
    RectangularMatrix<scalar> globalphiMatrix(modesNum * modesNum, modesNum * modesNum, Foam::Zero);
    RectangularMatrix<scalar> globalFMmatrix(modesNum * modesNum, 1, Foam::Zero);

    List<scalar> tempList;
    tempList.resize(mesh.C().size());

    // volumtric contribution
    RectangularMatrix<scalar> localphiMatrix(modesNum, modesNum, Foam::Zero);
    // localphiMatrix = gradModesM.T() & gradModesM;
    for (label row = 0; row < localphiMatrix.m(); ++row)
    {
        for (label column = 0; column < localphiMatrix.n(); ++column)
        {
            forAll(tempList, ListI)
            {
                tempList[ListI] = gradModesM(ListI, row) & gradModesM(ListI, column);
            }
            localphiMatrix(row, column) = gSum(tempList);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"localphiMatrix";
    writeMatrix(localphiMatrix, dataFile);

    // volumtric contribution
    RectangularMatrix<scalar> localphiMatrixTest(modesNum, modesNum, Foam::Zero);
    // localphiMatrix = gradModesM.T() & gradModesM;
    for (label row = 0; row < localphiMatrixTest.m(); ++row)
    {
        for (label column = 0; column < localphiMatrixTest.n(); ++column)
        {
            // localphiMatrix(row, column) = gSum(fieldModesList[row] * fieldModesList[column]);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"localphiMatrixTest";
    writeMatrix(localphiMatrixTest, dataFile);

    // // interface contribution
    // vector faceNormal(0, 0, 1);

    // // M11
    // RectangularMatrix<scalar> localInFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < localInFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < localInFaceMatrix.n(); ++column)
    //     {
    //         localInFaceMatrix(row, column) = -1/2 * (gradBoundaryModesMList[0](column, row) & faceNormal) 
    //                                               * boundaryModesMList[0](row, column)
    //                                          + 1/2 * epsilonPara * boundaryModesMList[0](column, row) 
    //                                               * (gradBoundaryModesMList[0](row, column) & faceNormal)
    //                                          + xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[0](column, row) 
    //                                               * boundaryModesMList[0](row, column)
    //                                          + xigema1/Foam::pow(elementVol, beta1) * (gradBoundaryModesMList[0](column, row) 
    //                                               & gradBoundaryModesMList[0](row, column));
    //     }
    // }

    // // M22
    // RectangularMatrix<scalar> localOutFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < localOutFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < localOutFaceMatrix.n(); ++column)
    //     {
    //         localOutFaceMatrix(row, column) =  1/2 * (gradBoundaryModesMList[1](column, row) & faceNormal)
    //                                                * boundaryModesMList[1](row, column)
    //                                           -1/2 * epsilonPara * boundaryModesMList[1](column, row) 
    //                                                * (gradBoundaryModesMList[1](row, column) & faceNormal)
    //                                           +xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[1](column, row) 
    //                                                * boundaryModesMList[1](row, column)
    //                                           +xigema1/Foam::pow(elementVol, beta1) * (gradBoundaryModesMList[1](column, row) 
    //                                                & gradBoundaryModesMList[1](row, column));
    //     }
    // }

    // // M12
    // RectangularMatrix<scalar> localIOFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < localIOFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < localIOFaceMatrix.n(); ++column)
    //     {
    //         localIOFaceMatrix(row, column) = -1/2 * (gradBoundaryModesMList[1](column, row) & faceNormal) 
    //                                               * boundaryModesMList[0](row, column)
    //                                          -1/2 * epsilonPara * boundaryModesMList[1](column, row) 
    //                                               * (gradBoundaryModesMList[0](row, column) & faceNormal)
    //                                          -xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[1](column, row) 
    //                                               * boundaryModesMList[0](row, column)
    //                                          -xigema1/Foam::pow(elementVol, beta1) * (gradBoundaryModesMList[1](column, row) 
    //                                               & gradBoundaryModesMList[0](row, column));
    //     }
    // }

    // // M21
    // RectangularMatrix<scalar> localOIFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < localOIFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < localOIFaceMatrix.n(); ++column)
    //     {
    //         localOIFaceMatrix(row, column) =  1/2 * (gradBoundaryModesMList[0](column, row) & faceNormal) 
    //                                               * boundaryModesMList[1](row, column)
    //                                          +1/2 * epsilonPara * boundaryModesMList[0](column, row) 
    //                                               * (gradBoundaryModesMList[1](row, column) & faceNormal)
    //                                          -xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[0](column, row) 
    //                                               * boundaryModesMList[1](row, column)
    //                                          -xigema1/Foam::pow(elementVol, beta1) * (gradBoundaryModesMList[0](column, row) 
    //                                               & gradBoundaryModesMList[1](row, column));
    //     }
    // }

    // // boundary penalty terms
    // // Min
    // RectangularMatrix<scalar> bounInFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < bounInFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < bounInFaceMatrix.n(); ++column)
    //     {
    //         bounInFaceMatrix(row, column) = -(gradBoundaryModesMList[0](column, row) & faceNormal) 
    //                                             * boundaryModesMList[0](row, column)
    //                                         +epsilonPara * boundaryModesMList[0](column, row) 
    //                                             * (gradBoundaryModesMList[0](row, column) & faceNormal)
    //                                         +xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[0](column, row) 
    //                                             * boundaryModesMList[0](row, column);
    //     }
    // }

    // // Mout
    // RectangularMatrix<scalar> bounOutFaceMatrix(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < bounOutFaceMatrix.m(); ++row)
    // {
    //     for (label column = 0; column < bounOutFaceMatrix.n(); ++column)
    //     {
    //         bounOutFaceMatrix(row, column) = -(gradBoundaryModesMList[1](column, row) & faceNormal) 
    //                                             * boundaryModesMList[1](row, column)
    //                                         +epsilonPara * boundaryModesMList[1](column, row) 
    //                                             * (gradBoundaryModesMList[1](row, column) & faceNormal)
    //                                         +xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[1](column, row) 
    //                                             * boundaryModesMList[1](row, column);
    //     }
    // }

    // // Fin
    // scalar Tin(273.0);
    // RectangularMatrix<scalar> bounInPenaltyMatrix(modesNum, 1, Foam::Zero);
    // for (label row = 0; row < bounInPenaltyMatrix.m(); ++row)
    // {
        
        
        
    //     for(label column = 0; column < gradBoundaryModesMList.n(); ++column)
    //     bounInPenaltyMatrix(row, 0) =  (gradBoundaryModesMList[0](0, row) & faceNormal) * Tin
    //                                   + xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[0](0, row);
    // }

    // // Fout
    // scalar Tout(573.0);
    // RectangularMatrix<scalar> bounOutPenaltyMatrix(modesNum, 1, Foam::Zero);
    // for (label row = 0; row < bounOutPenaltyMatrix.m(); ++row)
    // {
    //     // for(label column = 0; column < bounOutFaceMatrix.n(); ++column)
    //     bounOutPenaltyMatrix(row, 0) =   (gradBoundaryModesMList[1](0, row) & faceNormal) * Tout
    //                                   + xigema0/Foam::pow(elementVol, beta0) * boundaryModesMList[1](0, row);
    // }

    // // global matrix
    // // global phi matrix
    // for (label row = 0; row < modesNum; ++row)
    // {
    //     for (label column = 0; column < modesNum; ++column)
    //     {
    //         globalphiMatrix(row, column) = localphiMatrix(row, column) + localInFaceMatrix(row, column) 
    //                                         + bounInFaceMatrix(row, column);
    //     }
    // }

    // for(label elementI = 1; elementI < 4 ; ++elementI)
    // {
    //     for (label row = 0; row < modesNum; ++row)
    //     {
    //         for (label column = 0; column < modesNum; ++column)
    //         {
    //             globalphiMatrix(row+elementI*modesNum, column+elementI*modesNum) =  localphiMatrix(row, column) 
    //                                                                     + localOutFaceMatrix(row, column) 
    //                                                                     + bounInFaceMatrix(row, column);
    //         }
    //     }
    // }

    // for (label row = 0; row < modesNum; ++row)
    // {
    //     for (label column = 0; column < modesNum; ++column)
    //     {
    //         globalphiMatrix((modesNum-1)*modesNum, (modesNum-1)*modesNum) =  localphiMatrix(row, column) 
    //                                                                        + localOutFaceMatrix(row, column) 
    //                                                                        + bounOutFaceMatrix(row, column);
    //     }
    // }


    // global F matrix

    
    

    // solve matrix system


 


    

    return 0;
}


// ************************************************************************* //

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
    Info << "matchPatchID: " << matchPatchID << endl;

    // ptrlist to field value
    PtrList<volScalarField> fieldModesList;
    PtrList<FieldField<Foam::fvPatchField, scalar>> fieldBundaryModesList;
    
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
        fieldBundaryModesList.append(fieldValueMode.boundaryField().clone());
    }



    PtrList<volVectorField> gradfieldModesList;
    PtrList<FieldField<Foam::fvPatchField, vector>> gradfieldBundaryModesList;

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
        gradfieldBundaryModesList.append(fieldValueModegrad.boundaryField().clone());

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

    // // The following two statement are equal
    // // ---1 is matrix multiplation
    // List<scalar> testList(gradModesM.m());
    // for (label row = 0; row < gradModesM.m(); ++row)
    // {
    //     testList[row] = gradModesM(row, 0) & gradModesM(row, 1);
    // }
    // Info << "matrix: " << gSum(testList) << endl;
    // // ---2 is field value multiplation
    // Info << "fvm: " << gSum(volScalarField (gradfieldModesList[0] & gradfieldModesList[1])) << endl;

    // // The following two statement are equal
    // // ---1 is matrix multiplation
    // List<scalar> testList(gradBoundaryModesMList[0].m());
    // for (label row = 0; row < gradBoundaryModesMList[0].m(); ++row)
    // {
    //     testList[row] = gradBoundaryModesMList[0](row, 0) & gradBoundaryModesMList[1](row, 0);
    // }
    // Info << "matrix: " << gSum(testList) << endl;
    // // ---2 is field value multiplation
    // Info << "fvm: " << gSum(scalarField (gradfieldModesList[0].boundaryField()[0] & gradfieldModesList[0].boundaryField()[1])) << endl;

    // create Matrix system
    // initial global matrix
    RectangularMatrix<scalar> globalphiMatrix(modesNum * modesNum, modesNum * modesNum, Foam::Zero);
    RectangularMatrix<scalar> globalFMmatrix(modesNum * modesNum, 1, Foam::Zero);

    List<scalar> tempList;
    tempList.resize(mesh.C().size());

    // volumtric contribution
    RectangularMatrix<scalar> localphiMatrix(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < localphiMatrix.m(); ++row)
    {
        for (label column = 0; column < localphiMatrix.n(); ++column)
        {
            localphiMatrix(row, column) = gSum(scalarField (gradfieldModesList[row] & gradfieldModesList[column]));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"localphiMatrix";
    writeMatrix(localphiMatrix, dataFile);

    
    // Matrix IO
    // ---1, matrix write 
    // dataFile = mesh.time().path()/"SVD"/"localphiMatrixTest";
    // OFstream writeOF(dataFile);
    // localphiMatrix.writeMatrix(writeOF);
    // ---2, matrix read 
    // dataFile = mesh.time().path()/"SVD"/"localphiMatrixTest";
    // IFstream readOF(dataFile);
    // RectangularMatrix<scalar> localphiMatrix2(modesNum, modesNum, Foam::Zero);
    // localphiMatrix2.readMatrix(readOF);
    // Info << localphiMatrix2 << endl;


    // interface contribution
    vector faceNormal(0, 0, 1);

    // M11
    RectangularMatrix<scalar> M11(modesNum, modesNum, Foam::Zero);
    label bundaryPatch1(matchPatchID[0]);

    for (label row = 0; row < M11.m(); ++row)
    {
        for (label column = 0; column < M11.n(); ++column)
        {
            M11(row, column) = gSum(scalarField (
                                            - 1/2 * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & faceNormal) 
                                            + 1/2 * epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch1] & faceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema1/Foam::pow(elementVol, beta1) * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch1])));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M11";
    writeMatrix(M11, dataFile);

    // M22
    RectangularMatrix<scalar> M22(modesNum, modesNum, Foam::Zero);
    label bundaryPatch2(matchPatchID[1]);

    for (label row = 0; row < M22.m(); ++row)
    {
        for (label column = 0; column < M22.n(); ++column)
        {
            M22(row, column) = gSum(scalarField (
                                              1/2 * fieldBundaryModesList[row][bundaryPatch2] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch2] & faceNormal) 
                                            - 1/2 * epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch2] & faceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch2]
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            + xigema1/Foam::pow(elementVol, beta1) * (gradfieldBundaryModesList[row][bundaryPatch2]
                                                & gradfieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M22";
    writeMatrix(M22, dataFile);

    // M12
    RectangularMatrix<scalar> M12(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M12.m(); ++row)
    {
        for (label column = 0; column < M12.n(); ++column)
        {
            M12(row, column) = gSum(scalarField (
                                            - 1/2 * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch2] & faceNormal) 
                                            - 1/2 * epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch1] & faceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema1/Foam::pow(elementVol, beta1) * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M12";
    writeMatrix(M12, dataFile);

    // M21
    RectangularMatrix<scalar> M21(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M21.m(); ++row)
    {
        for (label column = 0; column < M21.n(); ++column)
        {
            M21(row, column) = gSum(scalarField (
                                              1/2 * fieldBundaryModesList[row][bundaryPatch2] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & faceNormal) 
                                            + 1/2 * epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch2] & faceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch2]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema1/Foam::pow(elementVol, beta1) * (gradfieldBundaryModesList[row][bundaryPatch2]
                                                & gradfieldBundaryModesList[column][bundaryPatch1])));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M21";
    writeMatrix(M21, dataFile);

    // boundary penalty terms
    // Min
    RectangularMatrix<scalar> Min(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < Min.m(); ++row)
    {
        for (label column = 0; column < Min.n(); ++column)
        {
            Min(row, column) = gSum(scalarField (
                                - fieldBundaryModesList[row][bundaryPatch1] 
                                    * (gradfieldBundaryModesList[column][bundaryPatch1] & faceNormal) 
                                + epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch1] & faceNormal) 
                                    * fieldBundaryModesList[column][bundaryPatch1]
                                + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch1]
                                    * fieldBundaryModesList[column][bundaryPatch1]));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"Min";
    writeMatrix(Min, dataFile);

    // Mout
    RectangularMatrix<scalar> Mout(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < Mout.m(); ++row)
    {
        for (label column = 0; column < Mout.n(); ++column)
        {
            Mout(row, column) = gSum(scalarField (
                                - fieldBundaryModesList[row][bundaryPatch2] 
                                    * (gradfieldBundaryModesList[column][bundaryPatch2] & faceNormal) 
                                + epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch2] & faceNormal) 
                                    * fieldBundaryModesList[column][bundaryPatch2]
                                + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch2]
                                    * fieldBundaryModesList[column][bundaryPatch2]));
        }
    }
    dataFile = mesh.time().path()/"SVD"/"Mout";
    writeMatrix(Mout, dataFile);

    // Fin
    scalar Tin(273.0);
    RectangularMatrix<scalar> Fin(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fin.m(); ++row)
    {
        Fin(row, 0) = gSum(scalarField (
                            epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch1] & faceNormal) 
                            * Tin
                            + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch1]
                            * Tin));
    }

    // Fout
    scalar Tout(573.0);
    RectangularMatrix<scalar> Fout(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fout.m(); ++row)
    {
        Fout(row, 0) = gSum(scalarField (
                            epsilonPara * (gradfieldBundaryModesList[row][bundaryPatch2] & faceNormal) 
                            * Tout
                            + xigema0/Foam::pow(elementVol, beta0) * fieldBundaryModesList[row][bundaryPatch2]
                            * Tout));
    }


    // global matrix
    // global phi matrix
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            globalphiMatrix(row, column) = localphiMatrix(row, column) + Min(row, column) + M11(row, column);
        }
    }

    for(label elementI = 1; elementI < 4 ; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalphiMatrix(row+elementI*modesNum, column+elementI*modesNum) =  localphiMatrix(row, column) 
                                                                        + M11(row, column) + M22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            globalphiMatrix(row+(modesNum-1)*modesNum, column+(modesNum-1)*modesNum) =  localphiMatrix(row, column) 
                                                                           + M22(row, column) + Mout(row, column);
        }
    }

    for(label elementI = 0; elementI < 4 ; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalphiMatrix(row+elementI*modesNum, column+(elementI+1)*modesNum) =  M12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < 4 ; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalphiMatrix(row+(elementI+1)*modesNum, column+elementI*modesNum) =  M21(row, column);
            }
        }
    }
    dataFile = mesh.time().path()/"SVD"/"globalphiMatrix";
    writeMatrix(globalphiMatrix, dataFile);


    // global F matrix
    for (label row = 0; row < modesNum; ++row)
    {
        globalFMmatrix(row, 0) = Fin(row, 0);
    }

    for (label row = 0; row < modesNum; ++row)
    {
        globalFMmatrix(row+(modesNum-1)*modesNum, 0) = Fout(row, 0);
    }
    

    // // solve matrix system
    // RectangularMatrix<scalar> calCoefficientM;
    // calCoefficientM = SVDinv(globalphiMatrix) * globalFMmatrix;
    // dataFile = mesh.time().path()/"SVD"/"calCoefficientM";
    // writeMatrix(calCoefficientM, dataFile);


    return 0;
}


// ************************************************************************* //

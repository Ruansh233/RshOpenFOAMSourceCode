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
    fileName dataFile;

    #include "readDGdict.H"

    // element volume
    scalar elementVol(gSum(mesh.V()));


    // ===========================================================
    // ------ read modesM and boundaryModesM from Block1  --------
    // ------ read grad and div of modes -------------------------
    // ===========================================================

    // create field value
    List<fileName> modeNames (modesNum);
    for (label iname = 0; iname < modesNum; ++iname)
    {
        modeNames[iname] = "Tmode" + name(iname + 1);
    }

    // set patch value
    List<word> boundaryWordRe ({".*_in", ".*_out", ".*_wall2"});
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
    Info << "matchPatchID: " << matchPatchID << endl << endl;

    // ptrlist to field value and grad field value
    PtrList<volScalarField> fieldModesList;
    PtrList<FieldField<Foam::fvPatchField, scalar>> fieldBundaryModesList;
    PtrList<volVectorField> gradfieldModesList;
    PtrList<FieldField<Foam::fvPatchField, vector>> gradfieldBundaryModesList;
    
    forAll(modeNames, No_)
    {       
        volScalarField fieldValueMode
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

        fieldModesList.append(fieldValueMode.clone());
        fieldBundaryModesList.append(fieldValueMode.boundaryField().clone());

        volVectorField fieldValueModegrad
        (
            IOobject
            (
                "grad" + modeNames[No_],
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );

        gradfieldModesList.append(fieldValueModegrad.clone());
        gradfieldBundaryModesList.append(fieldValueModegrad.boundaryField().clone());
    }

    
    // ===========================================================
    // -------- creating and solving the matrix system -----------
    // ===========================================================
    // create Matrix system
    // initial global matrix
    RectangularMatrix<scalar> globalphiMatrix(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    RectangularMatrix<scalar> globalFMmatrix(modesNum * elementNum, 1, Foam::Zero);

    List<scalar> tempList;
    tempList.resize(mesh.C().size());

    // Info<< "test: " << endl
    //     << gSum(scalarField(gradfieldModesList[0] & gradfieldModesList[1] * mesh.V())) << endl;

    // volumtric contribution
    RectangularMatrix<scalar> localphiMatrix(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < localphiMatrix.m(); ++row)
    {
        for (label column = 0; column < localphiMatrix.n(); ++column)
        {
            localphiMatrix(row, column) = heatConductivity * gSum(scalarField (gradfieldModesList[row] & gradfieldModesList[column]
                                                                               * mesh.V()));
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
    vector interfaceNormal(0, 0, 1);

    // M11
    RectangularMatrix<scalar> M11(modesNum, modesNum, Foam::Zero);
    label bundaryPatch1(matchPatchID[0]);  // in
    label bundaryPatch2(matchPatchID[1]);  // out
    label bundaryPatch3(matchPatchID[2]);  // wall

    // Info<< "test: " << endl
    //     << gSum(scalarField(gradfieldBundaryModesList[0][bundaryPatch2]
    //                         & gradfieldBundaryModesList[1][bundaryPatch2]) * mesh.boundary()[bundaryPatch2].magSf()) << endl;


    for (label row = 0; row < M11.m(); ++row)
    {
        for (label column = 0; column < M11.n(); ++column)
        {
            M11(row, column) = gSum(scalarField (
                                            - 0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch2] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch2] & interfaceNormal) 
                                            + 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            + xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            + xigema1 * (gradfieldBundaryModesList[row][bundaryPatch2]
                                                & gradfieldBundaryModesList[column][bundaryPatch2]))
                                                * mesh.boundary()[bundaryPatch2].magSf());
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M11";
    writeMatrix(M11, dataFile);

    // M22
    RectangularMatrix<scalar> M22(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M22.m(); ++row)
    {
        for (label column = 0; column < M22.n(); ++column)
        {
            M22(row, column) = gSum(scalarField (
                                              0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & interfaceNormal) 
                                            - 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema1 * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch1]))
                                                * mesh.boundary()[bundaryPatch1].magSf());
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
                                            - 0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch2] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & interfaceNormal) 
                                            - 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema1 * (gradfieldBundaryModesList[row][bundaryPatch2]
                                                & gradfieldBundaryModesList[column][bundaryPatch1]))
                                                * mesh.boundary()[bundaryPatch1].magSf());
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
                                              0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch2] & interfaceNormal) 
                                            + 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema1 * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch2]))
                                                * mesh.boundary()[bundaryPatch2].magSf());
        }
    }
    dataFile = mesh.time().path()/"SVD"/"M21";
    writeMatrix(M21, dataFile);

    // boundary penalty terms
    // Min
    RectangularMatrix<scalar> Min(modesNum, modesNum, Foam::Zero);
    vector inletfaceNormal(0, 0, -1);

    for (label row = 0; row < Min.m(); ++row)
    {
        for (label column = 0; column < Min.n(); ++column)
        {
            Min(row, column) = gSum(scalarField (
                                - heatConductivity * fieldBundaryModesList[row][bundaryPatch1] 
                                    * (gradfieldBundaryModesList[column][bundaryPatch1] & inletfaceNormal) 
                                + epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & inletfaceNormal) 
                                    * fieldBundaryModesList[column][bundaryPatch1]
                                + xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                                    * fieldBundaryModesList[column][bundaryPatch1])
                                    * mesh.boundary()[bundaryPatch1].magSf());
        }
    }
    dataFile = mesh.time().path()/"SVD"/"Min";
    writeMatrix(Min, dataFile);

    // Mout
    RectangularMatrix<scalar> Mout(modesNum, modesNum, Foam::Zero);
    vector outletfaceNormal(0, 0, 1);

    for (label row = 0; row < Mout.m(); ++row)
    {
        for (label column = 0; column < Mout.n(); ++column)
        {
            Mout(row, column) = gSum(scalarField (
                                - heatConductivity * fieldBundaryModesList[row][bundaryPatch2] 
                                    * (gradfieldBundaryModesList[column][bundaryPatch2] & outletfaceNormal) 
                                + epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & outletfaceNormal) 
                                    * fieldBundaryModesList[column][bundaryPatch2]
                                + xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                                    * fieldBundaryModesList[column][bundaryPatch2])
                                    * mesh.boundary()[bundaryPatch2].magSf());
        }
    }
    dataFile = mesh.time().path()/"SVD"/"Mout";
    writeMatrix(Mout, dataFile);

    // Fin
    RectangularMatrix<scalar> Fin(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fin.m(); ++row)
    {
        Fin(row, 0) = gSum(scalarField (
                            epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & inletfaceNormal) 
                            * Tin
                            + xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                            * Tin)
                            * mesh.boundary()[bundaryPatch1].magSf());
    }

    // Fout
    RectangularMatrix<scalar> Fout(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fout.m(); ++row)
    {
        Fout(row, 0) = gSum(scalarField (
                            epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & outletfaceNormal) 
                            * Tout
                            + xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                            * Tout)
                            * mesh.boundary()[bundaryPatch2].magSf());
    }

    // Fn
    RectangularMatrix<scalar> Fn(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fn.m(); ++row)
    {
        Fn(row, 0) = gSum(scalarField (
                            heatConductivity * fieldBundaryModesList[row][bundaryPatch3] * qn)
                            * mesh.boundary()[bundaryPatch3].magSf());
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

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
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
            globalphiMatrix(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  localphiMatrix(row, column) 
                                                                           + M22(row, column) + Mout(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalphiMatrix(row+elementI*modesNum, column+(elementI+1)*modesNum) =  M12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
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
        // globalFMmatrix(row, 0) = Fin(row, 0) + Fn(row, 0);
        globalFMmatrix(row, 0) = Fin(row, 0);
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            globalFMmatrix(row+elementI*modesNum, 0) = Fn(row, 0);
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        // globalFMmatrix(row+(elementNum-1)*modesNum, 0) = Fout(row, 0) + Fn(row, 0);
        globalFMmatrix(row+(elementNum-1)*modesNum, 0) = Fout(row, 0);
    }
    dataFile = mesh.time().path()/"SVD"/"globalFMmatrix";
    writeMatrix(globalFMmatrix, dataFile);


    // solve matrix system
    RectangularMatrix<scalar> tempCalCoefficientM;
    RectangularMatrix<scalar> calCoefficientM(modesNum, elementNum);
    tempCalCoefficientM = SVDinv(globalphiMatrix) * globalFMmatrix;
    for (label row = 0; row < calCoefficientM.m(); ++row)
    {
        for (label column = 0; column < calCoefficientM.n(); ++column)
        {
            calCoefficientM(row, column) =  tempCalCoefficientM(row+column*modesNum, 0);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"calCoefficientM";
    writeMatrix(calCoefficientM, dataFile);

    // read the matrix
    // read cell modes matrix
    dataFile = dataPath/"modesM";
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

    // calculate snapshots
    RectangularMatrix<scalar> calSnapshotsM;
    calSnapshotsM = modesM * calCoefficientM;
    dataFile = mesh.time().path()/"SVD"/"calSnapshotsM";
    writeMatrix(calSnapshotsM, dataFile);


    // The error matrix
    dataFile = dataPath/"snapshotsM";
    RectangularMatrix<scalar> snapshotsM(mesh.C().size(), elementNum);
    if(isFile(dataFile))
    {                
        IFstream dataStream(dataFile);
        word dataLine;
        label row(0);

        while(dataStream.getLine(dataLine) && dataLine != word::null)
        {
            IStringStream dataString (dataLine);
            token singleData;  // token stores the data read from IFstream 

            for(label elementI = 0; elementI < elementNum; ++elementI)
            {
                dataString.read(singleData);    
                snapshotsM(row, elementI) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }

    RectangularMatrix<scalar> errorM(mesh.C().size(), elementNum);
    for (label row = 0; row < errorM.m(); ++row)
    {
        for (label column = 0; column < errorM.n(); ++column)
        {
            errorM(row, column) =  (calSnapshotsM(row, column) - snapshotsM(row, column))/snapshotsM(row, column);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"errorM";
    writeMatrix(errorM, dataFile);    


    return 0;
}


// ************************************************************************* //

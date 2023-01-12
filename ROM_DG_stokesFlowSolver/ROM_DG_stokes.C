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

    fileName dataPath (runTime.globalPath()/"SVD");

    #include "readDGdict.H"

    // ===========================================================
    // ------ PtrList for modes, gradModes and ... ---------------
    // ===========================================================
    PtrList<volScalarField> pFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, scalar>> pFieldBundaryModesList;

    PtrList<volVectorField> uFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, vector>> uFieldBundaryModesList;
    PtrList<volTensorField> graduFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, tensor>> graduFieldBundaryModesList;
    PtrList<volScalarField> divuFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, scalar>> divuFieldBundaryModesList;

    // ===========================================================
    // ------ reading modes, gradModes from the ref case ---------
    // ------ creating mesh and time object for ref case ---------
    // ===========================================================
    // reference case name
    fileName refCaseName(DGdict.getWord("refCaseName"));
    // time object for reference cases
    Foam::Time runTimeRef
    (
        Foam::Time::controlDictName,
        args.rootPath(),
        refCaseName,
        "system",
        "constant"
    );

    // create new mesh object for reference cases
    fvMesh  refElementMesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            args.rootPath()/refCaseName/"constant",
            runTimeRef,
            IOobject::MUST_READ
        ),
        false
    );

    // modes name list
    List<fileName> pModeNames (modesNum);
    List<fileName> uModeNames (modesNum);
    forAll(pModeNames, nameI)
    {
        pModeNames[nameI] = "pMode" + name(nameI + 1);
        uModeNames[nameI] = "uMode" + name(nameI + 1);
    }

    // read pressure modes
    forAll(pModeNames, No_)
    {
        // pressure modes
        volScalarField pFieldValueMode
        (
            IOobject
            (
                pModeNames[No_],
                runTimeRef.timeName(),
                refElementMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh
        );
        pFieldModesList.append(pFieldValueMode.clone());
        pFieldBundaryModesList.append(pFieldValueMode.boundaryField().clone());
    }

    // read velocity modes
    forAll(uModeNames, No_)
    {
        // velocity modes
        volVectorField uFieldValueMode
        (
            IOobject
            (
                uModeNames[No_],
                runTimeRef.timeName(),
                refElementMesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            refElementMesh
        );
        uFieldModesList.append(uFieldValueMode.clone());
        uFieldBundaryModesList.append(uFieldValueMode.boundaryField().clone());

        // grad of velocity modes
        volTensorField uFieldValueModegrad
        (
            IOobject
            (
                "grad" + uModeNames[No_],
                runTimeRef.timeName(),
                refElementMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(uFieldValueMode)
        );

        graduFieldModesList.append(uFieldValueModegrad.clone());
        graduFieldBundaryModesList.append(uFieldValueModegrad.boundaryField().clone());

        // div of velocity modes
        volScalarField uFieldValueModediv
        (
            IOobject
            (
                "grad" + uModeNames[No_],
                runTimeRef.timeName(),
                refElementMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::div(uFieldValueMode)
        );
        divuFieldModesList.append(uFieldValueModediv.clone());
        divuFieldBundaryModesList.append(uFieldValueModediv.boundaryField().clone());
    }

    // ===========================================================
    // ------ The matrix system for momentum equation ------------
    // ===========================================================
    // create Matrix system
    // initial global matrix of momentum equations
    RectangularMatrix<scalar> MomGlobalAMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    RectangularMatrix<scalar> MomGlobalBMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    RectangularMatrix<scalar> MomGlobalFMat(modesNum * elementNum, 1, Foam::Zero);
    fileName dataFile;

    // ===========================================================
    // ------ The diffusion term ---------------------------------
    // ===========================================================
    // volumtric contribution, diffusion term
    RectangularMatrix<scalar> MomLocalAMat(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < MomLocalAMat.m(); ++row)
    {
        for (label column = 0; column < MomLocalAMat.n(); ++column)
        {
            MomLocalAMat(row, column) = gSum(scalarField (graduFieldModesList[row] && graduFieldModesList[column]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MomLocalAMat";
    writeMatrix(MomLocalAMat, dataFile);


    // patch ID for reference mesh of different interface
    label bundaryPatch1 (refElementMesh.boundary().findPatchID("block1_out"));
    label bundaryPatch2 (refElementMesh.boundary().findPatchID("block1_in"));


    // interface contribution
    vector interfaceNormal(0, 0, 1);


    // M11
    RectangularMatrix<scalar> M11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < M11.m(); ++row)
    {
        for (label column = 0; column < M11.n(); ++column)
        {
            M11(row, column) = gSum(scalarField (
                                            - 0.5 * (uFieldBundaryModesList[row][bundaryPatch1] 
                                                & graduFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            + 0.5 * epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal
                                                & uFieldBundaryModesList[column][bundaryPatch1])
                                            + xigema0 * (uFieldBundaryModesList[row][bundaryPatch1]
                                                & uFieldBundaryModesList[column][bundaryPatch1])
                                            + xigema1 * (graduFieldBundaryModesList[row][bundaryPatch1]
                                                && graduFieldBundaryModesList[column][bundaryPatch1])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M11";
    writeMatrix(M11, dataFile);

    // M22
    RectangularMatrix<scalar> M22(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M22.m(); ++row)
    {
        for (label column = 0; column < M22.n(); ++column)
        {
            M22(row, column) = 

                                gSum(scalarField (
                                            + 0.5 * (uFieldBundaryModesList[row][bundaryPatch2] 
                                                & graduFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal
                                                & uFieldBundaryModesList[column][bundaryPatch2])
                                            + xigema0 * (uFieldBundaryModesList[row][bundaryPatch2]
                                                & uFieldBundaryModesList[column][bundaryPatch2])
                                            + xigema1 * (graduFieldBundaryModesList[row][bundaryPatch2]
                                                && graduFieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M22";
    writeMatrix(M22, dataFile);

    // M12
    RectangularMatrix<scalar> M12(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M12.m(); ++row)
    {
        for (label column = 0; column < M12.n(); ++column)
        {
            M12(row, column) = gSum(scalarField (
                                            - 0.5 * (uFieldBundaryModesList[row][bundaryPatch1] 
                                                & graduFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal
                                                & uFieldBundaryModesList[column][bundaryPatch2])
                                            - xigema0 * (uFieldBundaryModesList[row][bundaryPatch1]
                                                & uFieldBundaryModesList[column][bundaryPatch2])
                                            - xigema1 * (graduFieldBundaryModesList[row][bundaryPatch1]
                                                && graduFieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M12";
    writeMatrix(M12, dataFile);

    // M21
    RectangularMatrix<scalar> M21(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < M21.m(); ++row)
    {
        for (label column = 0; column < M21.n(); ++column)
        {
            M21(row, column) = gSum(scalarField (
                                            + 0.5 * (uFieldBundaryModesList[row][bundaryPatch2] 
                                                & graduFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            + 0.5 * epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal
                                                & uFieldBundaryModesList[column][bundaryPatch1])
                                            - xigema0 * (uFieldBundaryModesList[row][bundaryPatch2]
                                                & uFieldBundaryModesList[column][bundaryPatch1])
                                            - xigema1 * (graduFieldBundaryModesList[row][bundaryPatch2]
                                                && graduFieldBundaryModesList[column][bundaryPatch1])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M21";
    writeMatrix(M21, dataFile);

    // boundary penalty terms
    // MtaoD, patch-inlet, bundaryPatch2
    RectangularMatrix<scalar> MtaoD(modesNum, modesNum, Foam::Zero);
    vector inletfaceNormal(0, 0, -1);

    for (label row = 0; row < MtaoD.m(); ++row)
    {
        for (label column = 0; column < MtaoD.n(); ++column)
        {
            MtaoD(row, column) = gSum(scalarField (
                                - (uFieldBundaryModesList[row][bundaryPatch2] 
                                    & graduFieldBundaryModesList[column][bundaryPatch2] & inletfaceNormal)
                                + epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch2] & inletfaceNormal
                                    & uFieldBundaryModesList[column][bundaryPatch2])
                                + xigema0 * (uFieldBundaryModesList[row][bundaryPatch2]
                                    & uFieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MtaoD";
    writeMatrix(MtaoD, dataFile);

    // FtaoD, patch-inlet, bundaryPatch2
    RectangularMatrix<scalar> FtaoD(modesNum, 1, Foam::Zero);
    for (label row = 0; row < FtaoD.m(); ++row)
    {
        FtaoD(row, 0) = gSum(scalarField (
                            epsilonPara * (graduFieldBundaryModesList[row][bundaryPatch2] & inletfaceNormal
                            & Uin)
                            + xigema0 * (uFieldBundaryModesList[row][bundaryPatch2]
                            & Uin)));
    }
    dataFile = runTime.globalPath()/"SVD"/"FtaoD";
    writeMatrix(FtaoD, dataFile);

    // ===========================================================
    // ------ The pressure gradient term -------------------------
    // ===========================================================
    // volumtric contribution, \nabla p
    RectangularMatrix<scalar> MomLocalBMat(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < MomLocalBMat.m(); ++row)
    {
        for (label column = 0; column < MomLocalBMat.n(); ++column)
        {
            MomLocalBMat(row, column) = gSum(scalarField (pFieldModesList[column] * divuFieldModesList[row]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MomLocalBMat";
    writeMatrix(MomLocalBMat, dataFile);

    // N11
    RectangularMatrix<scalar> N11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < N11.m(); ++row)
    {
        for (label column = 0; column < N11.n(); ++column)
        {
            // N11(row, column) = gSum(scalarField (
            //                                 + 0.5 * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)
            //                                 + 0.5 * epsilonPara * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)));
            
            N11(row, column) = gSum(scalarField ( + 0.5 * pFieldBundaryModesList[column][bundaryPatch1] 
                                                * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)));                                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"N11";
    writeMatrix(N11, dataFile);

    // N22
    RectangularMatrix<scalar> N22(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < N22.m(); ++row)
    {
        for (label column = 0; column < N22.n(); ++column)
        {
            // N22(row, column) = gSum(scalarField (
            //                                 - 0.5 * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)));

            N22(row, column) = gSum(scalarField ( - 0.5 * pFieldBundaryModesList[column][bundaryPatch2] 
                                                * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"N22";
    writeMatrix(N22, dataFile);

    // N12
    RectangularMatrix<scalar> N12(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < N12.m(); ++row)
    {
        for (label column = 0; column < N12.n(); ++column)
        {
            // N12(row, column) = gSum(scalarField (
            //                                 + 0.5 * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[column][bundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)));

            N12(row, column) = gSum(scalarField ( + 0.5 * pFieldBundaryModesList[column][bundaryPatch2] 
                                                * (uFieldBundaryModesList[row][bundaryPatch1] & interfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"N12";
    writeMatrix(N12, dataFile);

    // N21
    RectangularMatrix<scalar> N21(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < N21.m(); ++row)
    {
        for (label column = 0; column < N21.n(); ++column)
        {
            // N21(row, column) = gSum(scalarField (
            //                                 - 0.5 * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[column][bundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)));

            N21(row, column) = gSum(scalarField ( - 0.5 * pFieldBundaryModesList[column][bundaryPatch1] 
                                                * (uFieldBundaryModesList[row][bundaryPatch2] & interfaceNormal)));                                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"N21";
    writeMatrix(N21, dataFile);

    // boundary penalty terms
    // NtaoD, patch-inlet, bundaryPatch2
    RectangularMatrix<scalar> NtaoD(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < NtaoD.m(); ++row)
    {
        for (label column = 0; column < NtaoD.n(); ++column)
        {
            NtaoD(row, column) = gSum(scalarField (pFieldBundaryModesList[column][bundaryPatch2] 
                                                * (uFieldBundaryModesList[row][bundaryPatch2] & inletfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"NtaoD";
    writeMatrix(NtaoD, dataFile);

    // ===========================================================
    // ------ Assign value for global momentum matrix ------------
    // ===========================================================
    // create Matrix system
    // initial global matrix of momentum equations
    // RectangularMatrix<scalar> MomGlobalAMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    // RectangularMatrix<scalar> MomGlobalBMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    // RectangularMatrix<scalar> MomGlobalFMat(modesNum * elementNum, 1, Foam::Zero);

    // ===========================================================
    // --------------- MomGlobalAMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            MomGlobalAMat(row, column) = MomLocalAMat(row, column) + MtaoD(row, column) + M11(row, column);
        }
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalAMat(row+elementI*modesNum, column+elementI*modesNum) =  MomLocalAMat(row, column) 
                                                                        + M11(row, column) + M22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            MomGlobalAMat(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  MomLocalAMat(row, column) 
                                                                           + M22(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalAMat(row+elementI*modesNum, column+(elementI+1)*modesNum) =  M12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalAMat(row+(elementI+1)*modesNum, column+elementI*modesNum) =  M21(row, column);
            }
        }
    }
    dataFile = mesh.time().path()/"SVD"/"MomGlobalAMat";
    writeMatrix(MomGlobalAMat, dataFile);

    // ===========================================================
    // --------------- MomGlobalBMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            MomGlobalBMat(row, column) = MomLocalBMat(row, column) + N11(row, column);
        }
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalBMat(row+elementI*modesNum, column+elementI*modesNum) =  MomLocalBMat(row, column) 
                                                                        + N11(row, column) + N22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            MomGlobalBMat(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  MomLocalBMat(row, column) 
                                                                           + N22(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalBMat(row+elementI*modesNum, column+(elementI+1)*modesNum) =  N12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomGlobalBMat(row+(elementI+1)*modesNum, column+elementI*modesNum) =  N21(row, column);
            }
        }
    }
    dataFile = mesh.time().path()/"SVD"/"MomGlobalBMat";
    writeMatrix(MomGlobalBMat, dataFile);

    // ===========================================================
    // --------------- MomGlobalFMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        MomGlobalFMat(row, 0) = FtaoD(row, 0);
    }

    dataFile = mesh.time().path()/"SVD"/"MomGlobalFMat";
    writeMatrix(MomGlobalFMat, dataFile);
    

    // ===========================================================
    // ------ The matrix system for continuous equation ------------
    // ===========================================================
    // create Matrix system
    // initial global matrix of continuous equations
    RectangularMatrix<scalar> ConGlobalBMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    RectangularMatrix<scalar> ConGlobalFMat(modesNum * elementNum, 1, Foam::Zero);

    // ===========================================================
    // ------ The divergence of velocity -------------------------
    // ===========================================================
    // volumtric contribution, \nabla \cdot v
    RectangularMatrix<scalar> ConLocalBMat(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < ConLocalBMat.m(); ++row)
    {
        for (label column = 0; column < ConLocalBMat.n(); ++column)
        {
            ConLocalBMat(row, column) = gSum(scalarField (pFieldModesList[row] * divuFieldModesList[column]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"ConLocalBMat";
    writeMatrix(ConLocalBMat, dataFile);

    // K11
    RectangularMatrix<scalar> K11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < K11.m(); ++row)
    {
        for (label column = 0; column < K11.n(); ++column)
        {
            K11(row, column) = gSum(scalarField (
                                            + 0.5 * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            + 0.5 * epsilonPara * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            + xigema0 * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)));                                          
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"K11";
    writeMatrix(K11, dataFile);

    // K22
    RectangularMatrix<scalar> K22(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < K22.m(); ++row)
    {
        for (label column = 0; column < K22.n(); ++column)
        {
            K22(row, column) = gSum(scalarField (
                                            - 0.5 * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            + xigema0 * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"K22";
    writeMatrix(K22, dataFile);

    // K12
    RectangularMatrix<scalar> K12(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < K12.m(); ++row)
    {
        for (label column = 0; column < K12.n(); ++column)
        {
            K12(row, column) = gSum(scalarField (
                                            + 0.5 * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            - 0.5 * epsilonPara * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)
                                            - xigema0 * pFieldBundaryModesList[row][bundaryPatch2] 
                                                * (uFieldBundaryModesList[column][bundaryPatch1] & interfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"K12";
    writeMatrix(K12, dataFile);

    // K21
    RectangularMatrix<scalar> K21(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < K21.m(); ++row)
    {
        for (label column = 0; column < K21.n(); ++column)
        {
            K21(row, column) = gSum(scalarField (
                                            - 0.5 * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)
                                            - xigema0 * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & interfaceNormal)));                                  
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"K21";
    writeMatrix(K21, dataFile);

    // boundary penalty terms
    // KtaoD, patch-inlet, bundaryPatch2
    RectangularMatrix<scalar> KtaoD(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < KtaoD.m(); ++row)
    {
        for (label column = 0; column < KtaoD.n(); ++column)
        {
            KtaoD(row, column) = gSum(scalarField (epsilonPara * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & inletfaceNormal)
                                            + xigema0 * pFieldBundaryModesList[row][bundaryPatch1] 
                                                * (uFieldBundaryModesList[column][bundaryPatch2] & inletfaceNormal)));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"KtaoD";
    writeMatrix(KtaoD, dataFile);


    // ===========================================================
    // ------ The matrix system for momentum equation ------------
    // ===========================================================

    

    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //

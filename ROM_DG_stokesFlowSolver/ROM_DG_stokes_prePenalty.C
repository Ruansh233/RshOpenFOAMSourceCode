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
    PtrList<volVectorField> gradpFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, vector>> gradpFieldBundaryModesList;

    PtrList<volVectorField> uFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, vector>> uFieldBundaryModesList;
    PtrList<volTensorField> graduFieldModesList;
    PtrList<FieldField<Foam::fvPatchField, tensor>> graduFieldBundaryModesList;

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
    fvMesh refElementMesh
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

        // grad of pressure modes
        volVectorField gradpFieldValueMode
        (
            IOobject
            (
                "grad" + pModeNames[No_],
                runTimeRef.timeName(),
                refElementMesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fvc::grad(pFieldValueMode)
        );
        gradpFieldModesList.append(gradpFieldValueMode.clone());
        gradpFieldBundaryModesList.append(gradpFieldValueMode.boundaryField().clone());
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
            MomLocalAMat(row, column) = heatConductivity * gSum(scalarField (graduFieldModesList[row] && graduFieldModesList[column]
                                                            * mesh.V()));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MomLocalAMat";
    writeMatrix(MomLocalAMat, dataFile);

    // patch ID for reference mesh of different interface
    label boundaryPatch1 (refElementMesh.boundary().findPatchID("block1_out"));
    label boundaryPatch2 (refElementMesh.boundary().findPatchID("block1_in"));
    vector inletfaceNormal(0, 0, -1);
    vector outletfaceNormal(0, 0, 1);

    // interface contribution
    vector interfaceNormal(0, 0, 1);

    // M11
    RectangularMatrix<scalar> M11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < M11.m(); ++row)
    {
        for (label column = 0; column < M11.n(); ++column)
        {
            M11(row, column) = heatConductivity * gSum(scalarField (
                                            - 0.5 * (uFieldBundaryModesList[row][boundaryPatch1] 
                                                & graduFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            + 0.5 * epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal
                                                & uFieldBundaryModesList[column][boundaryPatch1])
                                            + xigema0 * (uFieldBundaryModesList[row][boundaryPatch1]
                                                & uFieldBundaryModesList[column][boundaryPatch1])
                                            + xigema1 * (graduFieldBundaryModesList[row][boundaryPatch1]
                                                && graduFieldBundaryModesList[column][boundaryPatch1]))
                                                * mesh.boundary()[boundaryPatch1].magSf());
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
            M22(row, column) = heatConductivity * gSum(scalarField (
                                            + 0.5 * (uFieldBundaryModesList[row][boundaryPatch2] 
                                                & graduFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal
                                                & uFieldBundaryModesList[column][boundaryPatch2])
                                            + xigema0 * (uFieldBundaryModesList[row][boundaryPatch2]
                                                & uFieldBundaryModesList[column][boundaryPatch2])
                                            + xigema1 * (graduFieldBundaryModesList[row][boundaryPatch2]
                                                && graduFieldBundaryModesList[column][boundaryPatch2]))
                                                * mesh.boundary()[boundaryPatch2].magSf());
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
            M12(row, column) = heatConductivity * gSum(scalarField (
                                            - 0.5 * (uFieldBundaryModesList[row][boundaryPatch1] 
                                                & graduFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            - 0.5 * epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal
                                                & uFieldBundaryModesList[column][boundaryPatch2])
                                            - xigema0 * (uFieldBundaryModesList[row][boundaryPatch1]
                                                & uFieldBundaryModesList[column][boundaryPatch2])
                                            - xigema1 * (graduFieldBundaryModesList[row][boundaryPatch1]
                                                && graduFieldBundaryModesList[column][boundaryPatch2]))
                                                * mesh.boundary()[boundaryPatch2].magSf());
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
            M21(row, column) = heatConductivity * gSum(scalarField (
                                            + 0.5 * (uFieldBundaryModesList[row][boundaryPatch2] 
                                                & graduFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            + 0.5 * epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal
                                                & uFieldBundaryModesList[column][boundaryPatch1])
                                            - xigema0 * (uFieldBundaryModesList[row][boundaryPatch2]
                                                & uFieldBundaryModesList[column][boundaryPatch1])
                                            - xigema1 * (graduFieldBundaryModesList[row][boundaryPatch2]
                                                && graduFieldBundaryModesList[column][boundaryPatch1]))
                                                * mesh.boundary()[boundaryPatch1].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M21";
    writeMatrix(M21, dataFile);

    // boundary penalty terms
    // MtaoD, patch-inlet, boundaryPatch2
    RectangularMatrix<scalar> MtaoD(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < MtaoD.m(); ++row)
    {
        for (label column = 0; column < MtaoD.n(); ++column)
        {
            MtaoD(row, column) = heatConductivity * gSum(scalarField (
                                - (uFieldBundaryModesList[row][boundaryPatch2] 
                                    & graduFieldBundaryModesList[column][boundaryPatch2] & inletfaceNormal)
                                + epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch2] & inletfaceNormal
                                    & uFieldBundaryModesList[column][boundaryPatch2])
                                + xigema0 * (uFieldBundaryModesList[row][boundaryPatch2]
                                    & uFieldBundaryModesList[column][boundaryPatch2]))
                                    * mesh.boundary()[boundaryPatch2].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MtaoD";
    writeMatrix(MtaoD, dataFile);

    // FtaoD, patch-inlet, boundaryPatch2
    RectangularMatrix<scalar> FtaoD(modesNum, 1, Foam::Zero);
    for (label row = 0; row < FtaoD.m(); ++row)
    {
        FtaoD(row, 0) = heatConductivity * gSum(scalarField (
                            epsilonPara * (graduFieldBundaryModesList[row][boundaryPatch2] & inletfaceNormal
                            & Uin)
                            + xigema0 * (uFieldBundaryModesList[row][boundaryPatch2]
                            & Uin))
                            * mesh.boundary()[boundaryPatch2].magSf());
    }
    dataFile = runTime.globalPath()/"SVD"/"FtaoD";
    writeMatrix(FtaoD, dataFile);

    // FtaoN, patch-outlelt, boundaryPatch1, it is 0


    // ===========================================================
    // ------ The pressure gradient term -------------------------
    // ===========================================================
    // volumtric contribution, \nabla p
    RectangularMatrix<scalar> MomLocalBMat(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < MomLocalBMat.m(); ++row)
    {
        for (label column = 0; column < MomLocalBMat.n(); ++column)
        {
            MomLocalBMat(row, column) = gSum(scalarField (gradpFieldModesList[column] & uFieldModesList[row]
                                                            * mesh.V()));
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
            //                                 + 0.5 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
            //                                 + 0.5 * epsilonPara * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)));
            
            // N11(row, column) = gSum(scalarField ( + 0.5 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)));   

            N11(row, column) = gSum(scalarField (
                                            + xigema2 * pFieldBundaryModesList[column][boundaryPatch1] 
                                                * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
                                            + xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal
                                                * (graduFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch1].magSf());                                 
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
            //                                 - 0.5 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)));

            // N22(row, column) = gSum(scalarField ( - 0.5 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)));

            N22(row, column) = gSum(scalarField (
                                            + xigema2 * pFieldBundaryModesList[column][boundaryPatch2] 
                                                * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
                                            + xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal
                                                * (graduFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch2].magSf());    
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
            //                                 + 0.5 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)));

            // N12(row, column) = gSum(scalarField ( + 0.5 * pFieldBundaryModesList[column][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)));

            N12(row, column) = gSum(scalarField (
                                            - xigema2 * pFieldBundaryModesList[column][boundaryPatch2] 
                                                * (uFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal)
                                            - xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal
                                                * (graduFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch2].magSf());   
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
            //                                 - 0.5 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)));

            // N21(row, column) = gSum(scalarField ( - 0.5 * pFieldBundaryModesList[column][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal))); 

            N21(row, column) = gSum(scalarField (
                                            - xigema2 * pFieldBundaryModesList[column][boundaryPatch1] 
                                                * (uFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal)
                                            - xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal
                                                * (graduFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch1].magSf());                                
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"N21";
    writeMatrix(N21, dataFile);

    // // boundary penalty terms
    // // NtaoD, patch-outlet (pressure outlet is 0), boundaryPatch1
    // RectangularMatrix<scalar> NtaoD(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < NtaoD.m(); ++row)
    // {
    //     for (label column = 0; column < NtaoD.n(); ++column)
    //     {
    //         NtaoD(row, column) = gSum(scalarField (xigema2 * pFieldBundaryModesList[column][boundaryPatch1] 
    //                                             * (uFieldBundaryModesList[row][boundaryPatch1] & outletfaceNormal))
    //                                             * mesh.boundary()[boundaryPatch1].magSf());
    //     }
    // }
    // dataFile = runTime.globalPath()/"SVD"/"NtaoD";
    // writeMatrix(NtaoD, dataFile);

    // // boundary penalty terms
    // // NtaoN, patch-inlet (zero pressure gradient inlet), boundaryPatch2
    // RectangularMatrix<scalar> NtaoN(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < NtaoN.m(); ++row)
    // {
    //     for (label column = 0; column < NtaoN.n(); ++column)
    //     {
    //         NtaoN(row, column) = gSum(scalarField (xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch2] & inletfaceNormal)
    //                                             * (graduFieldBundaryModesList[row][boundaryPatch2] & inletfaceNormal & inletfaceNormal))
    //                                             * mesh.boundary()[boundaryPatch2].magSf());
    //     }
    // }
    // dataFile = runTime.globalPath()/"SVD"/"NtaoN";
    // writeMatrix(NtaoN, dataFile);

    // GtaoD, patch-outlelt, boundaryPatch1, it is 0

    // GtaoN, patch-inlelt, boundaryPatch2, it is 0

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

    // MomGlobalAMat should contain heatConductivity
    // MomGlobalAMat = heatConductivity * MomGlobalAMat;
    dataFile = mesh.time().path()/"SVD"/"MomGlobalAMat";
    writeMatrix(MomGlobalAMat, dataFile);

    // ===========================================================
    // --------------- MomGlobalBMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            // MomGlobalBMat(row, column) = MomLocalBMat(row, column) + NtaoN(row, column) + N11(row, column);
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
            // MomGlobalBMat(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  MomLocalBMat(row, column) 
            //                                                                + N22(row, column) + NtaoD(row, column);
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
    RectangularMatrix<scalar> ConGlobalLMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
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
            ConLocalBMat(row, column) = gSum(scalarField (gradpFieldModesList[row] & uFieldModesList[column]
                                                            * mesh.V()));
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
            // K11(row, column) = gSum(scalarField (
            //                                 + 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
            //                                 + 0.5 * epsilonPara * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal))); 

            // K11(row, column) = gSum(scalarField (+ 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)));     

            K11(row, column) = gSum(scalarField (
                                            - 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            + xigema4 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            + xigema5 * (gradpFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal
                                                * (graduFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch1].magSf());         
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
            // K22(row, column) = gSum(scalarField (
            //                                 - 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
            //                                 + xigema0 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)));

            // K22(row, column) = gSum(scalarField (- 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)));

            K22(row, column) = gSum(scalarField (
                                            + 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            + xigema4 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            + xigema5 * (gradpFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal
                                                * (graduFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch2].magSf());  
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
            // K12(row, column) = gSum(scalarField (
            //                                 + 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)));

            // K12(row, column) = gSum(scalarField (+ 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)));

            K12(row, column) = gSum(scalarField (
                                            + 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            - xigema4 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
                                            - xigema5 * (gradpFieldBundaryModesList[row][boundaryPatch1] & interfaceNormal
                                                * (graduFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch2].magSf());  
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
            // K21(row, column) = gSum(scalarField (
            //                                 - 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
            //                                 - 0.5 * epsilonPara * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)
            //                                 - xigema0 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)));

            // K21(row, column) = gSum(scalarField (- 0.5 * pFieldBundaryModesList[row][boundaryPatch1] 
            //                                     * (uFieldBundaryModesList[column][boundaryPatch2] & interfaceNormal)));

            K21(row, column) = gSum(scalarField (
                                            - 0.5 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            - xigema4 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal)
                                            - xigema5 * (gradpFieldBundaryModesList[row][boundaryPatch2] & interfaceNormal
                                                * (graduFieldBundaryModesList[column][boundaryPatch1] & interfaceNormal & interfaceNormal)))
                                                * mesh.boundary()[boundaryPatch1].magSf());                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"K21";
    writeMatrix(K21, dataFile);

    // // boundary penalty terms
    // // KtaoD, patch-inlet, boundaryPatch2
    // RectangularMatrix<scalar> KtaoD(modesNum, modesNum, Foam::Zero);

    // for (label row = 0; row < KtaoD.m(); ++row)
    // {
    //     for (label column = 0; column < KtaoD.n(); ++column)
    //     {
    //         KtaoD(row, column) = gSum(scalarField (xigema4 * pFieldBundaryModesList[row][boundaryPatch2] 
    //                                             * (uFieldBundaryModesList[column][boundaryPatch2] & inletfaceNormal))
    //                                             * mesh.boundary()[boundaryPatch2].magSf());
    //     }
    // }
    // dataFile = runTime.globalPath()/"SVD"/"KtaoD";
    // writeMatrix(KtaoD, dataFile);

    // boundary penalty terms
    // KtaoN, patch-outlet, boundaryPatch1
    RectangularMatrix<scalar> KtaoN(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < KtaoN.m(); ++row)
    {
        for (label column = 0; column < KtaoN.n(); ++column)
        {
            KtaoN(row, column) = gSum(scalarField (
                                                - pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & outletfaceNormal))
                                                * mesh.boundary()[boundaryPatch1].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"KtaoN";
    writeMatrix(KtaoN, dataFile);

    // boundary penalty terms on the right side
    // QtaoD, patch-inlet, boundaryPatch2
    RectangularMatrix<scalar> QtaoD(modesNum, 1, Foam::Zero);

    for (label row = 0; row < QtaoD.m(); ++row)
    {
        for (label column = 0; column < QtaoD.n(); ++column)
        {
            QtaoD(row, column) = gSum(scalarField (
                                                  pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (Uin & inletfaceNormal))
                                                * mesh.boundary()[boundaryPatch2].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"QtaoD";
    writeMatrix(QtaoD, dataFile);

    // boundary penalty terms on the right side
    // QtaoN, patch-outlet, boundaryPatch1, it is zero


    // ===========================================================
    // ------------------ pressure penalty term ------------------
    // ===========================================================
    // 
    RectangularMatrix<scalar> L11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            L11(row, column) = gSum(scalarField (
                                                xigema6 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * pFieldBundaryModesList[column][boundaryPatch1])
                                                * mesh.boundary()[boundaryPatch1].magSf());                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"L11";
    writeMatrix(L11, dataFile);

    RectangularMatrix<scalar> L22(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            L22(row, column) = gSum(scalarField (
                                                xigema6 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * pFieldBundaryModesList[column][boundaryPatch2])
                                                * mesh.boundary()[boundaryPatch2].magSf());                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"L22";
    writeMatrix(L22, dataFile);

    RectangularMatrix<scalar> L12(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            L12(row, column) = gSum(scalarField (
                                                xigema6 * pFieldBundaryModesList[row][boundaryPatch1] 
                                                * pFieldBundaryModesList[column][boundaryPatch2])
                                                * mesh.boundary()[boundaryPatch2].magSf());                    
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"L12";
    writeMatrix(L12, dataFile);

    RectangularMatrix<scalar> L21(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            L21(row, column) = gSum(scalarField (
                                                xigema6 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * pFieldBundaryModesList[column][boundaryPatch1])
                                                * mesh.boundary()[boundaryPatch1].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"L21";
    writeMatrix(L21, dataFile);



    // ===========================================================
    // --------------- ConGlobalBMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            // ConGlobalBMat(row, column) = ConLocalBMat(row, column) + K11(row, column) + KtaoD(row, column);
            ConGlobalBMat(row, column) = ConLocalBMat(row, column) + K11(row, column);
        }
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalBMat(row+elementI*modesNum, column+elementI*modesNum) =  ConLocalBMat(row, column) 
                                                                        + K11(row, column) + K22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            ConGlobalBMat(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  ConLocalBMat(row, column) 
                                                                           + K22(row, column) + KtaoN(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalBMat(row+elementI*modesNum, column+(elementI+1)*modesNum) =  K12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalBMat(row+(elementI+1)*modesNum, column+elementI*modesNum) =  K21(row, column);
            }
        }
    }
    dataFile = mesh.time().path()/"SVD"/"ConGlobalBMat";
    writeMatrix(ConGlobalBMat, dataFile);


    // ===========================================================
    // --------------- ConGlobalLMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            ConGlobalLMat(row, column) = L11(row, column);
        }
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalLMat(row+elementI*modesNum, column+elementI*modesNum) =  L11(row, column) + L22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            ConGlobalLMat(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) = L22(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalLMat(row+elementI*modesNum, column+(elementI+1)*modesNum) =  L12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                ConGlobalLMat(row+(elementI+1)*modesNum, column+elementI*modesNum) =  L21(row, column);
            }
        }
    }
    dataFile = mesh.time().path()/"SVD"/"ConGlobalLMat";
    writeMatrix(ConGlobalLMat, dataFile);


    // ===========================================================
    // --------------- ConGlobalFMat assignment ------------------
    // ===========================================================
    for (label row = 0; row < modesNum; ++row)
    {
        ConGlobalFMat(row, 0) = QtaoD(row, 0);
    }

    dataFile = mesh.time().path()/"SVD"/"ConGlobalFMat";
    writeMatrix(ConGlobalFMat, dataFile);


    // ===========================================================
    // ------ The matrix system for stokes flow assignment -------
    // ------ Momentum equations and continuous equations --------
    // ===========================================================

    RectangularMatrix<scalar> GlobalAMat(modesNum * elementNum * 2, modesNum * elementNum * 2, Foam::Zero);
    RectangularMatrix<scalar> GlobalFMat(modesNum * elementNum * 2, 1, Foam::Zero);

    for (label row = 0; row < modesNum * elementNum; ++row)
    {
        for (label column = 0; column < modesNum * elementNum; ++column)
        {
            GlobalAMat(row, column) =  MomGlobalAMat(row, column);
        }
    }

    for (label row = 0; row < modesNum * elementNum; ++row)
    {
        for (label column = 0; column < modesNum * elementNum; ++column)
        {
            GlobalAMat(row, column+elementNum*modesNum) =  MomGlobalBMat(row, column);
        }
    }

    for (label row = 0; row < modesNum * elementNum; ++row)
    {
        for (label column = 0; column < modesNum * elementNum; ++column)
        {
            GlobalAMat(row+elementNum*modesNum, column) =  ConGlobalBMat(row, column);
        }
    }

    for (label row = 0; row < modesNum * elementNum; ++row)
    {
        for (label column = 0; column < modesNum * elementNum; ++column)
        {
            GlobalAMat(row+elementNum*modesNum, column+elementNum*modesNum) =  ConGlobalLMat(row, column);
        }
    }

    for (label row = 0; row < elementNum*modesNum; ++row)
    {
        GlobalFMat(row, 0) =  MomGlobalFMat(row, 0);
    }

    for (label row = 0; row < elementNum*modesNum; ++row)
    {
        GlobalFMat(row+elementNum*modesNum, 0) =  ConGlobalFMat(row, 0);
    }
    
    dataFile = mesh.time().path()/"SVD"/"GlobalAMat";
    writeMatrix(GlobalAMat, dataFile);
    
    dataFile = mesh.time().path()/"SVD"/"GlobalFMat";
    writeMatrix(GlobalFMat, dataFile);


    // ===========================================================
    // ---- solve The matrix system for stokes flow --------------
    // ---- output the calculated coefficient and snapshots ------
    // ===========================================================

    // solve matrix system
    RectangularMatrix<scalar> tempCalCoefficientM;
    RectangularMatrix<scalar> uCalCoefficientM(modesNum, elementNum);
    RectangularMatrix<scalar> pCalCoefficientM(modesNum, elementNum);
    tempCalCoefficientM = SVDinv(GlobalAMat) * GlobalFMat;

    for (label row = 0; row < uCalCoefficientM.m(); ++row)
    {
        for (label column = 0; column < uCalCoefficientM.n(); ++column)
        {
            uCalCoefficientM(row, column) = tempCalCoefficientM(row+column*modesNum, 0);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"uCalCoefficientM";
    writeMatrix(uCalCoefficientM.T(), dataFile);

    for (label row = 0; row < pCalCoefficientM.m(); ++row)
    {
        for (label column = 0; column < pCalCoefficientM.n(); ++column)
        {
            pCalCoefficientM(row, column) = tempCalCoefficientM(row+(column+elementNum)*modesNum, 0);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"pCalCoefficientM";
    writeMatrix(pCalCoefficientM.T(), dataFile);

    // ===========================================================
    // ------------ error of pressure field ----------------------
    // =========================================================== 
    // read the matrix
    // read cell pressure modes matrix
    dataFile = dataPath/"pModesM";
    RectangularMatrix<scalar> pModesM(mesh.C().size(), modesNum);

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
                pModesM(row, modesNo) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }

    // calculate pressure snapshots
    RectangularMatrix<scalar> pCalSnapshotsM;
    pCalSnapshotsM = pModesM * pCalCoefficientM;
    dataFile = mesh.time().path()/"SVD"/"pCalSnapshotsM";
    writeMatrix(pCalSnapshotsM, dataFile);

    // The snapshots matrix of pressure
    dataFile = dataPath/"pSnapshotsM";
    RectangularMatrix<scalar> pSnapshotsM(mesh.C().size(), elementNum);
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
                pSnapshotsM(row, elementI) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }

    // pressure error matrix
    RectangularMatrix<scalar> pErrorM(mesh.C().size(), elementNum);
    for (label row = 0; row < pErrorM.m(); ++row)
    {
        for (label column = 0; column < pErrorM.n(); ++column)
        {
            pErrorM(row, column) =  (pCalSnapshotsM(row, column) - pSnapshotsM(row, column))/pSnapshotsM(row, column);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"pErrorM";
    writeMatrix(pErrorM, dataFile);    


    // ===========================================================
    // ------------ error of velocity field ----------------------
    // =========================================================== 
    // read cell velocity modes matrix
    dataFile = dataPath/"uModesM";
    RectangularMatrix<scalar> uModesM(mesh.C().size() * 3, modesNum);

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
                uModesM(row, modesNo) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }

    // calculate velocity snapshots
    RectangularMatrix<scalar> uCalSnapshotsM;
    uCalSnapshotsM = uModesM * uCalCoefficientM;
    dataFile = mesh.time().path()/"SVD"/"uCalSnapshotsM";
    writeMatrix(uCalSnapshotsM, dataFile);

    // The snapshots matrix of pressure
    dataFile = dataPath/"uSnapshotsM";
    RectangularMatrix<scalar> uSnapshotsM(mesh.C().size()*3, elementNum);
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
                uSnapshotsM(row, elementI) = singleData.scalarToken();
            }   
            ++row;
        }                       
    }  
    else
    {
        Info << "file: " << dataFile << " is not exist!" << endl;
        // break;
    }

    // pressure error matrix
    RectangularMatrix<scalar> uErrorM(mesh.C().size() * 3, elementNum);
    for (label row = 0; row < uErrorM.m(); ++row)
    {
        for (label column = 0; column < uErrorM.n(); ++column)
        {
            uErrorM(row, column) =  (uCalSnapshotsM(row, column) - uSnapshotsM(row, column))/uSnapshotsM(row, column);
        }
    }
    dataFile = mesh.time().path()/"SVD"/"uErrorM";
    writeMatrix(uErrorM, dataFile);    


    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //

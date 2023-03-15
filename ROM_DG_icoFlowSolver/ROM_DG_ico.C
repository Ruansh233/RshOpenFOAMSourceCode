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
    RectangularMatrix<scalar> MomGlobalCMat(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
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

    // boundary penalty terms
    // NtaoD, patch-outlet (pressure outlet is 0), boundaryPatch1
    RectangularMatrix<scalar> NtaoD(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < NtaoD.m(); ++row)
    {
        for (label column = 0; column < NtaoD.n(); ++column)
        {
            NtaoD(row, column) = gSum(scalarField (xigema2 * pFieldBundaryModesList[column][boundaryPatch1] 
                                                * (uFieldBundaryModesList[row][boundaryPatch1] & outletfaceNormal))
                                                * mesh.boundary()[boundaryPatch1].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"NtaoD";
    writeMatrix(NtaoD, dataFile);

    // boundary penalty terms
    // NtaoN, patch-inlet (zero pressure gradient inlet), boundaryPatch2
    RectangularMatrix<scalar> NtaoN(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < NtaoN.m(); ++row)
    {
        for (label column = 0; column < NtaoN.n(); ++column)
        {
            NtaoN(row, column) = gSum(scalarField (xigema3 * (gradpFieldBundaryModesList[column][boundaryPatch2] & inletfaceNormal)
                                                * (graduFieldBundaryModesList[row][boundaryPatch2] & inletfaceNormal & inletfaceNormal))
                                                * mesh.boundary()[boundaryPatch2].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"NtaoN";
    writeMatrix(NtaoN, dataFile);

    // ===========================================================
    // ------------- The convection term -------------------------
    // ===========================================================
    // \vec{u} \cdot \nabla \vec{u}, containing tensor product
    // volumtric contribution, u'u:\nabla u and uu':\nabla u 
    RectangularMatrix<scalar> MomLocalCMat(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                MomLocalCMat(row + modesNum * projI, column) = gSum(scalarField ( 
                                                                - ((uFieldModesList[row] * uFieldModesList[column])
                                                                && graduFieldModesList[projI])
                                                                * mesh.V()));
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"MomLocalCMat";
    writeMatrix(MomLocalCMat, dataFile); 

    // surface contribution, J11
    RectangularMatrix<scalar> J11(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                J11(row + modesNum * projI, column) = gSum(scalarField (
                                                    0.5 * (uFieldBundaryModesList[row][boundaryPatch1] 
                                                    * uFieldBundaryModesList[column][boundaryPatch1]) 
                                                    & uFieldBundaryModesList[projI][boundaryPatch1] & interfaceNormal)
                                                    * mesh.boundary()[boundaryPatch1].magSf());                             
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"J11";
    writeMatrix(J11, dataFile);

    // surface contribution, J22
    RectangularMatrix<scalar> J22(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                J22(row + modesNum * projI, column) = gSum(scalarField (
                                                    - 0.5 * ((uFieldBundaryModesList[row][boundaryPatch2] 
                                                    * uFieldBundaryModesList[column][boundaryPatch2]) 
                                                    & uFieldBundaryModesList[projI][boundaryPatch2] & interfaceNormal))
                                                    * mesh.boundary()[boundaryPatch2].magSf());                             
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"J22";
    writeMatrix(J22, dataFile);

    // surface contribution, J12
    RectangularMatrix<scalar> J12(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                J12(row + modesNum * projI, column) = gSum(scalarField (
                                                    0.5 * ((uFieldBundaryModesList[row][boundaryPatch2] 
                                                    * uFieldBundaryModesList[column][boundaryPatch2]) 
                                                    & uFieldBundaryModesList[projI][boundaryPatch1] & interfaceNormal))
                                                    * mesh.boundary()[boundaryPatch2].magSf());                             
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"J12";
    writeMatrix(J12, dataFile);

    // surface contribution, J21
    RectangularMatrix<scalar> J21(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                J21(row + modesNum * projI, column) = gSum(scalarField (
                                                    - 0.5 * ((uFieldBundaryModesList[row][boundaryPatch1] 
                                                    * uFieldBundaryModesList[column][boundaryPatch1]) 
                                                    & uFieldBundaryModesList[projI][boundaryPatch2] & interfaceNormal))
                                                    * mesh.boundary()[boundaryPatch1].magSf());                             
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"J21";
    writeMatrix(J21, dataFile);

    // surface contribution, JtaoN
    // patch outlet, boundaryPatch1
    RectangularMatrix<scalar> JtaoN(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                JtaoN(row + modesNum * projI, column) = gSum(scalarField (
                                                    0.5 * ((uFieldBundaryModesList[row][boundaryPatch1] 
                                                    * uFieldBundaryModesList[column][boundaryPatch1]) 
                                                    & uFieldBundaryModesList[projI][boundaryPatch1] & outletfaceNormal))
                                                    * mesh.boundary()[boundaryPatch1].magSf());                             
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"JtaoN";
    writeMatrix(JtaoN, dataFile);

    // volumetric contribution on right side, GLocal
    // volumtric contribution, u'u:\nabla u and uu':\nabla u 
    RectangularMatrix<scalar> GLocal(modesNum * modesNum, modesNum, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                GLocal(row + modesNum * projI, column) = gSum(scalarField ( 
                                                        - (uFieldModesList[row] & graduFieldModesList[projI]
                                                        & uFieldModesList[column])
                                                        * mesh.V()));
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"GLocal";
    writeMatrix(GLocal, dataFile); 

    // surface contribution, GtaoD
    // patch inlet, boundaryPatch2
    RectangularMatrix<scalar> GtaoD(modesNum, 1, Foam::Zero);
    for (label projI = 0; projI < modesNum; ++projI)
    {
        GtaoD(projI, 0) = gSum(scalarField ( - ((Uin * Uin) 
                                            & uFieldBundaryModesList[projI][boundaryPatch2] & inletfaceNormal))
                                            * mesh.boundary()[boundaryPatch2].magSf()); 
    }
    dataFile = runTime.globalPath()/"SVD"/"GtaoD";
    writeMatrix(GtaoD, dataFile);

    // ===========================================================
    // ------ The matrix system for continuous equation ------------
    // ===========================================================

    // ===========================================================
    // ------ The divergence of velocity -------------------------
    // ===========================================================
    // volumtric contribution, \nabla \cdot v
    RectangularMatrix<scalar> ConLocalEMat(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < ConLocalEMat.m(); ++row)
    {
        for (label column = 0; column < ConLocalEMat.n(); ++column)
        {
            ConLocalEMat(row, column) = gSum(scalarField (gradpFieldModesList[row] & uFieldModesList[column]
                                                            * mesh.V()));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"ConLocalEMat";
    writeMatrix(ConLocalEMat, dataFile);

    // K11
    RectangularMatrix<scalar> K11(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < K11.m(); ++row)
    {
        for (label column = 0; column < K11.n(); ++column)
        {
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

    // boundary penalty terms
    // KtaoD, patch-inlet, boundaryPatch2
    RectangularMatrix<scalar> KtaoD(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < KtaoD.m(); ++row)
    {
        for (label column = 0; column < KtaoD.n(); ++column)
        {
            KtaoD(row, column) = gSum(scalarField (xigema4 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (uFieldBundaryModesList[column][boundaryPatch2] & inletfaceNormal))
                                                * mesh.boundary()[boundaryPatch2].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"KtaoD";
    writeMatrix(KtaoD, dataFile);

    // boundary penalty terms
    // KtaoN, patch-outlet, boundaryPatch1
    RectangularMatrix<scalar> KtaoN(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < KtaoN.m(); ++row)
    {
        for (label column = 0; column < KtaoN.n(); ++column)
        {
            KtaoN(row, column) = gSum(scalarField (
                                                - pFieldBundaryModesList[row][boundaryPatch1] 
                                                * (uFieldBundaryModesList[column][boundaryPatch1] & outletfaceNormal)
                                                + xigema5 * (gradpFieldBundaryModesList[row][boundaryPatch1] & outletfaceNormal
                                                * (graduFieldBundaryModesList[column][boundaryPatch1] & outletfaceNormal & outletfaceNormal)))
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
                                                * (Uin & inletfaceNormal)
                                                + xigema4 * pFieldBundaryModesList[row][boundaryPatch2] 
                                                * (Uin & inletfaceNormal))
                                                * mesh.boundary()[boundaryPatch2].magSf());
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"QtaoD";
    writeMatrix(QtaoD, dataFile);

    // boundary penalty terms on the right side
    // QtaoN, patch-outlet, boundaryPatch1, it is zero

    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //

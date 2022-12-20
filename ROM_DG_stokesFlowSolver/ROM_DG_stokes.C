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
    fileName refCaseName(svdDict.getWord("refCaseName"));
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

        // grad of pressure modes
        volVectorField pFieldValueModegrad
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
        gradpFieldModesList.append(pFieldValueModegrad.clone());
        gradpFieldBundaryModesList.append(pFieldValueModegrad.boundaryField().clone());
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
        graduFieldBundaryModesList.append(uFieldValueMode.boundaryField().clone());
    }

    // ===========================================================
    // ------ The matrix system ----------------------------------
    // ===========================================================
    // create Matrix system
    // initial global matrix
    RectangularMatrix<vector> globalAMatrix(modesNum * elementNum, modesNum * elementNum, Foam::Zero);
    RectangularMatrix<vector> globalBMmatrix(modesNum * elementNum, 1, Foam::Zero);

    // volumtric contribution
    RectangularMatrix<vector> localAMatrix(modesNum, modesNum, Foam::Zero);
    for (label row = 0; row < localAMatrix.m(); ++row)
    {
        for (label column = 0; column < localAMatrix.n(); ++column)
        {
            localAMatrix(row, column) = kineticViscosity * gSum(vectorField (gradfieldModesList[row] && gradfieldModesList[column]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"localAMatrix";
    writeMatrix(localAMatrix, dataFile);

    // interface contribution
    vector interfaceNormal(0, 0, 1);

    // M11
    RectangularMatrix<scalar> M11(modesNum, modesNum, Foam::Zero);
    label bundaryPatch1(matchPatchID[0]);
    label bundaryPatch2(matchPatchID[1]);
    label bundaryPatch3(matchPatchID[2]);


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
                                                & gradfieldBundaryModesList[column][bundaryPatch2])));
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
            M22(row, column) = gSum(scalarField (
                                              0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & interfaceNormal) 
                                            - 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            + xigema1 * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch1])));
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
                                            - 0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch2] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch1] & interfaceNormal) 
                                            - 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                                                * fieldBundaryModesList[column][bundaryPatch1]
                                            - xigema1 * (gradfieldBundaryModesList[row][bundaryPatch2]
                                                & gradfieldBundaryModesList[column][bundaryPatch1])));
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
                                              0.5 * heatConductivity * fieldBundaryModesList[row][bundaryPatch1] 
                                                * (gradfieldBundaryModesList[column][bundaryPatch2] & interfaceNormal) 
                                            + 0.5 * epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & interfaceNormal) 
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                                                * fieldBundaryModesList[column][bundaryPatch2]
                                            - xigema1 * (gradfieldBundaryModesList[row][bundaryPatch1]
                                                & gradfieldBundaryModesList[column][bundaryPatch2])));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"M21";
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
                                    * fieldBundaryModesList[column][bundaryPatch1]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"Min";
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
                                    * fieldBundaryModesList[column][bundaryPatch2]));
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"Mout";
    writeMatrix(Mout, dataFile);

    // Fin
    RectangularMatrix<scalar> Fin(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fin.m(); ++row)
    {
        Fin(row, 0) = gSum(scalarField (
                            epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch1] & inletfaceNormal) 
                            * Tin
                            + xigema0 * fieldBundaryModesList[row][bundaryPatch1]
                            * Tin));
    }

    // Fout
    RectangularMatrix<scalar> Fout(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fout.m(); ++row)
    {
        Fout(row, 0) = gSum(scalarField (
                            epsilonPara * heatConductivity * (gradfieldBundaryModesList[row][bundaryPatch2] & outletfaceNormal) 
                            * Tout
                            + xigema0 * fieldBundaryModesList[row][bundaryPatch2]
                            * Tout));
    }

    // Fn
    RectangularMatrix<scalar> Fn(modesNum, 1, Foam::Zero);
    for (label row = 0; row < Fn.m(); ++row)
    {
        Fn(row, 0) = gSum(scalarField (
                            heatConductivity * fieldBundaryModesList[row][bundaryPatch3] * qn));
    }


    // global matrix
    // global phi matrix
    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            globalAMatrix(row, column) = localAMatrix(row, column) + Min(row, column) + M11(row, column);
        }
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalAMatrix(row+elementI*modesNum, column+elementI*modesNum) =  localAMatrix(row, column) 
                                                                        + M11(row, column) + M22(row, column);
            }
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        for (label column = 0; column < modesNum; ++column)
        {
            globalAMatrix(row+(elementNum-1)*modesNum, column+(elementNum-1)*modesNum) =  localAMatrix(row, column) 
                                                                           + M22(row, column) + Mout(row, column);
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalAMatrix(row+elementI*modesNum, column+(elementI+1)*modesNum) =  M12(row, column);
            }
        }
    }

    for(label elementI = 0; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            for (label column = 0; column < modesNum; ++column)
            {
                globalAMatrix(row+(elementI+1)*modesNum, column+elementI*modesNum) =  M21(row, column);
            }
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"globalAMatrix";
    writeMatrix(globalAMatrix, dataFile);


    // global F matrix
    for (label row = 0; row < modesNum; ++row)
    {
        // globalBMmatrix(row, 0) = Fin(row, 0) + Fn(row, 0);
        globalBMmatrix(row, 0) = Fin(row, 0);
    }

    for(label elementI = 1; elementI < elementNum - 1; ++elementI)
    {
        for (label row = 0; row < modesNum; ++row)
        {
            globalBMmatrix(row+elementI*modesNum, 0) = Fn(row, 0);
        }
    }

    for (label row = 0; row < modesNum; ++row)
    {
        // globalBMmatrix(row+(elementNum-1)*modesNum, 0) = Fout(row, 0) + Fn(row, 0);
        globalBMmatrix(row+(elementNum-1)*modesNum, 0) = Fout(row, 0);
    }
    dataFile = runTime.globalPath()/"SVD"/"globalBMmatrix";
    writeMatrix(globalBMmatrix, dataFile);


    // solve matrix system
    RectangularMatrix<scalar> tempCalCoefficientM;
    RectangularMatrix<scalar> calCoefficientM(modesNum, elementNum);
    tempCalCoefficientM = SVDinv(globalAMatrix) * globalBMmatrix;
    for (label row = 0; row < calCoefficientM.m(); ++row)
    {
        for (label column = 0; column < calCoefficientM.n(); ++column)
        {
            calCoefficientM(row, column) =  tempCalCoefficientM(row+column*modesNum, 0);
        }
    }
    dataFile = runTime.globalPath()/"SVD"/"calCoefficientM";
    writeMatrix(calCoefficientM, dataFile);


    // calculate snapshots
    RectangularMatrix<scalar> calSnapshotsM;
    calSnapshotsM = modesM * calCoefficientM;
    dataFile = runTime.globalPath()/"SVD"/"calSnapshotsM";
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
    dataFile = runTime.globalPath()/"SVD"/"errorM";
    writeMatrix(errorM, dataFile);    


    return 0;
}


// ************************************************************************* //

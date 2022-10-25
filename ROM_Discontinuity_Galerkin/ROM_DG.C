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

    // // matchPatchID
    // List<label> matchPatchID(2*modesNum);
    // dataFile = dataPath + "/boundaryPatchID";
    // label row(0);

    // if(isFile(dataFile))
    // {                
    //     IFstream dataStream(dataFile);
    //     while(dataStream.read(matchPatchID[row]))
    //     {
    //         ++row;
    //     }            
    // }  
    // else
    // {
    //     Info << "file: " << dataFile << " is not exist!" << endl;
    //     // break;
    // }
    // // Info << matchPatchID << endl;


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
    }

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

        forAll(mesh.C(), cellI)
        {
            gradModesM[cellI][No_] = fieldValueModegrad[cellI];
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

    // create Matrix system
    // initial global matrix
    RectangularMatrix<scalar> globalphiMatrix(modesNum * modesNum, modesNum * modesNum, Foam::Zero);
    RectangularMatrix<scalar> globalFM(modesNum * modesNum, 1, Foam::Zero);

    // volumtric contribution
    RectangularMatrix<scalar> localphiMatrix(modesNum, modesNum, Foam::Zero);

    for (label row = 0; row < localphiMatrix.m(); ++row)
    {
        for (label column = 0; column < localphiMatrix.n(); ++column)
        {
            localphiMatrix(row, column) = gradModesM(column, row) & gradModesM(row, column);
        }
    }

    // interface contribution
    // M11

    // M22

    // M12

    // M12

    // boundary penalty terms

    // global matrix

    

    // solve matrix system


 


    

    return 0;
}


// ************************************************************************* //

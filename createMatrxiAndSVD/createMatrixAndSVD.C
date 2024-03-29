/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
// #include "emptyPolyPatch.H"
#include "IFstream.H"
// #include "stringOps.H"
#include "SVD.H"
// #include "writeMatrix.H"

int main(int argc, char *argv[])
{
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (mesh.time().path()/"0.5");

    fileName fieldName("p");    

    volScalarField p
    (
        IOobject
        (
            fieldName,
            dataPath,
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // fileName dataFile (dataPath/"pMatrix");
    // OFstream outputFile_(dataFile);

    RectangularMatrix<scalar> testM(mesh.C().size(), 2);

    forAll(mesh.C(), cellI)
    {
        testM[cellI][0] = p[cellI];
        testM[cellI][1] = p[cellI];
    }

    // for (label row = 0; row < testM.m(); ++row)
    // {
    //     for (label column = 0; column < testM.n(); ++column)
    //     {
    //         outputFile_.width(16);
    //         outputFile_ << testM[row][column];
    //     }
    //     outputFile_ << endl;
    // }
    
    // forAll(testM.subColumn(0), rows)
    // {        
    //     forAll(testM.subRow(0), columns)
    //     {
    //         outputFile_.width(16);
    //         outputFile_ << testM[rows][columns];
    //     }
    //     outputFile_ << endl;
    // }

    SVD UField(testM);
    // Info << "UField.U().m() " << UField.U().m() << endl
    //      << "UField.U().n() " << UField.U().n() << endl
    //      << "UField.V().m() " << UField.V().m() << endl
    //      << "UField.V().n() " << UField.V().n() << endl
    //      << "UField.S().size() " << UField.S().size() << endl;

    // Info << "UField.U() " << endl << UField.U() << endl
    //      << "UField.V() " << endl << UField.V() << endl;
    // Info << "UField.S() " << endl << UField.S() << endl; 
    // Info << "inv(UField.S()) " << endl << inv(UField.S()) << endl;

    // for (label row = 0; row < UField.S().size(); ++row)
    // {
    //     outputFile_ << UField.S()[row];
    //     outputFile_ << endl;
    // }

	autoPtr<OFstream> outputFilePtr;
    // // Write stuff
    // outputFilePtr() << "# This is a header" << endl;
    // outputFilePtr() << "0 1 2 3 4 5" << endl;

    fileName dataFile (dataPath/"diagS");
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < UField.S().size(); ++row)
    {
        outputFilePtr() << UField.S()[row];
        outputFilePtr() << endl;
    }

    dataFile = dataPath/"UMatrix";
    outputFilePtr.reset(new OFstream(dataFile));
    for (label row = 0; row < UField.U().m(); ++row)
    {
        for (label column = 0; column < UField.U().n(); ++column)
        {
            outputFilePtr().width(16);
            outputFilePtr() << UField.U()[row][column];
        }
        outputFilePtr() << endl;
    }


    // RectangularMatrix<scalar> testM2;
    // multiply(testM2, UField.U(), UField.S(), UField.V().T());
    // multiply(testM2, UField.U(), UField.S(), UField.V());
    // testM2 = UField.U() * inv(UField.S()) * UField.V().T();
    // testM2 = testM * testM;

    // Info << "testM " << endl << testM <<endl;
    // Info << "testM2 " << endl << testM2 <<endl;


    // RectangularMatrix<scalar> U()



    // fileName outputFile (dataPath/fieldName + "_out");
    // OFstream outputStream(outputFile);

    // // outputStream << "testM[0][0]" << endl;
    // // outputStream << testM[0][0];
    // outputStream << testM;






    // List<word> modeNumber ({"0", "1", "2", "3", "4"});
    // List<word> modeNumber ({"0", "1", "2"});
    // List<word> fieldName ({"Ux", "Uy", "U", "magU", "gradp0", "gradp1", "gradp2"});

    // forAll(fieldName, nameNo)
    // {
    //     Info << "FieldName: " << fieldName[nameNo] << endl;

    //     forAll(modeNumber, No_)
    //     {
    //         fileName dataFile (dataPath/fieldName[nameNo] + "_mode" + modeNumber[No_]);
    //         fileName modeFieldName(fieldName[nameNo] + "_mode" + modeNumber[No_]);

    //         volScalarField U_mode
    //         (
    //             IOobject
    //             (
    //                 modeFieldName,
    //                 mesh.time().timeName(),
    //                 mesh,
    //                 IOobject::NO_READ,
    //                 IOobject::AUTO_WRITE
    //             ),
    //             mesh,
    //             dimVelocity
    //         );

    //         IFstream dataStream(dataFile);

    //         forAll(mesh.cellZones()[0], cellI)
    //         {
    //             forAll(mesh.cellZones(), zoneI)
    //             {
    //                 label cell = mesh.cellZones()[zoneI][cellI];
    //                 dataStream.read(U_mode[cell]);
    //             }
    //         }

    //         U_mode.write();
    //     }

    // }

    // test Matrix multiplication 
    RectangularMatrix<scalar> testM1(2, 2);
    RectangularMatrix<scalar> testM2(2, 2);

    testM1[0][0] = 1;
    testM1[1][0] = 1;
    testM1[0][1] = 0;
    testM1[1][1] = 1;

    testM2[0][0] = 1;
    testM2[1][0] = 2;
    testM2[0][1] = 1;
    testM2[1][1] = 2;

    RectangularMatrix<scalar> testM3(testM1 * testM2);
    Info << "testM3: " << testM3 << endl;



    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

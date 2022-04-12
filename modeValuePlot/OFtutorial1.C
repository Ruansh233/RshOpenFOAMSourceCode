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
#include "emptyPolyPatch.H"
#include "IFstream.H"
#include "stringOps.H"

int main(int argc, char *argv[])
{
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    fileName dataPath (mesh.time().path()/"postProcessing");

    List<word> modeNumber ({"0", "1", "2", "3", "4"});
    // List<word> modeNumber ({"0", "1", "2"});
    List<word> fieldName ({"Ux", "Uy", "Uz", "magU", "gradpx", "gradpy", "gradpz"});

    forAll(fieldName, nameNo)
    {
        Info << "FieldName: " << fieldName[nameNo] << endl;

        forAll(modeNumber, No_)
        {
            fileName dataFile (dataPath/fieldName[nameNo] + "_mode" + modeNumber[No_]);
            fileName modeFieldName(fieldName[nameNo] + "_mode" + modeNumber[No_]);

            if(isFile(dataFile))
            {
                volScalarField fieldValueMode
                (
                    IOobject
                    (
                        modeFieldName,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    dimVelocity
                );
                
                IFstream dataStream(dataFile);
                forAll(mesh.cellZones()[0], cellI)
                {
                    forAll(mesh.cellZones(), zoneI)
                    {
                        label cell = mesh.cellZones()[zoneI][cellI];
                        dataStream.read(fieldValueMode[cell]);
                    }
                }

                fieldValueMode.write();
            }     
            
        }

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

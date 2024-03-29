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
#include "wallDist.H"

int main(int argc, char *argv[])
{   
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;
    
    volTensorField volDfield
    (
        IOobject
        (
            "nonUniformDF_D",
            mesh.time().caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/sqr(dimLength),
        tensorField (mesh.C().size(), Zero)
    );


    volTensorField volFfield
    (
        IOobject
        (
            "nonUniformDF_F",
            mesh.time().caseConstant(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength,
        tensorField (mesh.C().size(), Zero)
    );  


    const word dictName("fvOptions");

    // Create and input-output object - this holds the path to the dict and its name
    IOdictionary customDict
    (
        IOobject
        (
            dictName, // name of the file
            mesh.time().caseConstant(), // path to where the file is
            mesh, // reference to the mesh needed by the constructor
            IOobject::MUST_READ // indicate that reading this dictionary is compulsory
        )
    );

    wordRe zoneName_ (customDict.subDict("blockage").subDict("explicitPorositySourceCoeffs").lookup("cellZone"));

    labelList cellZoneIDs_ = mesh.cellZones().indices(zoneName_);

    forAll(cellZoneIDs_, i)
    {
        forAll(mesh.cellZones()[cellZoneIDs_[i]], j)
        {
            const label cellI = mesh.cellZones()[cellZoneIDs_[i]][j];

            if(y_[cellI] > 0.5e-3)
            {
                volDfield[cellI].xx() = y_[cellI]*1.0e6;
                volDfield[cellI].yy() = y_[cellI]*1.0e6;
                volDfield[cellI].zz() = y_[cellI]*1.0e6;

                volFfield[cellI].xx() = y_[cellI]*1.0e5;
                volFfield[cellI].yy() = y_[cellI]*1.0e5;
                volFfield[cellI].zz() = y_[cellI]*1.0e5;
            }   
        }
    }

    volDfield.write();
    volFfield.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

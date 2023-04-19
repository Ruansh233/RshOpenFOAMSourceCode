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

    volScalarField y_ = wallDist(mesh).y();

    y_.write();

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

    forAll(volDfield, cellI)
    {
        if(y_[cellI] > 0.5e-3)
        {
            volDfield[cellI].xx() = 1.0e4;
            volDfield[cellI].yy() = 1.0e4;
            volDfield[cellI].zz() = 1.0e4;
        } 
    }
    volDfield.write();


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

    forAll(volFfield, cellI)
    {
        if(y_[cellI] > 0.5e-3)
        {
            volFfield[cellI].xx() = 1.0e3;
            volFfield[cellI].yy() = 1.0e3;
            volFfield[cellI].zz() = 1.0e3;
        }   
    }
    volFfield.write();


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

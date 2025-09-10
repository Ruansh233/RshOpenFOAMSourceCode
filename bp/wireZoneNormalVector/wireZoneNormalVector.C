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

int main(int argc, char *argv[])
{   
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    Info<< "Start\n" << endl;

    #include "readDict.H"
    
    volVectorField nnV
    (
        IOobject
        (
            "nnV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    );

    volVectorField ntV
    (
        IOobject
        (
            "ntV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    );

    volVectorField npnV
    (
        IOobject
        (
            "npnV",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimless
    );

    nnV = vector::zero;
    ntV = vector::zero;
    npnV = vector::zero;

    const scalar rodWireR(0.5*(rodD+wireD));
    const scalar inclAng(Foam::atan((2*constant::mathematical::pi*rodWireR)/wireH));

    scalar theta;

    const label wireZoneID (mesh.cellZones().findZoneID(cellZoneName));

    forAll(mesh.cellZones()[wireZoneID], i)
    {
        const label cellI (mesh.cellZones()[wireZoneID][i]);

        theta = (2*constant::mathematical::pi*mesh.C()[cellI].z()-refHight)/wireH+refAngle;
            
        nnV[cellI].x() = Foam::cos(inclAng)*Foam::cos(theta-constant::mathematical::pi/2);
        nnV[cellI].y() = Foam::cos(inclAng)*Foam::sin(theta-constant::mathematical::pi/2);
        nnV[cellI].z() = Foam::sin(inclAng);

        ntV[cellI].x() = Foam::sin(inclAng)*Foam::cos(theta+constant::mathematical::pi/2);
        ntV[cellI].y() = Foam::sin(inclAng)*Foam::sin(theta+constant::mathematical::pi/2);
        ntV[cellI].z() = Foam::cos(inclAng);

        npnV[cellI].x() = -Foam::cos(theta);
        npnV[cellI].y() = -Foam::sin(theta);
        npnV[cellI].z() = 0.0; 
    }

    forAll(nnV.boundaryField(), patchI)
    {
        if(nnV.boundaryField()[patchI].type() == "calculated")
        {
            nnV.boundaryFieldRef().set(patchI, 
            fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], nnV));

            ntV.boundaryFieldRef().set(patchI, 
                fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], ntV));

            npnV.boundaryFieldRef().set(patchI,
                fvPatchField<vector>::New("zeroGradient", mesh.boundary()[patchI], npnV));
        }
    }

    // // Rsh, 2024-01-25, replaced by mesh.time().write() for bwunicluster running
    nnV.write();
    ntV.write();
    npnV.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

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
            "nonUniDcoeff",
            // mesh.time().caseConstant(),
            mesh.time().timeName(),
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
            "nonUniFcoeff",
            // mesh.time().caseConstant(),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless/dimLength,
        tensorField (mesh.C().size(), Zero)
    );  


    const scalar rodR(8.2e-3);
    const scalar wireR(2.2e-3);
    const scalar rodWireP(0.5*(rodR+wireR));
    const scalar wireH(300.0e-3);

    const scalar x0(0.0);
    const scalar y0(0.0);
    
    scalar xw;
    scalar yw;
    scalar theta;

    const scalar D0(500);
    const scalar F0(200);

    forAll(mesh.C(), cellI)
    {
        theta = mesh.C()[cellI].z()/wireH*2*3.14;

        scalar cellx = mesh.C()[cellI].x();
        scalar celly = mesh.C()[cellI].y();

        xw = x0 + rodWireP*Foam::cos(theta);
        yw = y0 + rodWireP*Foam::sin(theta);

        if(Foam::sqrt(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)) <  wireR/2)
        {
            volDfield[cellI].xx() = D0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volDfield[cellI].yy() = D0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volDfield[cellI].zz() = D0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));

            volFfield[cellI].xx() = F0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volFfield[cellI].yy() = F0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volFfield[cellI].zz() = F0*Foam::exp(-5.0e6*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
        }
    }

    forAll(volDfield.boundaryField(), patchI)
    {
        if(volDfield.boundaryField()[patchI].type() == "calculated")
        {
            volDfield.boundaryFieldRef().set(patchI, 
            fvPatchField<tensor>::New("zeroGradient", mesh.boundary()[patchI], volDfield));

            volFfield.boundaryFieldRef().set(patchI, 
                fvPatchField<tensor>::New("zeroGradient", mesh.boundary()[patchI], volFfield));
        }
    }

    volDfield.write();
    volFfield.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

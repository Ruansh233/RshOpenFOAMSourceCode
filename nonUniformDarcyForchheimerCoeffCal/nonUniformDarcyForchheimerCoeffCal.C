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
    const scalar wireInfluR(3.3e-3);
    const scalar rodWireP(0.5*(rodR+wireInfluR));
    const scalar wireH(300.0e-3);

    const scalar x0(0.0);
    const scalar y0(0.0);
    
    scalar xw;
    scalar yw;
    scalar theta;

    const scalar D0(500);
    const scalar F0(105);

    const scalar xigma = -2.0e6;

    scalar wireResistance = 0.0;
    scalar wireVolume = 0.0;

    forAll(mesh.C(), cellI)
    {
        theta = mesh.C()[cellI].z()/wireH*2*constant::mathematical::pi;

        scalar cellx = mesh.C()[cellI].x();
        scalar celly = mesh.C()[cellI].y();

        xw = x0 + rodWireP*Foam::cos(theta);
        yw = y0 + rodWireP*Foam::sin(theta);

        if(Foam::sqrt(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)) <  wireInfluR/2)
        {
            volDfield[cellI].xx() = D0*Foam::exp(xigma*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volDfield[cellI].yy() = volDfield[cellI].xx();
            volDfield[cellI].zz() = volDfield[cellI].xx();

            volFfield[cellI].xx() = F0*Foam::exp(xigma*(Foam::sqr(cellx-xw) + Foam::sqr(celly-yw)));
            volFfield[cellI].yy() = volFfield[cellI].xx();
            volFfield[cellI].zz() = volFfield[cellI].xx();

            wireResistance += volFfield[cellI].xx() * mesh.V()[cellI];
            wireVolume += mesh.V()[cellI];
        }
    }

    Info << "wireResistance is, " << wireResistance << endl;

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

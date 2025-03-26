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
    argList::addOption //
    (
        "DirIndex",
        "labelList",
        "index of Dirichlet boundary patches, e.g., '(0)' or '(0 1)'."
    );
    
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
	// runTime and mesh are instances of objects (or classes).

    volVectorField U_orig // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            mesh.time().path()/"0.orig",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    labelList DirIndex (args.getList<label>("DirIndex"));

    forAll(DirIndex, indexI)
    {
        label patchI_ (DirIndex[indexI]);

        fvPatchField<vector>& inletUorig = U_orig.boundaryFieldRef()[patchI_];

        // loop over all hub faces
        forAll(inletUorig, faceI)
        {
            label adjcellID = mesh.boundary()[patchI_].patch().faceCells()[faceI];

            scalar velz = 50*(0.01 + 2*mesh.C()[adjcellID].x() + mesh.C()[adjcellID].y());
            if (velz >= 0.1 && velz <= 1.0)
            {
                inletUorig[faceI] = vector(0, 0, velz);
            }
            else if (velz > 1.0)
            {
                inletUorig[faceI] = vector(0, 0, 1.0);
            }
            else
            {
                inletUorig[faceI] = vector(0, 0, 0.0);
            }
        }
    }
    
    U_orig.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

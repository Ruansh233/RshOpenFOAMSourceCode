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
#include "IFstream.H"
#include "stringOps.H"


int main(int argc, char *argv[])
{
    argList::addOption //
    (
        "DirIndex",
        "labelList",
        "index of Dirichlet boundary patches, e.g., '(0)' or '(0 1)'."
    );

    argList::addOption //
    (
        "maxValueVector",
        "vector",
        "the maximum value of the distribution, e.g., '(0 0 1)'."
    );

    argList::addOption //
    (
        "order",
        "scalar",
        "order of the profile, e.g., '2'."
    );

    argList::addOption //
    (
        "minRatio",
        "scalar",
        "the minimum value of the distribution, e.g. '0.02'. The default is 0.02."
    );
    
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"
	// runTime and mesh are instances of objects (or classes).

    volScalarField wallDistance_ = wallDist(mesh, "wall").y();

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

    // read the order of the distribution
    scalar order;
    if (args.readIfPresent("order", order))
    {
        // code
    }
    else
    {
        Info<< "\t order is not specified, set it to 2." << endl;
        order = 2;
    }

    // read the minimum value of the distribution
    scalar minRatio;
    if (args.readIfPresent("minRatio", minRatio))
    {
        // code
    }
    else
    {
        Info<< "\t minRatio is not specified, set it to 0.02." << endl;
        minRatio = 0.02;
    }
    
    // read the maximum value of the distribution
    vector maxValueVector (0, 0, 0);
    if (args.found("maxValueVector"))
    {
        maxValueVector = args.get<vector>("maxValueVector");
    }
    else
    {
        FatalErrorIn("nonUniformDirichletBCs")
            << "maxValueVector is not specified." << endl
            << "Please specify it." << endl
            << exit(FatalError);
    }

    forAll(DirIndex, indexI)
    {
        label patchI_ (DirIndex[indexI]);
        fvPatchField<vector>& inletUorig = U_orig.boundaryFieldRef()[patchI_];
        vector averageValue (0, 0, 0);

        // loop over all hub faces
        forAll(mesh.boundary()[patchI_].Cf(), faceI)
        {
            label adjcellID = mesh.boundary()[patchI_].patch().faceCells()[faceI];
            scalar wallDistanceValue = wallDistance_[adjcellID];
            scalar ratioFace = Foam::pow(wallDistanceValue/max(wallDistance_).value(), order);

            if (ratioFace < minRatio)
            {
                inletUorig[faceI] = minRatio * maxValueVector;
                averageValue += inletUorig[faceI] * mesh.boundary()[patchI_].magSf()[faceI];
            }
            else
            {
                inletUorig[faceI] = ratioFace * maxValueVector;
                averageValue += inletUorig[faceI] * mesh.boundary()[patchI_].magSf()[faceI];
            }
        }

        averageValue = averageValue / gSum(mesh.boundary()[patchI_].magSf());
        Info<< "\t Patch "<< patchI_ << " is called: " << mesh.boundary()[patchI_].name() << endl
            << "\t The order of the distribution is: " << order << endl
            << "\t The maximum value of the distribution is: " << maxValueVector 
            << ". The minimum value of the distribution is: " << minRatio << endl
            << "\t The averageValue: " << averageValue << endl;    
    }
    
    U_orig.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

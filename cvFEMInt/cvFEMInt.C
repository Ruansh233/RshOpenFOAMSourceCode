/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    

Group
    

Description
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "writeMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addOption // string variable
    (
        "dict",
        "name",
        "alternative cvFEMIntDict"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    // -----------------------------------------------------------
    // --------- select the reference cell and surfaces ----------
    // -----------------------------------------------------------

    wordRe cvNameRex ("cv[1-2,11]*_Zone");
    cvNameRex.compile();
    const List<label> selectedCVs(mesh.cellZones().indices(cvNameRex));
    Info<< "selectedCVs: " << selectedCVs << endl;

    wordRe cvFaceNameRex ("cv[1-2,11]*_face");
    cvFaceNameRex.compile();
    const List<label> selectedCVFaces(mesh.faceZones().indices(cvFaceNameRex));
    Info<< "selectedCVFaces: " << selectedCVFaces << endl;


    // -----------------------------------------------------------
    // - integral of pressure gradient and save the z direction --
    // -----------------------------------------------------------
    
    List<word> vectorFields ({"grad(pMode1)", "grad(pMode2)"});
    Info<< "vectorFields: " << vectorFields << endl;

    forAll(vectorFields, i)
    {
        volVectorField vectorFieldRead
        (
            IOobject
            (
                vectorFields[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        forAll(selectedCVs, cvI)
        {
            scalar vInt (0);
            
            forAll(mesh.cellZones()[selectedCVs[cvI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCVs[cvI]][cellI];
                vInt += (mesh.V()[cell] * vectorFieldRead[cell]).z();
            }

            Info<< vectorFields[i] 
                << ", cvI, " << selectedCVs[cvI] << ", "
                << vInt << endl;
        }
    }
    
    // -----------------------------------------------------------
    // - integral of pressure div.grad and save the z direction --
    // -----------------------------------------------------------
    
    List<word> scalarFields ({"div(grad(pMode1))", "div(grad(pMode2))"});
    Info<< "scalarFields: " << scalarFields << endl;

    forAll(scalarFields, i)
    {
        volScalarField scalarFieldRead
        (
            IOobject
            (
                scalarFields[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        );
        
        forAll(selectedCVs, cvI)
        {
            scalar vInt (0);
            
            forAll(mesh.cellZones()[selectedCVs[cvI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCVs[cvI]][cellI];
                vInt += mesh.V()[cell] * scalarFieldRead[cell];
            }

            Info<< scalarFields[i]
                << ", cvI, " << selectedCVs[cvI] << ", "
                << vInt << endl;
        }
    }

    // -----------------------------------------------------------
    // ---------- integral of laplacian of velocity  -------------
    // -----------------------------------------------------------

    vectorFields = {"laplacianuMode1", "laplacianuMode2"};
    const scalar mu = 1.0e-6;
    Info<< "vectorFields: " << vectorFields << endl;

    forAll(vectorFields, i)
    {
        volVectorField vectorFieldRead
        (
            IOobject
            (
                vectorFields[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        forAll(selectedCVs, cvI)
        {
            scalar vInt (0);
            
            forAll(mesh.cellZones()[selectedCVs[cvI]], cellI)
            {
                label cell = mesh.cellZones()[selectedCVs[cvI]][cellI];
                vInt += (mesh.V()[cell] * vectorFieldRead[cell]).z();
            }

            Info<< vectorFields[i] 
                << ", cvI, " << selectedCVs[cvI] << ", "
                << -mu*vInt << endl;
        }
    }

    // -----------------------------------------------------------
    // ---------- integral of velocity in inlet patch ------------
    // -----------------------------------------------------------

    const label inletPatchID (mesh.boundaryMesh().findPatchID("block1_in"));
    Info<< "inletPatchID: " << inletPatchID << endl;
    scalar Uin (0.2);

    vectorFields = {"uMode1", "uMode2"};

    forAll(vectorFields, i)
    {
        volVectorField vectorFieldRead
        (
            IOobject
            (
                vectorFields[i],
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        
        scalar sInt (0);
        forAll(mesh.boundary()[inletPatchID], faceI)
        {
            sInt += mesh.boundary()[inletPatchID].Sf()[faceI] 
                    & vectorFieldRead.boundaryField()[inletPatchID][faceI];            
        }
        Info<< vectorFields[i] 
                << ", sInlet_unknown, "
                << sInt << endl;

        sInt = 0;
        forAll(mesh.boundary()[inletPatchID], faceI)
        {
            sInt += mesh.boundary()[inletPatchID].magSf()[faceI];
        }
        Info<< vectorFields[i] 
                << ", sInlet_known, "
                << sInt * Uin << endl;
    }

    Info<< "End" << endl;
}
    
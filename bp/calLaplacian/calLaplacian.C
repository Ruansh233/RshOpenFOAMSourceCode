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
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// void calcLap(GeometricField<scalar, fvPatchField, volMesh> field)
// {
//     volScalarField scalarFieldWrite
//     (
//         IOobject
//         (
//             "Laplacian"+field,
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh,
//         fvc::laplacian(field)
//     );
//     scalarFieldWrite.write();
// }

// void calcLap(GeometricField<vector, fvPatchField, volMesh> field)
// {
//     volVectorField vectorFieldWrite
//     (
//         IOobject
//         (
//             "Laplacian"+field,
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         mesh,
//         fvc::laplacian(field)
//     );
//     vectorFieldWrite.write();
// }


int main(int argc, char *argv[])
{
    argList::addOption // string variable
    (
        "scalarfields",
        "list",
        "scalar fields to be processed"
    );

    argList::addOption // string variable
    (
        "vectorfields",
        "list",
        "vector fields to be processed"
    );
    
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"    
    
    // -----------------------------------------------------------
    // ----------------- read and process scalar fields ----------
    // -----------------------------------------------------------

    List<word> scalarFieldsName;
    List<word> vectorFieldsName;    
    
    
    if (!args.readListIfPresent<word>("scalarfields", scalarFieldsName) 
        && !args.readListIfPresent<word>("vectorfields", vectorFieldsName))
    {
        Info<< "No fields to be processed" << endl;
        return 0;
    }
    

    if (args.readListIfPresent<word>("scalarfields", scalarFieldsName))
    {
        Info<< "scalarFieldsName: " << scalarFieldsName << endl;
        
        forAll(scalarFieldsName, fieldI)
        {
            volScalarField scalarFieldRead
            (
                IOobject
                (
                    scalarFieldsName[fieldI],
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            
            volScalarField scalarFieldWrite
            (
                IOobject
                (
                    "laplacian" + scalarFieldsName[fieldI],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(scalarFieldRead)
            );

            scalarFieldWrite.write();
        }
    }
    
    
    // -----------------------------------------------------------
    // ----------------- read and process vector fields ----------
    // -----------------------------------------------------------

    if (args.readListIfPresent<word>("vectorfields", vectorFieldsName))
    {
        Info<< "vectorFieldsName: " << vectorFieldsName << endl;
        
        forAll(vectorFieldsName, fieldI)
        {
            volVectorField vectorFieldRead
            (
                IOobject
                (
                    vectorFieldsName[fieldI],
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );
            
            volVectorField vectorFieldWrite
            (
                IOobject
                (
                    "laplacian" + vectorFieldsName[fieldI],
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(vectorFieldRead)
            );

            vectorFieldWrite.write();
        }
    }
    

    Info<< "End of program" << endl;

    return 0;
}
    
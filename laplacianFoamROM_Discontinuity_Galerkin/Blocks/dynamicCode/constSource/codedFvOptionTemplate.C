/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

#include "codedFvOptionTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
#include "fvMatrix.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 1270e31d939fb5c5f582428782dce9c90bf3cbb1
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void constSource_1270e31d939fb5c5f582428782dce9c90bf3cbb1(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(constSourceFvOptionscalarSource, 0);
addRemovableToRunTimeSelectionTable
(
    option,
    constSourceFvOptionscalarSource,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constSourceFvOptionscalarSource::
constSourceFvOptionscalarSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    if (false)
    {
        printMessage("Construct constSource fvOption from dictionary");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

constSourceFvOptionscalarSource::
~constSourceFvOptionscalarSource()
{
    if (false)
    {
        printMessage("Destroy constSource");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
constSourceFvOptionscalarSource::correct
(
    GeometricField<scalar, fvPatchField, volMesh>& fld
)
{
    if (false)
    {
        Info<< "constSourceFvOptionscalarSource::correct()\n";
    }

//{{{ begin code
    
//}}} end code
}


void
constSourceFvOptionscalarSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "constSourceFvOptionscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 29 "/home/ruan/Desktop/Shenhui_ruan/OpenFOAM/OpenFOAMsourcecode/laplacianFoamROM_Discontinuity_Galerkin/Blocks/constant/fvOptions.codedSource"
const scalarField& V = mesh_.V();
        scalarField& heSource = eqn.source();

        // Retrieve the x component of the cell centres
        const scalarField& cellx = mesh_.C().component(0);
        const scalarField& celly = mesh_.C().component(1);
        const scalarField& cellz = mesh_.C().component(2);

        // Apply the source
        forAll(mesh_.C(), cellI)
        {
            // // cell volume specific source
            // scalar relx = cellx[cellI]*constant::mathematical::pi;
            // scalar rely = celly[cellI]*constant::mathematical::pi;
            scalar relz = cellz[cellI]/(gMax(cellz)-gMin(cellz))*constant::mathematical::pi;

            // heSource[cellI] = 10.0*(sin(relx) + sin(rely) + sin(relz))*V[cellI] + 100.0*V[cellI];

            heSource[cellI] = 5.0e2*cos(relz)*V[cellI];

        };
//}}} end code
}


void
constSourceFvOptionscalarSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "constSourceFvOptionscalarSource::addSup()\n";
    }

//{{{ begin code
    #line 29 "/home/ruan/Desktop/Shenhui_ruan/OpenFOAM/OpenFOAMsourcecode/laplacianFoamROM_Discontinuity_Galerkin/Blocks/constant/fvOptions.codedSource"
const scalarField& V = mesh_.V();
        scalarField& heSource = eqn.source();

        // Retrieve the x component of the cell centres
        const scalarField& cellx = mesh_.C().component(0);
        const scalarField& celly = mesh_.C().component(1);
        const scalarField& cellz = mesh_.C().component(2);

        // Apply the source
        forAll(mesh_.C(), cellI)
        {
            // // cell volume specific source
            // scalar relx = cellx[cellI]*constant::mathematical::pi;
            // scalar rely = celly[cellI]*constant::mathematical::pi;
            scalar relz = cellz[cellI]/(gMax(cellz)-gMin(cellz))*constant::mathematical::pi;

            // heSource[cellI] = 10.0*(sin(relx) + sin(rely) + sin(relz))*V[cellI] + 100.0*V[cellI];

            heSource[cellI] = 5.0e2*cos(relz)*V[cellI];

        };
//}}} end code
}


void
constSourceFvOptionscalarSource::constrain
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (false)
    {
        Info<< "constSourceFvOptionscalarSource::constrain()\n";
    }

//{{{ begin code
    
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// ************************************************************************* //


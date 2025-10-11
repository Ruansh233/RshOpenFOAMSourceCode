/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#ifndef MapMeshes_H
#define MapMeshes_H

#include "MapVolFields.H"
#include "MapConsistentVolFields.H"
#include "mapLagrangian.H"
#include "UnMapped.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<template<class> class CombineOp>
void MapConsistentMesh
(
    Time& runTimeSource,
    Time& runTimeTarget,
    const instantList& timeSourceDirs,
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const wordRes& selectedFields
)
{
    // Create the interpolation scheme
    meshToMesh0 meshToMesh0Interp(meshSource, meshTarget);

    forAll(timeSourceDirs, timeI)
    {
        const word& sourceTimeName = timeSourceDirs[timeI].name();

        if (sourceTimeName == "0" || sourceTimeName == "constant")
        {
            continue;
        }

        // Set source and target time
        runTimeSource.setTime(timeSourceDirs[timeI], timeI);
        runTimeTarget.setTime(timeSourceDirs[timeI], timeI);

        Info<< nl
            << "Consistently creating and mapping fields for time "
            << meshSource.time().timeName() << nl << endl;

        {
            // Search for list of objects for this time
            IOobjectList objects(meshSource, meshSource.time().timeName());

            // Map volFields
            // ~~~~~~~~~~~~~
            MapConsistentVolFields<scalar>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<scalar>(),
                selectedFields
            );
            MapConsistentVolFields<vector>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<vector>(),
                selectedFields
            );
            MapConsistentVolFields<sphericalTensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<sphericalTensor>(),
                selectedFields
            );
            MapConsistentVolFields<symmTensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<symmTensor>(),
                selectedFields
            );
            MapConsistentVolFields<tensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<tensor>(),
                selectedFields
            );
        }

        {
            // Search for list of target objects for this time
            IOobjectList objects(meshTarget, meshTarget.time().timeName());

            // Mark surfaceFields as unmapped
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            UnMapped<surfaceScalarField>(objects);
            UnMapped<surfaceVectorField>(objects);
            UnMapped<surfaceSphericalTensorField>(objects);
            UnMapped<surfaceSymmTensorField>(objects);
            UnMapped<surfaceTensorField>(objects);

            // Mark pointFields as unmapped
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            UnMapped<pointScalarField>(objects);
            UnMapped<pointVectorField>(objects);
            UnMapped<pointSphericalTensorField>(objects);
            UnMapped<pointSymmTensorField>(objects);
            UnMapped<pointTensorField>(objects);
        }

        mapLagrangian(meshToMesh0Interp);
    }
}


template<template<class> class CombineOp>
void MapSubMesh
(
    Time& runTimeSource,
    Time& runTimeTarget,
    const instantList& timeSourceDirs,
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const meshToMesh0::order& mapOrder,
    const wordRes& selectedFields
)
{
    // Create the interpolation scheme
    meshToMesh0 meshToMesh0Interp
    (
        meshSource,
        meshTarget,
        patchMap,
        cuttingPatches
    );

    forAll(timeSourceDirs, timeI)
    {
        const word& sourceTimeName = timeSourceDirs[timeI].name();

        if (sourceTimeName == "0" || sourceTimeName == "constant")
        {
            continue;
        }

        // Set source and target time
        runTimeSource.setTime(timeSourceDirs[timeI], timeI);
        runTimeTarget.setTime(timeSourceDirs[timeI], timeI);

        Info<< nl
            << "Mapping fields for time " << meshSource.time().timeName()
            << nl << endl;

        {
            // Search for list of source objects for this time
            IOobjectList objects(meshSource, meshSource.time().timeName());

            // Map volFields
            // ~~~~~~~~~~~~~
            MapVolFields<scalar>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<scalar>(),
                selectedFields
            );
            MapVolFields<vector>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<vector>(),
                selectedFields
            );
            MapVolFields<sphericalTensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<sphericalTensor>(),
                selectedFields
            );
            MapVolFields<symmTensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<symmTensor>(),
                selectedFields
            );
            MapVolFields<tensor>
            (
                objects,
                meshToMesh0Interp,
                mapOrder,
                CombineOp<tensor>(),
                selectedFields
            );
        }

        {
            // Search for list of target objects for this time
            IOobjectList objects(meshTarget, meshTarget.time().timeName());

            // Mark surfaceFields as unmapped
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            UnMapped<surfaceScalarField>(objects);
            UnMapped<surfaceVectorField>(objects);
            UnMapped<surfaceSphericalTensorField>(objects);
            UnMapped<surfaceSymmTensorField>(objects);
            UnMapped<surfaceTensorField>(objects);

            // Mark pointFields as unmapped
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            UnMapped<pointScalarField>(objects);
            UnMapped<pointVectorField>(objects);
            UnMapped<pointSphericalTensorField>(objects);
            UnMapped<pointSymmTensorField>(objects);
            UnMapped<pointTensorField>(objects);
        }

        mapLagrangian(meshToMesh0Interp);
    }
}


template<template<class> class CombineOp>
void MapConsistentSubMesh
(
    Time& runTimeSource,
    Time& runTimeTarget,
    const instantList& timeSourceDirs,
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const wordRes& selectedFields
)
{
    HashTable<word> patchMap;
    HashTable<label> cuttingPatchTable;

    forAll(meshTarget.boundary(), patchi)
    {
        if (!isA<processorFvPatch>(meshTarget.boundary()[patchi]))
        {
            patchMap.insert
            (
                meshTarget.boundary()[patchi].name(),
                meshTarget.boundary()[patchi].name()
            );
        }
        else
        {
            cuttingPatchTable.insert
            (
                meshTarget.boundaryMesh()[patchi].name(),
                -1
            );
        }
    }

    MapSubMesh<CombineOp>
    (
        runTimeSource,
        runTimeTarget,
        timeSourceDirs,
        meshSource,
        meshTarget,
        patchMap,
        cuttingPatchTable.toc(),
        mapOrder,
        selectedFields
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

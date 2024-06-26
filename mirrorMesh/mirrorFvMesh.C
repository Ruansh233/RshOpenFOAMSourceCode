/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "mirrorFvMesh.H"
#include "Time.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mirrorFvMesh::mirrorFvMesh(const IOobject& io)
:
    mirrorFvMesh
    (
        io,
        IOdictionary
        (
            IOobject
            (
                "mirrorMeshDict",
                io.time().system(),
                io.time(),
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        )
    )
{}


Foam::mirrorFvMesh::mirrorFvMesh
(
    const IOobject& io,
    const IOdictionary& mirrorDict
)
:
    fvMesh(io)
{
    plane mirrorPlane(mirrorDict);

    const scalar planeTolerance
    (
        mirrorDict.get<scalar>("planeTolerance")
    );

    const pointField& oldPoints = points();
    const faceList& oldFaces = faces();
    const cellList& oldCells = cells();
    const label nOldInternalFaces = nInternalFaces();
    const polyPatchList& oldPatches = boundaryMesh();

    Info<< "Mirroring mesh at origin:" << mirrorPlane.origin()
        << " normal:" << mirrorPlane.normal() << nl;

    // Mirror the points
    Info<< "Mirroring points. Old points: " << oldPoints.size();

    pointField newPoints(2*oldPoints.size());
    label nNewPoints = 0;

    labelList mirrorPointLookup(oldPoints.size(), -1);

    if (!removeOriginal)
    {
        // Grab the old points
        for (const point& pt : oldPoints)
        {
            newPoints[nNewPoints] = pt;
            ++nNewPoints;
        }
    }

    forAll(oldPoints, pointi)
    {
        const scalar alpha =
            mirrorPlane.normalIntersect
            (
                oldPoints[pointi],
                mirrorPlane.normal()
            );

        // Check plane on tolerance
        if (mag(alpha) > planeTolerance)
        {
            // The point gets mirrored
            newPoints[nNewPoints] =
                oldPoints[pointi] + 2.0*alpha*mirrorPlane.normal();

            // remember the point correspondence
            mirrorPointLookup[pointi] = nNewPoints;
            nNewPoints++;
        }
        else
        {
            if (!removeOriginal)
            {
                // The point is on the plane and does not get mirrored
                // Adjust plane location
                newPoints[nNewPoints] =
                    oldPoints[pointi] + alpha*mirrorPlane.normal();

                mirrorPointLookup[pointi] = pointi;
            }
            else
            {
                // As the old points are not grapped, 
                // the point is on the plane is drectly assigned to the new points
                newPoints[nNewPoints] = oldPoints[pointi];

                mirrorPointLookup[pointi] = nNewPoints;
                nNewPoints++;
            }
        }
    }

    // Reset the size of the point list
    Info<< " New points: " << nNewPoints << endl;
    newPoints.setSize(nNewPoints);

    // Construct new to old map
    pointMapPtr_.reset(new labelList(newPoints.size()));
    labelList& pointMap = pointMapPtr_();

    if (!removeOriginal)
    {
        // Insert old points
        forAll(oldPoints, oldPointi)
        {
            pointMap[oldPointi] = oldPointi;
        }
    }
    forAll(mirrorPointLookup, oldPointi)
    {
        pointMap[mirrorPointLookup[oldPointi]] = oldPointi;
    }


    Info<< "Mirroring faces. Old faces: " << oldFaces.size();

    // Algorithm:
    // During mirroring, the faces that were previously boundary faces
    // in the mirror plane may become internal faces. In order to
    // deal with the ordering of the faces, the algorithm is split
    // into two parts.  For original faces, the internal faces are
    // distributed to their owner cells.  Once all internal faces are
    // distributed, the boundary faces are visited and if they are in
    // the mirror plane they are added to the master cells (the future
    // boundary faces are not touched).  After the first phase, the
    // internal faces are collected in the cell order and numbering
    // information is added.  Then, the internal faces are mirrored
    // and the face numbering data is stored for the mirrored section.
    // Once all the internal faces are mirrored, the boundary faces
    // are added by mirroring the faces patch by patch.

    // Distribute internal faces
    labelListList newCellFaces(oldCells.size());

    const labelUList& oldOwnerStart = lduAddr().ownerStartAddr();

    if (!removeOriginal)
    {
        forAll(newCellFaces, celli)
        {
            labelList& curFaces = newCellFaces[celli];

            const label s = oldOwnerStart[celli];
            const label e = oldOwnerStart[celli + 1];

            curFaces.setSize(e - s);

            forAll(curFaces, i)
            {
                curFaces[i] = s + i;
            }
        }
    }

    // Distribute boundary faces.  Remember the faces that have been inserted
    // as internal
    boolListList insertedBouFace(oldPatches.size());

    forAll(oldPatches, patchi)
    {
        const polyPatch& curPatch = oldPatches[patchi];

        if (curPatch.coupled())
        {
            WarningInFunction
                << "Found coupled patch " << curPatch.name() << endl
                << "    Mirroring faces on coupled patches destroys"
                << " the ordering. This might be fixed by running a dummy"
                << " createPatch afterwards." << endl;
        }

        boolList& curInsBouFace = insertedBouFace[patchi];

        curInsBouFace.setSize(curPatch.size());
        curInsBouFace = false;

        // Get faceCells for face insertion
        const labelUList& curFaceCells = curPatch.faceCells();
        const label curStart = curPatch.start();

        forAll(curPatch, facei)
        {
            // Find out if the mirrored face is identical to the
            // original.  If so, the face needs to become internal and
            // added to its owner cell
            const face& origFace = curPatch[facei];

            face mirrorFace(origFace.size());
            forAll(mirrorFace, pointi)
            {
                mirrorFace[pointi] = mirrorPointLookup[origFace[pointi]];
            }

            // No internal faces are added if the original cells are
            // removed
            if (!removeOriginal)
            {
                if (origFace == mirrorFace)
                {
                    // The mirror is identical to current face.  This will
                    // become an internal face
                    const label oldSize = newCellFaces[curFaceCells[facei]].size();

                    newCellFaces[curFaceCells[facei]].setSize(oldSize + 1);
                    newCellFaces[curFaceCells[facei]][oldSize] = curStart + facei;

                    curInsBouFace[facei] = true;
                }
            }
        }
    }

    // Construct the new list of faces.  Boundary faces are added
    // last, cush that each patch is mirrored separately.  The
    // addressing is stored in two separate arrays: first for the
    // original cells (face order has changed) and then for the
    // mirrored cells.
    labelList masterFaceLookup(oldFaces.size(), -1);
    labelList mirrorFaceLookup(oldFaces.size(), -1);

    faceList newFaces(2*oldFaces.size());
    label nNewFaces = 0;

    if (!removeOriginal)
    {
        // Insert original (internal) faces
        forAll(newCellFaces, celli)
        {
            const labelList& curCellFaces = newCellFaces[celli];

            forAll(curCellFaces, cfI)
            {
                newFaces[nNewFaces] = oldFaces[curCellFaces[cfI]];
                masterFaceLookup[curCellFaces[cfI]] = nNewFaces;

                nNewFaces++;
            }
        }
    }
    
    // Mirror internal faces
    for (label facei = 0; facei < nOldInternalFaces; facei++)
    {
        const face& oldFace = oldFaces[facei];
        face& nf = newFaces[nNewFaces];
        nf.setSize(oldFace.size());

        nf[0] = mirrorPointLookup[oldFace[0]];

        for (label i = 1; i < oldFace.size(); i++)
        {
            nf[i] = mirrorPointLookup[oldFace[oldFace.size() - i]];
        }

        mirrorFaceLookup[facei] = nNewFaces;
        nNewFaces++;
    }

    // Mirror boundary faces patch by patch
    labelList newPatchSizes(boundary().size(), Zero);
    labelList newPatchStarts(boundary().size(), -1);

    forAll(boundaryMesh(), patchi)
    {
        const label curPatchSize = boundaryMesh()[patchi].size();
        const label curPatchStart = boundaryMesh()[patchi].start();
        const boolList& curInserted = insertedBouFace[patchi];

        newPatchStarts[patchi] = nNewFaces;

        if (!removeOriginal)
        {
            // Master side
            for (label facei = 0; facei < curPatchSize; facei++)
            {
                // Check if the face has already been added.  If not, add it and
                // insert the numbering details.
                if (!curInserted[facei])
                {
                    newFaces[nNewFaces] = oldFaces[curPatchStart + facei];

                    masterFaceLookup[curPatchStart + facei] = nNewFaces;
                    nNewFaces++;
                }
            }
        }

        // Mirror side
        for (label facei = 0; facei < curPatchSize; facei++)
        {
            // Check if the face has already been added.  If not, add it and
            // insert the numbering details.
            if (!curInserted[facei])
            {
                const face& oldFace = oldFaces[curPatchStart + facei];
                face& nf = newFaces[nNewFaces];
                nf.setSize(oldFace.size());

                nf[0] = mirrorPointLookup[oldFace[0]];

                for (label i = 1; i < oldFace.size(); i++)
                {
                    nf[i] = mirrorPointLookup[oldFace[oldFace.size() - i]];
                }

                mirrorFaceLookup[curPatchStart + facei] = nNewFaces;
                nNewFaces++;
            }
            else
            {
                if (!removeOriginal)
                {
                    // Grab the index of the master face for the mirror side
                    mirrorFaceLookup[curPatchStart + facei] =
                        masterFaceLookup[curPatchStart + facei];
                }
            }
        }

        // If patch exists, grab the name and type of the original patch
        if (nNewFaces > newPatchStarts[patchi])
        {
            newPatchSizes[patchi] =
                nNewFaces - newPatchStarts[patchi];
        }
    }

    // Tidy up the lists
    newFaces.setSize(nNewFaces);
    Info<< " New faces: " << nNewFaces << endl;

    Info<< "Mirroring patches. Old patches: " << boundary().size()
        << " New patches: " << boundary().size() << endl;

    if (!removeOriginal)
    {
        Info<< "Mirroring cells. Old cells: " << oldCells.size()
            << " New cells: " << 2*oldCells.size() << endl;
    }
    else
    {
        Info<< "Mirroring cells. Old cells: " << oldCells.size()
            << " New cells: " << oldCells.size() << endl;
    }

    cellList newCells(2*oldCells.size());
    label nNewCells = 0;

    // Construct new to old cell map
    cellMapPtr_.reset(new labelList(newCells.size()));
    labelList& cellMap = cellMapPtr_();

    if (!removeOriginal)
    {
        // Grab the original cells.  Take care of face renumbering.
        forAll(oldCells, celli)
        {
            const cell& oc = oldCells[celli];

            cell& nc = newCells[nNewCells];
            nc.setSize(oc.size());

            forAll(oc, i)
            {
                nc[i] = masterFaceLookup[oc[i]];
            }

            cellMap[nNewCells] = celli;

            nNewCells++;
        }
    }

    // Mirror the cells
    forAll(oldCells, celli)
    {
        const cell& oc = oldCells[celli];

        cell& nc = newCells[nNewCells];
        nc.setSize(oc.size());

        forAll(oc, i)
        {
            nc[i] = mirrorFaceLookup[oc[i]];
        }

        cellMap[nNewCells] = celli;

        nNewCells++;
    }

    if (removeOriginal)
    {
        // Reset the size of the cell list
        newCells.setSize(nNewCells);
    }

    // Mirror the cell shapes
    Info<< "Mirroring cell shapes." << endl;

    Info<< nl << "Creating new mesh" << endl;
    mirrorMeshPtr_ = autoPtr<fvMesh>::New
    (
        io,
        std::move(newPoints),
        std::move(newFaces),
        std::move(newCells)
    );

    fvMesh& pMesh = mirrorMeshPtr_();

    // Add the boundary patches
    polyPatchList newPatches(newPatchSizes.size());

    forAll(newPatches, patchi)
    {
        newPatches.set
        (
            patchi,
            boundaryMesh()[patchi].clone
            (
                pMesh.boundaryMesh(),
                patchi,
                newPatchSizes[patchi],
                newPatchStarts[patchi]
            )
        );
    }

    pMesh.addPatches(newPatches);
}


// ************************************************************************* //

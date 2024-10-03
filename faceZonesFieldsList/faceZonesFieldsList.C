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
#include "Time.H"
#include "timeSelector.H"
#include "IFstream.H"

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "fields",
        "wordList",
        "Specify the fields to be exported."
    );

    argList::addOption
    (
        "faceZones",
        "wordList",
        "Specify the fields to be exported."
    );

    timeSelector::addOptions(false, true);
    
    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"   

    const fileName dataPath ("boundaryFields");
    if(!isDir(runTime.globalPath()/dataPath))
        mkDir(runTime.globalPath()/dataPath);

    const wordList fieldsName (args.getList<word>("fields"));
    const wordList faceZonesName (args.getList<word>("faceZones"));

    instantList timeDirs = timeSelector::select0(runTime, args);

    Info<< "Read fields, " << fieldsName << ", in faceZones, "
        << faceZonesName << " for times, " << timeDirs << "." << nl;

    List<label> faceZonesID (faceZonesName.size());

    forAll(faceZonesName, i)
    {
        faceZonesID[i] = mesh.faceZones().findZoneID(faceZonesName[i]);
    }
    
    if(fieldsName.found("U"))
    {
        PtrList<List<vector>> fieldsList;

        forAll(timeDirs, timei)
        {
            volVectorField U // note that velocity is a vector field
            (
                IOobject
                (
                    "U",
                    timeDirs[timei].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            auto surfaceField =
                autoPtr<surfaceVectorField>::New
                (fvc::interpolate(U));

            forAll(faceZonesID, i)
            {
                List<vector> boundary_field (mesh.faceZones()[faceZonesID[i]].size());

                forAll(mesh.faceZones()[faceZonesID[i]], j)
                {
                    label faceID = mesh.faceZones()[faceZonesID[i]][j];
                    if (faceID < mesh.boundary()[0].start())
                    {
                        boundary_field[j] = surfaceField()[faceID];
                    }
                    else
                    {
                        forAll(mesh.boundary(), k)
                        {
                            label respectID = faceID - mesh.boundary()[k].start();
                            if(0 <= respectID && respectID < mesh.boundary()[k].size())
                            {
                                boundary_field[j] = U.boundaryField()[k][faceID - mesh.boundary()[k].start()];
                            }
                        }
                    }
                }

                fieldsList.append(boundary_field.clone());  
            }

            autoPtr<OFstream> outputFilePtr;
            outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"U_boundary"));

            outputFilePtr() << fieldsList; 
        }             
    }

    if(fieldsName.found("gradU"))
    {
        PtrList<List<tensor>> fieldsList;

        forAll(timeDirs, timei)
        {
            volTensorField gradU // note that velocity gradient is a tensor field
            (
                IOobject
                (
                    "gradU",
                    timeDirs[timei].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            auto surfaceField =
                autoPtr<surfaceTensorField>::New
                (fvc::interpolate(gradU));

            forAll(faceZonesID, i)
            {
                List<tensor> boundary_field (mesh.faceZones()[faceZonesID[i]].size());

                forAll(mesh.faceZones()[faceZonesID[i]], j)
                {
                    label faceID = mesh.faceZones()[faceZonesID[i]][j];
                    if (faceID < mesh.boundary()[0].start())
                    {
                        boundary_field[j] = surfaceField()[faceID];
                    }
                    else
                    {
                        forAll(mesh.boundary(), k)
                        {
                            label respectID = faceID - mesh.boundary()[k].start();
                            if(0 <= respectID && respectID < mesh.boundary()[k].size())
                            {
                                boundary_field[j] = gradU.boundaryField()[k][faceID - mesh.boundary()[k].start()];
                            }
                        }
                    }
                }

                fieldsList.append(boundary_field.clone());  
            }

            autoPtr<OFstream> outputFilePtr;
            outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"gradU_boundary"));

            outputFilePtr() << fieldsList; 
        }             
    }

    if(fieldsName.found("p"))
    {
        PtrList<List<scalar>> fieldsList;

        forAll(timeDirs, timei)
        {
            volScalarField p // note that velocity is a vector field
            (
                IOobject
                (
                    "p",
                    timeDirs[timei].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            auto surfaceField =
                autoPtr<surfaceScalarField>::New
                (fvc::interpolate(p));

            forAll(faceZonesID, i)
            {
                List<scalar> boundary_field (mesh.faceZones()[faceZonesID[i]].size());

                forAll(mesh.faceZones()[faceZonesID[i]], j)
                {
                    label faceID = mesh.faceZones()[faceZonesID[i]][j];
                    if (faceID < mesh.boundary()[0].start())
                    {
                        boundary_field[j] = surfaceField()[faceID];
                    }
                    else
                    {
                        forAll(mesh.boundary(), k)
                        {
                            label respectID = faceID - mesh.boundary()[k].start();
                            if(0 <= respectID && respectID < mesh.boundary()[k].size())
                            {
                                boundary_field[j] = p.boundaryField()[k][faceID - mesh.boundary()[k].start()];
                            }
                        }
                    }
                }

                fieldsList.append(boundary_field.clone());  
            }

            autoPtr<OFstream> outputFilePtr;
            outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"p_boundary"));

            outputFilePtr() << fieldsList; 
        }             
    }

    if(fieldsName.found("gradp"))
    {
        PtrList<List<vector>> fieldsList;

        forAll(timeDirs, timei)
        {
            volVectorField gradp // note that velocity is a vector field
            (
                IOobject
                (
                    "gradp",
                    timeDirs[timei].name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                mesh
            );

            auto surfaceField =
                autoPtr<surfaceVectorField>::New
                (fvc::interpolate(gradp));

            forAll(faceZonesID, i)
            {
                List<vector> boundary_field (mesh.faceZones()[faceZonesID[i]].size());

                forAll(mesh.faceZones()[faceZonesID[i]], j)
                {
                    label faceID = mesh.faceZones()[faceZonesID[i]][j];
                    if (faceID < mesh.boundary()[0].start())
                    {
                        boundary_field[j] = surfaceField()[faceID];
                    }
                    else
                    {
                        forAll(mesh.boundary(), k)
                        {
                            label respectID = faceID - mesh.boundary()[k].start();
                            if(0 <= respectID && respectID < mesh.boundary()[k].size())
                            {
                                boundary_field[j] = gradp.boundaryField()[k][faceID - mesh.boundary()[k].start()];
                            }
                        }
                    }
                }

                fieldsList.append(boundary_field.clone());  
            }

            autoPtr<OFstream> outputFilePtr;
            outputFilePtr.reset(new OFstream(runTime.globalPath()/dataPath/"gradp_boundary"));

            outputFilePtr() << fieldsList; 
        }             
    }
}


// ************************************************************************* //

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
#include "regionProperties.H"

int main(int argc, char *argv[])
{   
    argList::addOption
    (
        "pointsFile",
        "word",
        "Specify the relative directory of the points file. (e.g. constant/polyMesh/points)"
    );

    #include "addAllRegionOptions.H"

    // Initialise OF case
    #include "setRootCase.H"

    // These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // Handle -allRegions, -regions, -region
    #include "getAllRegionOptions.H"

    Info<< "Start\n" << endl;

    // read the file path
    fileName pointsFile("points");
    const bool fileFound (args.readIfPresent <fileName> ("pointsFile", pointsFile));
    
    // Extract the directory and file name
    fileName pathOnly = pointsFile.path();
    fileName fileNameOnly = pointsFile.name();

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir =
        (
            regionName == polyMesh::defaultRegion ? word::null : regionName
        );
        const fileName meshDir = regionDir/polyMesh::meshSubDir;

        if (regionNames.size() > 1)
        {
            Info<< "region=" << regionName << nl;
        }

        points.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    plusPostRANS

Description

    calculates y+ and u+ fields for wall-bounded flows computed with
    one of the available low-Re RANS (no wall function!) turbulence 
    models. More specifically it

    :: 	calculates and outputs y+ (avg., min., max.) based on the 
	velocity gradient at the wall  

    ::	calculates and outputs the wall averaged friction velocity 

    ::  writes fields of y+ and U+ to the corresponding time directory

\*---------------------------------------------------------------------------*/
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
// #include "RASModel.H"
// #include "LESModel.H"
#include "nearWallDist.H"
#include "wallDist.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    Foam::argList::addOption
    (
        "tolerance",
        "matchTolerance",
        "specify a matching tolerance fraction of cell-to-face distance and y (default 1.0e-6)"
    );

    argList::addBoolOption
    (
        "noWrite",
        "do not write the y+, u+, and uTau fields and only calculate the min, max, and average. "
    );

    #include "setRootCase.H"

    scalar matchTol;

    bool noWrite = args.found("noWrite");

    if (!args.readIfPresent("tolerance", matchTol))
    {
        matchTol = 1.0e-6;
        Info<< "Using default cell-to-face distance to y match tolerance fraction of " << matchTol << nl << endl;
    }
    else
    {
        Info<< "Using cell distance to y match tolerance fraction of " << matchTol << nl << endl;    
    }


    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        wallDist y(mesh);

	    #include "createFields.H"

        forAll(mesh.boundary(), patchi)
        {
            if (isA<wallFvPatch>(mesh.boundary()[patchi]))
            {
                uTau.boundaryFieldRef()[patchi] = sqrt(turbulence->nu()() * mag(U.boundaryField()[patchi].snGrad()));
            }
        }

        //go through all the cells in the mesh
        forAll(uTau, cellI)
        {
            //loop over all the patches
            forAll(mesh.boundary()[patchi], patchi)
            {            
                //if this patch is a wall...
                if(isA<wallFvPatch>(mesh.boundary()[patchi]))
                {
                    forAll(mesh.boundary()[patchi], facei)
                    {
                        //calculated distance from the current cell to a face on a wall patch
                        scalar cellFaceDist ;

                        cellFaceDist = mag(mesh.C()[cellI] - mesh.Cf().boundaryField()[patchi][facei]);

                        //if the fraction difference is less than or equal to the match tolerance, search no further.
                        if( abs(cellFaceDist - y.y()[cellI]) <= matchTol)
                        { 
                            uTau[cellI] = uTau.boundaryField()[patchi][facei];	
                            break;
                        }
                    }//end for loop over faces 
                }//end if statement checking isA wall 
            }//end loop over patches
        }//end loop over uTau cells

        //true y+ for arbitrary geometry
        yPlus = y.y() * uTau / turbulence->nu();

        //dummy variable holding velocity units
        dimensionedScalar UUnits ( "UUnits", dimensionSet(0,1,-1,0,0), 1.0 ); 

        //true uPlus over arbitrary geometry
        uPlus = U / stabilise(uTau, SMALL*UUnits);//used to fix divide by zero error if uTau is zero
            
        if(noWrite)
        {
            Info << "  noWrite option selected, nothing written." << nl <<endl;
        }
        else
        {
            Info << "  Writing yPlus and uPlus to corresponding fields." << nl <<endl;
            yPlus.write();
            uPlus.write();
            uTau.write();
        }

    }//end loop over time directories

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //

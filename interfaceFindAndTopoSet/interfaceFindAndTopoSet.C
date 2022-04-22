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
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // Rsh, test code to find the center of boundary faces
    // label patchI(0);
    List<vector> patchCenterList(mesh.boundaryMesh().size());
    
    // Rsh, find the center of each patch, and add into patchCenterList
    forAll(mesh.boundaryMesh(), patchI)
    {
        vector patchCenterSum = vector::zero;
        vector patchCenterCoordinate = vector::zero;
        forAll(mesh.boundaryMesh()[patchI], patchFaceI)
        {
            patchCenterSum = patchCenterSum + mesh.boundary()[patchI].Cf()[patchFaceI];
        }
        patchCenterCoordinate = patchCenterSum / mesh.boundary()[patchI].Cf().size(); 
        patchCenterList[patchI] = patchCenterCoordinate;
    }
    
    List<word> interfacePatchListA(mesh.boundaryMesh().size());
    List<word> interfacePatchListB(mesh.boundaryMesh().size());
    label interfacePairNumber(0);

    scalar toleranceOfDistance(1.0e-10);
    scalar distanceOfPatchCenter(0.0);
    List<label> patchCenterChoosed;
    bool choosedFlag;

    // Rsh, find matched patch center, and add them into interfacePatchListA and interfacePatchListB respectively 
    forAll(patchCenterList, patchCenterFirstI)
    {
        vector firstPatchCenter (patchCenterList[patchCenterFirstI]);
      
        // Rsh, avoid repeating selection of patches which is already in pair
        choosedFlag = false;
        forAll(patchCenterChoosed, choosedI)
        {
            if (patchCenterFirstI == patchCenterChoosed[choosedI])
                choosedFlag = true;
        }

        if (choosedFlag)
            continue;

        forAll(patchCenterList, patchCenterSecondI)
        {         
            if (patchCenterSecondI == patchCenterFirstI)
                continue;

            else
            {
                vector secondPatchCenter (patchCenterList[patchCenterSecondI]);
                distanceOfPatchCenter = mag(firstPatchCenter - secondPatchCenter);

                if (distanceOfPatchCenter < toleranceOfDistance)
                {
                    patchCenterChoosed.append(patchCenterSecondI);
                    interfacePatchListA[interfacePairNumber] = mesh.boundary()[patchCenterFirstI].name();
                    interfacePatchListB[interfacePairNumber] = mesh.boundary()[patchCenterSecondI].name();

                    interfacePairNumber += 1;                  
                    continue;                   
                }

            }
                              
        }
        
    }

    interfacePatchListA.setSize(interfacePairNumber);
    interfacePatchListB.setSize(interfacePairNumber);


    Info << "end\n" << endl;

}

// ************************************************************************* //

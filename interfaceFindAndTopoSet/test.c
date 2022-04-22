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
#include "namedDictionary.H"
#include "dictionaryEntry.H"
#include "dictionaryListEntry.H"
#include "IFstream.H"
#include "keyType.H"


int main(int argc, char *argv[])
{
    #include "setRootCase.H"

	// These two create the time system (instance called runTime) and fvMesh (instance called mesh).
    #include "createTime.H"
    #include "createMesh.H"

    // IO dictionary
    IOdictionary vectorDict
    (
        IOobject
        (
            "vectorDict",
            runTime.system(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    // IOdictionary topoSetDict
    // (
    //     IOobject
    //     (
    //         "topoSetDict",
    //         runTime.system(),
    //         runTime,
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     )
    // );

    List<dictionary> testDictList;

    dictionary testDict;

    // List<namedDictionary> testDictList;

    testDict.add("test1", "test1");
    testDict.add("test2", "test2");    

    // testDict.lookupEntry("test1", keyType::option::LITERAL);
    autoPtr<entry> testEntry;

    testEntry.reset(new entry testEntry(testDict.lookupEntry("test1", keyType::option::LITERAL)));



    // entry testEntry(testDict.lookupEntry("test1", keyType::option::LITERAL));

    // Info << testDict.lookupEntry("test1", keyType::option::LITERAL).keyword() << endl;

    // Info << testEntry << endl;
    

    vectorDict.subDictOrAdd("test3");
    // vectorDict.subDict("actions").add("name1", "Ruan");

    // List<namedDictionary> actionEntries(topoSetDict.lookup("actions"));

    // Info << "name: " << actionEntries << endl;

    // PtrList<dictionary> actionEntries(topoSetDict.lookup("actions"));

    // word testword(topoSetDict.subDict("actions").lookup("name"));

    // word testword(topoSetDict.lookup("test"));
    // Info << "name: " << testword << endl;

    // Info << endl << actionEntries << endl;

    // vectorDict.add(testEn);

    vectorDict.subDict("test3").add("name1", "Ruan");
    vectorDict.subDict("test3").add("name2", testDict);

    // Info << testDict << endl;
    // Info << testDict.Keyword() <<endl;

    // vectorDict.add("test", vectorDict);

    // testDictList.append(testDict);
    // testDictList.append(testDict);
    // testDictList.append(testDict);
    // testDictList.append(testDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);
    // testDictList.append(vectorDict);


    // vectorDict.add("Dict1", testDictList);

    // Info << testDictList << endl;

    // keyType typeName("testE");
    // entry testE(typename);

    // fileName dataPath (mesh.time().system());
    // fileName dataFile (dataPath/"topoSetDict");
    // IFstream dataStream(dataFile);

    // dictionaryListEntry dictE(vectorDict, dataStream);

    // Info << dictE << endl;
    // Info << dictE.startLineNumber() << endl;

    // dictE.readEntry();

    // fileName outputDir = mesh.time().constant();
    // // Info << mesh.time().path() << nl;
    // // // Creathe the directory
    // mkDir(outputDir);

    // // // File pointer to direct the output to
	// autoPtr<OFstream> outputFilePtr;
    // // // Open the file in the newly created directory
    // outputFilePtr.reset(new OFstream(outputDir/"vectorDict2", IOstreamOption(), true));

    // outputFilePtr().beginBlock ("hello");
    // outputFilePtr().endBlock();
    // outputFilePtr().endEntry();


    // outputFilePtr() << testDictList;

    // vectorDict.subDict("Dict1").add("testlist", vectorDict);

    vectorDict.regIOobject::write();

    // vectorDict.add("test4", testDict);

//     // Rsh, test code to find the center of boundary faces
//     // label patchI(0);
//     List<vector> patchCenterList(mesh.boundaryMesh().size());

//     // Info << mesh.boundaryMesh().size();
//     // label numberOfPatch(0);
    
//     forAll(mesh.boundaryMesh(), patchI)
//     {
//         vector patchCenterSum = vector::zero;
//         vector patchCenterCoordinate = vector::zero;
//         forAll(mesh.boundaryMesh()[patchI], patchFaceI)
//         {
//             // Info << "The center of Face " << patchFaceI << " in patch " << patchI << ": "
//             //      << ": " << mesh.boundary()[patchI].name() 
//             //      << " is " << mesh.boundary()[patchI].Cf()[patchFaceI]
//             //      << endl;
//             patchCenterSum = patchCenterSum + mesh.boundary()[patchI].Cf()[patchFaceI];
//             // boundaryFacesCenterList[patchFaceI] = mesh.boundary()[patchI].Cf()[patchFaceI];
//         }
//         patchCenterCoordinate = patchCenterSum / mesh.boundary()[patchI].Cf().size(); 
//         // Info << "patchCenterCoordinate.y() " << patchCenterCoordinate.y() << endl;
//         // Info << endl << "The center of patch " << mesh.boundary()[patchI].name() << " is "  
//         //      << "boundaryFacesCenterList[patchFaceI]: " 
//         //      << patchCenterCoordinate
//         //      << endl;
//         patchCenterList[patchI] = patchCenterCoordinate;
//         // numberOfPatch = numberOfPatch + 1;
//     }
//     // Info << "mesh.boundaryMesh().size() " << mesh.boundaryMesh().size() << endl;
//     // Info << patchCenterList;
//     // Info << "numberOfPatch " << numberOfPatch;
//     Info << endl;

    
//     List<word> interfacePatchListA(mesh.boundaryMesh().size());
//     List<word> interfacePatchListB(mesh.boundaryMesh().size());
//     label interfacePairNumber(0);

//     scalar toleranceOfDistance(1.0e-10);
//     scalar distanceOfPatchCenter(0.0);
//     List<label> patchCenterChoosed;
//     bool choosedFlag;
//     // label patchCenterLast(-1);

//     vector firstPatchCenter = vector::zero;
//     vector secondPatchCenter = vector::zero;

//     forAll(patchCenterList, patchCenterFirstI)
//     {
//         firstPatchCenter = patchCenterList[patchCenterFirstI];
      
//         choosedFlag = false;
//         forAll(patchCenterChoosed, choosedI)
//         {
//             if (patchCenterFirstI == patchCenterChoosed[choosedI])
//                 choosedFlag = true;
//         }

//         if (choosedFlag)
//             continue;
        
//         // if (patchCenterFirstI > 7)
//         //     break;

//         // Info << "firstPatchCenter " << firstPatchCenter << endl;
//         forAll(patchCenterList, patchCenterSecondI)
//         {
//             // if (patchCenterSecondI > 10)
//             // break;
            
//             if (patchCenterSecondI == patchCenterFirstI)
//                 continue;
//             else
//             {
//                 // patchCenterSecondI = patchCenterSecondI + 1;
//                 secondPatchCenter = patchCenterList[patchCenterSecondI];
//                 // Info << "secondPatchCenter " << secondPatchCenter << endl;
//                 // Info << firstPatchCenter - secondPatchCenter << endl;
//                 // Info << mag(firstPatchCenter - secondPatchCenter) << endl;
//                 distanceOfPatchCenter = mag(firstPatchCenter - secondPatchCenter);
//                 // Info << "distanceOfPatchCenter " << mag(firstPatchCenter - secondPatchCenter) << endl;
//                 // Info << "distanceOfPatchCenter " << distanceOfPatchCenter << endl;

//                 // Info << "patch " << mesh.boundary()[patchCenterFirstI].name() 
//                 //      << " and " << mesh.boundary()[patchCenterSecondI].name() 
//                 //      << " have center: " << distanceOfPatchCenter << " m "
//                 //      << endl;

//                 if (distanceOfPatchCenter<toleranceOfDistance)
//                 {
//                     // Info << "patch " << patchCenterFirstI << ": " << mesh.boundary()[patchCenterFirstI].name() 
//                     //      << " and " << patchCenterSecondI << ": " << mesh.boundary()[patchCenterSecondI].name() 
//                     //      << " have same center. "
//                     //      << endl;

//                     patchCenterChoosed.append(patchCenterSecondI);
//                     // writeFile << mesh.boundary()[patchCenterFirstI].name()
//                     //           << "     "
//                     //           << mesh.boundary()[patchCenterSecondI].name()
//                     //           << endl;
//                     interfacePatchListA[interfacePairNumber] = mesh.boundary()[patchCenterFirstI].name();
//                     interfacePatchListB[interfacePairNumber] = mesh.boundary()[patchCenterSecondI].name();

//                     interfacePairNumber = interfacePairNumber + 1;
//                     // patchCenterLast = patchCenterSecondI;
                    
//                     continue;
                    
//                 }

//             }
                              
//         }
        
//     }

//     // Info << patchCenterChoosed;

//     interfacePatchListA.setSize(interfacePairNumber);
//     interfacePatchListB.setSize(interfacePairNumber);

//     // Info << "interfacePatchListA" << interfacePatchListA << endl
//     //      << "interfacePatchListB" << interfacePatchListB << endl;

//     // fileName outputFile("patchPair");
//     // std::ofstream writeFile(runTime.path()/runTime.system()/outputFile);
//     // forAll(interfacePatchListA, patchListI)
//     // {
//     //     writeFile << interfacePatchListA[patchListI] << "    " << interfacePatchListB[patchListI] << nl;
//     // }
    
    
//     fileName outputFile("createPatchDict_interfaces");
//     std::ofstream writeFile(runTime.path()/runTime.system()/outputFile);

//     // Rsh, file head
//     const word fileHead = R"deli(
// /*--------------------------------*- C++ -*----------------------------------*\
//   =========                 |
//   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
//    \\    /   O peration     | Website:  https://openfoam.org
//     \\  /    A nd           | Version:  dev
//      \\/     M anipulation  |
// \*---------------------------------------------------------------------------*/
// FoamFile
// {
//     format      ascii;
//     class       dictionary;
//     object      createPatchDict;
// }

// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// pointSync false;

// patches
// (
//     )deli";

//     writeFile << fileHead;

//     // Rsh, interface pair
//     const word emptyDict = R"deli(
//     {
//         name                interface_pair1;
//         patchInfo
//         {
//             type            cyclic;
//             neighbourPatch  interface_pair2;
//             transform       coincidentFullMatch;
//         }
//         constructFrom patches;
//         patches (pair1);
//     }

//     {
//         name                interface_pair2;
//         patchInfo
//         {
//             type            cyclic;
//             neighbourPatch  interface_pair1;
//             transform       coincidentFullMatch;
//         }
//         constructFrom patches;
//         patches (pair2);
//     })deli";

//     forAll(interfacePatchListA, patchListI)
//     {
//         word interfacePairDict(emptyDict);
//         interfacePairDict = replaceStringPart(interfacePairDict, "pair1", interfacePatchListA[patchListI]);
//         interfacePairDict = replaceStringPart(interfacePairDict, "pair2", interfacePatchListB[patchListI]);

//         writeFile << interfacePairDict << nl;
//     }

//     const word endOfFile = R"deli(
// );

// // ************************************************************************* //
//     )deli";

//     writeFile << endOfFile;
//     writeFile.close();

//     // system("cat system/createPatchDict");
//     // system("createPatch");

    return 0;

}

// ************************************************************************* //

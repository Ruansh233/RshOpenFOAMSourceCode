/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "addToRunTimeSelectionTable.H"
#include "nonUniformDarcyForchheimer.H"
#include "geometricOneField.H"
#include "fvMatrices.H"
#include "pointIndList.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace porosityModels
    {
        defineTypeNameAndDebug(nonUniformDarcyForchheimer, 0);
        addToRunTimeSelectionTable(porosityModel, nonUniformDarcyForchheimer, mesh);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityModels::nonUniformDarcyForchheimer::nonUniformDarcyForchheimer
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& cellZoneName
)
:
    porosityModel(name, modelType, mesh, dict, cellZoneName),
    dXYZ_("d", dimless/sqr(dimLength), coeffs_),
    fXYZ_("f", dimless/dimLength, coeffs_),
    D_(cellZoneIDs_.size()),
    F_(cellZoneIDs_.size()),
    Dfield_(mesh_.C().size(), Zero),
    Ffield_(mesh_.C().size(), Zero),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    muName_(coeffs_.getOrDefault<word>("mu", "thermo:mu")),
    nuName_(coeffs_.getOrDefault<word>("nu", "nu"))
{
    // Rsh, 2023-04-18
    volTensorField volDfield
    (
        IOobject
        (
            "nonUniDcoeff",
            mesh.time().timeName(),
            // mesh.time().caseConstant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );
    volTensorField volFfield
    (
        IOobject
        (
            "nonUniFcoeff",
            mesh.time().timeName(),
            // mesh.time().caseConstant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    );   

    forAll(Dfield_, celli)
    {
        Dfield_[celli] = volDfield[celli];
        Ffield_[celli] = volFfield[celli];
    }
    // Rsh, 2023-04-18
    
    adjustNegativeResistance(dXYZ_);
    adjustNegativeResistance(fXYZ_);

    calcTransformModelData();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porosityModels::nonUniformDarcyForchheimer::calcTransformModelData()
{
    // // Rsh, 2023-04-18
    // // The Darcy coefficient as a tensor
    // tensor darcyCoeff(Zero);
    // darcyCoeff.xx() = dXYZ_.value().x();
    // darcyCoeff.yy() = dXYZ_.value().y();
    // darcyCoeff.zz() = dXYZ_.value().z();

    // // The Forchheimer coefficient as a tensor
    // // - the leading 0.5 is from 1/2*rho
    // tensor forchCoeff(Zero);
    // forchCoeff.xx() = 0.5*fXYZ_.value().x();
    // forchCoeff.yy() = 0.5*fXYZ_.value().y();
    // forchCoeff.zz() = 0.5*fXYZ_.value().z();

    // if (csys().uniform())
    // {
    //     forAll(cellZoneIDs_, zonei)
    //     {
    //         D_[zonei].resize(1);
    //         F_[zonei].resize(1);

    //         D_[zonei] = csys().transform(darcyCoeff);
    //         F_[zonei] = csys().transform(forchCoeff);

    //         // Rsh, 2023-04-25
    //         // virtual scalar transform (const scalar &input) const
    //         // With constant rotation tensor.
    //     }
    // }
    // else
    // {
    //     forAll(cellZoneIDs_, zonei)
    //     {
    //         const pointUIndList cc
    //         (
    //             mesh_.cellCentres(),
    //             mesh_.cellZones()[cellZoneIDs_[zonei]]
    //         );

    //         D_[zonei] = csys().transform(cc, darcyCoeff);
    //         F_[zonei] = csys().transform(cc, forchCoeff);

    //         // Rsh, 2023-04-25
    //         // virtual tmp< Field< scalar > > transform (const pointUIndList &global, const scalar &input) const
    //     }
    // }

    // OFstream D1_out(mesh_.time().caseConstant()+"D1");
    // D_[0].writeEntry("D1", D1_out);

    // OFstream F1_out(mesh_.time().caseConstant()+"F1");
    // F_[0].writeEntry("F1", F1_out);

    // Info<< "D_1" << D_ << endl
    //     << "F_1" << F_ << endl;

    forAll(cellZoneIDs_, zonei)
    {
        D_[zonei].resize(mesh_.cellZones()[cellZoneIDs_[zonei]].size());
        F_[zonei].resize(mesh_.cellZones()[cellZoneIDs_[zonei]].size());

        forAll(mesh_.cellZones()[cellZoneIDs_[zonei]], i)
        {
            const label celli = mesh_.cellZones()[cellZoneIDs_[zonei]][i];
            
            D_[zonei][i] = Dfield_[celli];
            F_[zonei][i] = Ffield_[celli];
        }            
    }

    // OFstream D2_out(mesh_.time().caseConstant()+"/D2");
    // D_[0].writeEntry("D2", D2_out);

    // OFstream F2_out(mesh_.time().caseConstant()+"/F2");
    // F_[0].writeEntry("F2", F2_out);

    // Info<< "D_2\n" << D_ << endl
    //     << "F_2\n" << F_ << endl;

    // // Rsh, 2023-04-18

    // // Rsh, 2023-04-25
    // if (csys().uniform())
    // {
    //     forAll(cellZoneIDs_, zonei)
    //     {
    //         D_[zonei].resize(1);
    //         F_[zonei].resize(1);

    //         D_[zonei] = csys().transform(darcyCoeff);
    //         F_[zonei] = csys().transform(forchCoeff);

    //         // Rsh, 2023-04-25
    //         // virtual scalar transform (const scalar &input) const
    //         // With constant rotation tensor.
    
    //         // void 	shallowCopy (const UList< T > &list)
    //         // Copy the pointer and size held by the given UList.
    //         // void 	deepCopy (const UList< T > &list)
    //         // Copy elements of the given UList. Sizes must match!
    //     }
    // }
    // else
    // {
    //     forAll(cellZoneIDs_, zonei)
    //     {
    //         const pointUIndList cc
    //         (
    //             mesh_.cellCentres(),
    //             mesh_.cellZones()[cellZoneIDs_[zonei]]
    //         );

    //         D_[zonei] = csys().transform(cc, darcyCoeff);
    //         F_[zonei] = csys().transform(cc, forchCoeff);

    //         // Rsh, 2023-04-25
    //         // virtual tmp< Field< scalar > > transform (const pointUIndList &global, const scalar &input) const
    //     }
    // }
    // // Rsh, 2023-04-25


    if (debug && mesh_.time().writeTime())
    {
        volTensorField Dout
        (
            IOobject
            (
                typeName + ":D",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(dXYZ_.dimensions(), Zero)
        );
        volTensorField Fout
        (
            IOobject
            (
                typeName + ":F",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedTensor(fXYZ_.dimensions(), Zero)
        );


        forAll(cellZoneIDs_, zonei)
        {
            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zonei]];

            if (csys().uniform())
            {
                UIndirectList<tensor>(Dout, cells) = D_[zonei].first();
                UIndirectList<tensor>(Fout, cells) = F_[zonei].first();
            }
            else
            {
                UIndirectList<tensor>(Dout, cells) = D_[zonei];
                UIndirectList<tensor>(Fout, cells) = F_[zonei];
            }
        }

        Dout.write();
        Fout.write();
    }
}


void Foam::porosityModels::nonUniformDarcyForchheimer::calcForce
(
    const volVectorField& U,
    const volScalarField& rho,
    const volScalarField& mu,
    vectorField& force
) const
{
    scalarField Udiag(U.size(), Zero);
    vectorField Usource(U.size(), Zero);
    const scalarField& V = mesh_.V();

    apply(Udiag, Usource, V, rho, mu, U);

    force = Udiag*U - Usource;
}


void Foam::porosityModels::nonUniformDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn
) const
{
    const volVectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>(rhoName);

        if (mesh_.foundObject<volScalarField>(muName))
        {
            const auto& mu = mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, rho, mu, U);
        }
        else
        {
            const auto& nu = mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, rho, rho*nu, U);
        }
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const auto& nu = mesh_.lookupObject<volScalarField>(nuName);

            apply(Udiag, Usource, V, geometricOneField(), nu, U);
        }
        else
        {
            const auto& rho = mesh_.lookupObject<volScalarField>(rhoName);
            const auto& mu = mesh_.lookupObject<volScalarField>(muName);

            apply(Udiag, Usource, V, geometricOneField(), mu/rho, U);
        }
    }
}


void Foam::porosityModels::nonUniformDarcyForchheimer::correct
(
    fvVectorMatrix& UEqn,
    const volScalarField& rho,
    const volScalarField& mu
) const
{
    const vectorField& U = UEqn.psi();
    const scalarField& V = mesh_.V();
    scalarField& Udiag = UEqn.diag();
    vectorField& Usource = UEqn.source();

    apply(Udiag, Usource, V, rho, mu, U);
}


void Foam::porosityModels::nonUniformDarcyForchheimer::correct
(
    const fvVectorMatrix& UEqn,
    volTensorField& AU
) const
{
    const volVectorField& U = UEqn.psi();

    word rhoName(IOobject::groupName(rhoName_, U.group()));
    word muName(IOobject::groupName(muName_, U.group()));
    word nuName(IOobject::groupName(nuName_, U.group()));

    if (UEqn.dimensions() == dimForce)
    {
        const auto& rho = mesh_.lookupObject<volScalarField>(rhoName);
        const auto& mu = mesh_.lookupObject<volScalarField>(muName);

        apply(AU, rho, mu, U);
    }
    else
    {
        if (mesh_.foundObject<volScalarField>(nuName))
        {
            const auto& nu = mesh_.lookupObject<volScalarField>(nuName);

            apply(AU, geometricOneField(), nu, U);
        }
        else
        {
            const auto& rho = mesh_.lookupObject<volScalarField>(rhoName);
            const auto& mu = mesh_.lookupObject<volScalarField>(muName);

            apply(AU, geometricOneField(), mu/rho, U);
        }
    }
}


bool Foam::porosityModels::nonUniformDarcyForchheimer::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);

    return true;
}


// ************************************************************************* //
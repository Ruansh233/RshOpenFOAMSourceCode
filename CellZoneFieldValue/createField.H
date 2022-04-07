
Info << "fieldName: " << fieldName_ << endl;

if (fieldName_ == "U")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(U.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}
else if (fieldName_ == "Ux")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField field_(U.component(vector::X));
    #include "writeData.H"
}
else if (fieldName_ == "Uy")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField field_(U.component(vector::Y));
    #include "writeData.H"
}
else if (fieldName_ == "Uz")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField field_(U.component(vector::Z));
    #include "writeData.H"
}
else if (fieldName_ == "magU")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volScalarField field_(mag(U));
    #include "writeData.H"
}
else if (fieldName_ == "gradUx")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField gradUx(fvc::grad(U.component(vector::X)));

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(gradUx.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}
else if (fieldName_ == "gradUy")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField gradUy(fvc::grad(U.component(vector::Y)));

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(gradUy.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}
else if (fieldName_ == "gradUz")
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField gradUz(fvc::grad(U.component(vector::Z)));

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(gradUz.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}



else if (fieldName_ == "p")
{
    volScalarField field_ // note that pressure is a scalar field
	(
		IOobject
		(
		    "p", // name of the field
		    runTimeName_, // name of the current time, i.e. the time folder to read from
		    mesh,
		    IOobject::MUST_READ, // always gets imported, will throw an error if the field is missing
		    IOobject::NO_WRITE // will get saved automatically when the controlDict parameters will request it
		),
		mesh // initialises the field to match the size of the mesh with default (0) values
	);

    #include "writeData.H"
}

else if (fieldName_ == "gradp")
{
    volScalarField p // note that velocity is a vector field
    (
        IOobject
        (
            "p",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField gradp(fvc::grad(p));

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(gradp.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}

else
{
    volVectorField U // note that velocity is a vector field
    (
        IOobject
        (
            "U",
            runTimeName_,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    for(label i = 0; i <= 2; i++)
    {
        volScalarField field_(U.component(i));
        word fileName_(fieldName_ + name(i));
        #include "writeData.H"
    }

    // writeDataToFile(fieldName_ + "x", field_.component(vector::X), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "y", field_.component(vector::Y), cellZonesName_, mesh, outputDir);
    // writeDataToFile(fieldName_ + "z", field_.component(vector::Z), cellZonesName_, mesh, outputDir);
}

    
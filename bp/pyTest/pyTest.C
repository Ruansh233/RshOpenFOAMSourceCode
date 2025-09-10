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
#include "IFstream.H"
#include <pybind11/embed.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

int main(int argc, char *argv[])
{
    // OpenFOAM setup
    argList::noParallel();
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    List<scalar> data(100);
    forAll(data, i)
    {
        data[i] = scalar(i) * 0.5;
    }

    // Compute the sum
    Info<< "Sum of data: " << sum(data) << endl;

    std::vector<double> vecData(data.size());
    forAll(data, i)
    {
        vecData[i] = data[i];
    }

    py::array_t<double> np_array(vecData.size(), vecData.data());

    py::module myscript = py::module::import("add");

    try
    {
        py::object py_result = myscript.attr("add_numbers")(np_array);
        double result = py_result.cast<double>();
        Info << "Result from Python: " << result << nl << endl;
    }
    catch (const std::exception &e)
    {
        Info << "Python error: " << e.what() << nl << endl;
    }

    return 0;
}

// ************************************************************************* //

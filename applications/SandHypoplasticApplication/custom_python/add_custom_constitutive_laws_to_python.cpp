//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/path_a_fortran_dll/sand_hypoplastic_fortran_dll_law.h"
#include "custom_constitutive/path_b_cpp_port/sand_hypoplastic_cpp_law.h"
#include "custom_constitutive/path_b_cpp_port/sand_hypoplastic_cpp_law_2d_plane_strain.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    py::class_<SandHypoplasticFortranDllLaw, SandHypoplasticFortranDllLaw::Pointer, ConstitutiveLaw>(m, "SandHypoplasticFortranDllLaw")
        .def(py::init<>());

    py::class_<SandHypoplasticCppLaw, SandHypoplasticCppLaw::Pointer, ConstitutiveLaw>(m, "SandHypoplasticCppLaw")
        .def(py::init<>());

    py::class_<SandHypoplasticCppLaw2DPlaneStrain, SandHypoplasticCppLaw2DPlaneStrain::Pointer, ConstitutiveLaw>(m, "SandHypoplasticCppLaw2DPlaneStrain")
        .def(py::init<>());
}

} // namespace Python
} // namespace Kratos
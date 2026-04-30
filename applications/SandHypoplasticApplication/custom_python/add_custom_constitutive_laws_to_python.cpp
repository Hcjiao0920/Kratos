//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

// Project includes
#include "custom_python/add_custom_constitutive_laws_to_python.h"
#include "custom_constitutive/path_a_fortran_dll/sand_hypoplastic_fortran_dll_law.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void AddCustomConstitutiveLawsToPython(pybind11::module& m)
{
    py::class_<SandHypoplasticFortranDllLaw, SandHypoplasticFortranDllLaw::Pointer, ConstitutiveLaw>(m, "SandHypoplasticFortranDllLaw")
        .def(py::init<>());
}

} // namespace Python
} // namespace Kratos
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

#if defined(KRATOS_PYTHON)

// Project includes
#include "includes/define_python.h"
#include "sand_hypoplastic_application.h"
#include "sand_hypoplastic_application_variables.h"
#include "custom_python/add_custom_constitutive_laws_to_python.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

PYBIND11_MODULE(KratosSandHypoplasticApplication, m)
{
    py::class_<KratosSandHypoplasticApplication, KratosSandHypoplasticApplication::Pointer, KratosApplication>(m, "KratosSandHypoplasticApplication")
        .def(py::init<>());

    AddCustomConstitutiveLawsToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SAND_HYPO_FORTRAN_DLL_PATH)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SAND_HYPOPLASTIC_PROPS_16)
}

} // namespace Python
} // namespace Kratos

#endif // KRATOS_PYTHON defined
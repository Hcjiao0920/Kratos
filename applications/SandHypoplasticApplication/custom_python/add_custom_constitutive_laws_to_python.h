//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

#if !defined(KRATOS_SAND_HYPOPLASTIC_ADD_CUSTOM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED)
#define KRATOS_SAND_HYPOPLASTIC_ADD_CUSTOM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED

// Project includes
#include "includes/define_python.h"

namespace Kratos
{
namespace Python
{
void AddCustomConstitutiveLawsToPython(pybind11::module& m);
} // namespace Python
} // namespace Kratos

#endif // KRATOS_SAND_HYPOPLASTIC_ADD_CUSTOM_CONSTITUTIVE_LAWS_TO_PYTHON_H_INCLUDED
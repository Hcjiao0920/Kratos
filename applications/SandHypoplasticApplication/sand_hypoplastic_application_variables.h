//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

#if !defined(KRATOS_SAND_HYPOPLASTIC_APPLICATION_VARIABLES_H_INCLUDED)
#define KRATOS_SAND_HYPOPLASTIC_APPLICATION_VARIABLES_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/variables.h"

namespace Kratos
{
    KRATOS_DEFINE_APPLICATION_VARIABLE(SAND_HYPOPLASTIC_APPLICATION, std::string, SAND_HYPO_FORTRAN_DLL_PATH)
    KRATOS_DEFINE_APPLICATION_VARIABLE(SAND_HYPOPLASTIC_APPLICATION, Vector, SAND_HYPOPLASTIC_PROPS_16)
}

#endif // KRATOS_SAND_HYPOPLASTIC_APPLICATION_VARIABLES_H_INCLUDED
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

// Project includes
#include "sand_hypoplastic_application.h"

namespace Kratos
{

KratosSandHypoplasticApplication::KratosSandHypoplasticApplication()
    : KratosApplication("SandHypoplasticApplication")
{
}

void KratosSandHypoplasticApplication::Register()
{
    std::cout << "Initializing KratosSandHypoplasticApplication..." << std::endl;

    KRATOS_REGISTER_VARIABLE(SAND_HYPO_FORTRAN_DLL_PATH)
    KRATOS_REGISTER_VARIABLE(SAND_HYPOPLASTIC_PROPS_16)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SandHypoplasticFortranDllLaw", mSandHypoplasticFortranDllLaw)
}

} // namespace Kratos
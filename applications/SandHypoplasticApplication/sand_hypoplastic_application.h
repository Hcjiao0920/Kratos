//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics

#if !defined(KRATOS_SAND_HYPOPLASTIC_APPLICATION_H_INCLUDED)
#define KRATOS_SAND_HYPOPLASTIC_APPLICATION_H_INCLUDED

// Project includes
#include "includes/kratos_application.h"
#include "sand_hypoplastic_application_variables.h"
#include "custom_constitutive/path_a_fortran_dll/sand_hypoplastic_fortran_dll_law.h"
#include "custom_constitutive/path_b_cpp_port/sand_hypoplastic_cpp_law.h"

namespace Kratos
{

class KRATOS_API(SAND_HYPOPLASTIC_APPLICATION) KratosSandHypoplasticApplication : public KratosApplication
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(KratosSandHypoplasticApplication);

    KratosSandHypoplasticApplication();
    ~KratosSandHypoplasticApplication() override = default;

    void Register() override;

    std::string Info() const override
    {
        return "KratosSandHypoplasticApplication";
    }

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

private:
    const SandHypoplasticFortranDllLaw mSandHypoplasticFortranDllLaw;
    const SandHypoplasticCppLaw        mSandHypoplasticCppLaw;
};

} // namespace Kratos

#endif // KRATOS_SAND_HYPOPLASTIC_APPLICATION_H_INCLUDED

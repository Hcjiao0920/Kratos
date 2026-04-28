//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Path A — dynamic-link Fortran UMAT bridge for the
//  von Wolffersdorff (1996) sand hypoplastic + Niemunis-Herle (1997)
//  intergranular strain model. The Fortran kernel lives in
//  custom_constitutive/path_a_fortran_dll/fortran/sand_hypo_is.dll
//  and exposes the standard Abaqus UMAT entry point `umat_`.
//
//  This class only handles the C++ <-> Fortran bridge: deformation-gradient
//  history, polar-decomposition based co-rotational wrap, Voigt-order swap
//  between Kratos MPM ([11,22,33,12,23,13]) and Abaqus ([11,22,33,12,13,23]),
//  Fortran-style column-major DDSDDE transpose, and the `inittension`
//  tensile-elastic fallback. The constitutive physics lives entirely in
//  the Fortran kernel.

#if !defined(KRATOS_SAND_HYPOPLASTIC_FORTRAN_DLL_LAW_H_INCLUDED)
#define       KRATOS_SAND_HYPOPLASTIC_FORTRAN_DLL_LAW_H_INCLUDED

// System includes
#include <cstddef>
#include <string>

// Project includes
#include "includes/constitutive_law.h"

// Application includes
#include "custom_constitutive/hyperelastic_3D_law.hpp"

namespace Kratos
{

// Abaqus UMAT 37-argument signature. The Fortran kernel is built with
// gfortran (MinGW-w64 UCRT POSIX SEH 15.2.0) using the default ABI; on
// x86-64 Windows there is a single calling convention so no __stdcall
// decoration is needed (verified in poc_fortran_bridge/).
//
// Hidden trailing size_t = declared length of `cmname` (character*80).
// Omitting it corrupts the stack on gfortran-built DLLs.
//
// `ddsdde` is a single contiguous 36-double column-major buffer; the bridge
// transposes it into a row-major Kratos `Matrix` after the call.
using SandHypoFortranUmatPtr = void (*)(
    double* stress, double* statev, double* ddsdde,
    double* sse, double* spd, double* scd,
    double* rpl, double* ddsddt, double* drplde, double* drpldt,
    double* stran, double* dstran, double* time, double* dtime,
    double* temp, double* dtemp, double* predef, double* dpred,
    char*   cmname,
    int*    ndi, int* nshr, int* ntens, int* nstatv,
    double* props, int* nprops,
    double* coords, double* drot, double* pnewdt, double* celent,
    double* dfgrd0, double* dfgrd1,
    int*    noel, int* npt, int* layer, int* kspt, int* kstep, int* kinc,
    std::size_t cmname_len);

class KRATOS_API(MPM_APPLICATION) SandHypoplasticFortranDllLaw
    : public HyperElastic3DLaw
{
public:
    using BaseType = HyperElastic3DLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SandHypoplasticFortranDllLaw);

    SandHypoplasticFortranDllLaw();
    SandHypoplasticFortranDllLaw(const SandHypoplasticFortranDllLaw& rOther);
    ~SandHypoplasticFortranDllLaw() override = default;

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override { return 3; }
    SizeType GetStrainSize() const override   { return 6; }

    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Deformation_Gradient; }
    StressMeasure GetStressMeasure() override { return StressMeasure_Kirchhoff; }

    void GetLawFeatures(Features& rFeatures) override;

    // Expose mStateVarsFinalized via INTERNAL_VARIABLES so regression tests
    // can verify intergranular-strain co-rotation. The slot semantics are
    // documented next to mStateVarsFinalized below; in particular slots
    // 0..5 are intergranular strain components in Kratos Voigt order
    // [del_11, del_22, del_33, del_12, del_23, del_13].
    bool Has(const Variable<Vector>& rThisVariable) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    // Direct injection of mStateVarsFinalized via INTERNAL_VARIABLES. Used
    // by regression tests that need to set up specific intergranular-strain
    // states without driving them through a real loading history.
    void SetValue(const Variable<Vector>& rThisVariable,
                  const Vector& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) const override;

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    void ResetMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues) override;

    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    // Helper used by both finalize overrides after the post-Calculate
    // recompute: copy trial state into the finalized slots and roll the
    // F-history forward. Kept outside the overrides so the three-stage
    // (set flag -> Calculate -> reset flag -> commit) shape stays
    // visible at a glance.
    void CommitTrialAndRollFHistory(Parameters& rValues);

    std::string Info() const override { return "SandHypoplasticFortranDllLaw"; }
    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:
    // Cauchy stress at the end of the previous (committed) step, in Kratos
    // Voigt order [11,22,33,12,23,13]. The bridge converts to Abaqus order
    // before each UMAT call. Stored Cauchy (not Kirchhoff) because that is
    // the natural UMAT input/output.
    Vector mStressVectorCauchyFinalized;

    // Full UMAT statev array (15 doubles per project_umat_model.md):
    //   0..5   intergranular strain (continuum sign, KRATOS Voigt order
    //          [del_11, del_22, del_33, del_12, del_23, del_13]) — the
    //          bridge swaps idx 4<->5 just before/after the UMAT call to
    //          translate to Abaqus order. Pre-rotated (R · IS · R^T) on the
    //          way into UMAT, returned in NEW lab frame (Abaqus convention).
    //   6      void ratio e
    //   7..14  postprocessing slots written by the kernel: excess pore
    //          pressure, mean p, nfev, phi_mob, rho, dtsub, padding
    // The first 7 slots are integrated by the RKF23 driver; entries 7..14
    // are written by the kernel for diagnostics only.
    Vector mStateVarsFinalized;

    // Trial state — written by CalculateMaterialResponseKirchhoff during
    // Newton iterations, committed to the *Finalized members in
    // FinalizeMaterialResponseKirchhoff. Same units, same Voigt order, same
    // length conventions as the *Finalized counterparts.
    Vector mStressVectorCauchyTrial;
    Vector mStateVarsTrial;

    // True once the DLL has been resolved successfully for at least one
    // instance of this law in the current process.
    bool mUmatResolved = false;

private:
    // Process-lifetime DLL handle. LoadLibrary is called once per process
    // (lazy on first InitializeMaterial); never FreeLibrary'd. Stored as
    // void* so this header does not need <windows.h>; the .cpp casts to
    // HMODULE on Windows or to a dlopen handle on POSIX.
    static SandHypoFortranUmatPtr s_umat_ptr;
    static void*                  s_dll_handle;
    static std::string            s_loaded_dll_path;

    // Idempotent: each instance calls this in InitializeMaterial; the first
    // call resolves the DLL, later calls fast-return.
    void EnsureUmatLoaded(const Properties& rMaterialProperties);

    // OS-specific loader bodies; both look up the symbol "umat_"
    // (gfortran's default name-mangled export).
    bool LoadUmatWindows(const std::string& rDllPath);
    bool LoadUmatLinux  (const std::string& rDllPath);

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, HyperElastic3DLaw)
        rSerializer.save("StressVectorCauchyFinalized", mStressVectorCauchyFinalized);
        rSerializer.save("StateVarsFinalized",          mStateVarsFinalized);
        rSerializer.save("StressVectorCauchyTrial",     mStressVectorCauchyTrial);
        rSerializer.save("StateVarsTrial",              mStateVarsTrial);
        rSerializer.save("UmatResolved",                mUmatResolved);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, HyperElastic3DLaw)
        rSerializer.load("StressVectorCauchyFinalized", mStressVectorCauchyFinalized);
        rSerializer.load("StateVarsFinalized",          mStateVarsFinalized);
        rSerializer.load("StressVectorCauchyTrial",     mStressVectorCauchyTrial);
        rSerializer.load("StateVarsTrial",              mStateVarsTrial);
        rSerializer.load("UmatResolved",                mUmatResolved);
    }

}; // class SandHypoplasticFortranDllLaw

} // namespace Kratos

#endif // KRATOS_SAND_HYPOPLASTIC_FORTRAN_DLL_LAW_H_INCLUDED

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Path A — sand_hypoplastic_fortran_dll_law.cpp
//
//  Step 3 SKELETON. The DLL loader plumbing is fully wired; the
//  CalculateMaterialResponseKirchhoff bridge is intentionally stubbed
//  (returns zero Kirchhoff stress + identity 6x6 tangent and emits a
//  one-time KRATOS_WARNING). The full F -> R_inc -> DSTRAN -> UMAT ->
//  undo-rotation algorithm, the Voigt sigma_xz<->sigma_yz swap, the
//  column-major DDSDDE transpose, and the inittension tensile-elastic
//  fallback are all Step 4 of memory/project_path_a_fortran_dll.md.

// System includes
#include <algorithm>
#include <cmath>
#include <cstring>
#include <sstream>
#include <vector>

// Project includes (must come BEFORE the platform DLL-loader block, because
// includes/define.h is what sets KRATOS_COMPILED_IN_{WINDOWS,LINUX,OS} from
// the compiler's own _WIN32 / __linux__ / __APPLE__ predefined macros).
#include "includes/define.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_constitutive/path_a_fortran_dll/sand_hypoplastic_fortran_dll_law.h"
#include "mpm_application_variables.h"

// Platform DLL loader
#ifdef KRATOS_COMPILED_IN_WINDOWS
    #ifndef NOMINMAX
        #define NOMINMAX
        #include <windows.h>
        #undef NOMINMAX
    #else
        #include <windows.h>
    #endif
#elif defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
    #include <dlfcn.h>
#endif

namespace Kratos
{

// Process-lifetime statics. Once the DLL is loaded for any one material in
// the process, every other particle re-uses the same function pointer.
SandHypoFortranUmatPtr SandHypoplasticFortranDllLaw::s_umat_ptr        = nullptr;
void*                  SandHypoplasticFortranDllLaw::s_dll_handle      = nullptr;
std::string            SandHypoplasticFortranDllLaw::s_loaded_dll_path = "";

// ----- Construction / cloning ----------------------------------------------

SandHypoplasticFortranDllLaw::SandHypoplasticFortranDllLaw()
    : HyperElastic3DLaw(),
      mStressVectorCauchyFinalized(ZeroVector(6)),
      mStateVarsFinalized(ZeroVector(15)),
      mStressVectorCauchyTrial(ZeroVector(6)),
      mStateVarsTrial(ZeroVector(15))
{
}

SandHypoplasticFortranDllLaw::SandHypoplasticFortranDllLaw(
    const SandHypoplasticFortranDllLaw& rOther)
    : HyperElastic3DLaw(rOther),
      mStressVectorCauchyFinalized(rOther.mStressVectorCauchyFinalized),
      mStateVarsFinalized(rOther.mStateVarsFinalized),
      mStressVectorCauchyTrial(rOther.mStressVectorCauchyTrial),
      mStateVarsTrial(rOther.mStateVarsTrial),
      mUmatResolved(rOther.mUmatResolved)
{
}

ConstitutiveLaw::Pointer SandHypoplasticFortranDllLaw::Clone() const
{
    return Kratos::make_shared<SandHypoplasticFortranDllLaw>(*this);
}

// ----- Law features --------------------------------------------------------

bool SandHypoplasticFortranDllLaw::Has(const Variable<Vector>& rThisVariable)
{
    if (rThisVariable == INTERNAL_VARIABLES) return true;
    return BaseType::Has(rThisVariable);
}

Vector& SandHypoplasticFortranDllLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        if (rValue.size() != mStateVarsFinalized.size())
            rValue.resize(mStateVarsFinalized.size(), false);
        noalias(rValue) = mStateVarsFinalized;
        return rValue;
    }
    return BaseType::GetValue(rThisVariable, rValue);
}

void SandHypoplasticFortranDllLaw::SetValue(const Variable<Vector>& rThisVariable,
                                            const Vector& rValue,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        KRATOS_ERROR_IF(rValue.size() != 15)
            << "SandHypoplasticFortranDllLaw::SetValue(INTERNAL_VARIABLES) requires "
               "a 15-element vector (6 IS components in Kratos Voigt order, void "
               "ratio at index 6, 8 postprocessing slots); got size "
            << rValue.size() << "." << std::endl;
        if (mStateVarsFinalized.size() != 15) mStateVarsFinalized.resize(15, false);
        if (mStateVarsTrial.size() != 15)     mStateVarsTrial.resize(15, false);
        noalias(mStateVarsFinalized) = rValue;
        noalias(mStateVarsTrial) = rValue;
        return;
    }
    BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
}

void SandHypoplasticFortranDllLaw::GetLawFeatures(Features& rFeatures)
{
    // FINITE_STRAINS is mandatory for MPM compatibility:
    // mpm_updated_lagrangian.cpp's strain-measure check rejects laws that
    // declare INFINITESIMAL_STRAINS (this is why GMA's SmallStrainUMAT3DLaw
    // cannot be used directly for MPM — see docs/kratos_vs_abaqus_kinematics.md).
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(ISOTROPIC);

    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize    = GetStrainSize();
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

// ----- Check ---------------------------------------------------------------

int SandHypoplasticFortranDllLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPO_FORTRAN_DLL_PATH))
        << "SandHypoplasticFortranDllLaw requires Properties[SAND_HYPO_FORTRAN_DLL_PATH] "
        << "(absolute or build-relative path to sand_hypo_is.dll)." << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(MATERIAL_PARAMETERS))
        << "SandHypoplasticFortranDllLaw requires Properties[MATERIAL_PARAMETERS] "
        << "containing the 16 von Wolffersdorff + Niemunis-Herle parameters "
        << "(see project_umat_model.md for the slot layout)." << std::endl;

    const Vector& r_params = rMaterialProperties[MATERIAL_PARAMETERS];
    KRATOS_ERROR_IF(r_params.size() < 16)
        << "MATERIAL_PARAMETERS must have >= 16 entries; got " << r_params.size() << "."
        << std::endl;

    // Finite-value hygiene: NaN/Inf in MATERIAL_PARAMETERS would slip past
    // the magnitude/sign comparisons below (NaN compares false to anything),
    // then propagate into the UMAT and contaminate stress/state. Codex
    // adversarial review 2026-04-29 round 11 Finding 2 — mirrors the J<=0
    // guard pattern already used in CalculateMaterialResponseKirchhoff.
    for (std::size_t i = 0; i < r_params.size() && i < 16; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(r_params[i]))
            << "MATERIAL_PARAMETERS[" << i << "] is non-finite (got "
            << r_params[i] << "). NaN/Inf in material parameters cannot be "
               "validated by domain checks and would propagate through the "
               "UMAT into committed material history." << std::endl;
    }

    // Mirror the Fortran kernel's check_parms_h domain checks BEFORE any UMAT
    // call. The kernel's xit_h on error=10 executes Fortran STOP, which
    // terminates the entire Kratos process; surfacing those rejections as
    // KRATOS_ERROR here makes them recoverable. Codex adversarial review
    // 2026-04-29 round 10 Finding 1 / TD-2 closure. The 8 checked params
    // match the Fortran source 1-for-1; other slots (hs, en, ed0, ec0, ei0,
    // alpha, beta, e0_param) have no explicit Fortran domain check, so we
    // pass them through. Slot indices match project_umat_model.md.
    KRATOS_ERROR_IF(r_params[0] <= 0.0)
        << "MATERIAL_PARAMETERS[0] (phi, friction angle in degrees) must be > 0; got "
        << r_params[0] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[1] < 0.0)
        << "MATERIAL_PARAMETERS[1] (p_t, tension cut-off) must be >= 0; got "
        << r_params[1] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[9] < 0.0)
        << "MATERIAL_PARAMETERS[9] (m_R, intergranular-strain reversal stiffness multiplier) "
           "must be >= 0; got " << r_params[9] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[10] < 0.0)
        << "MATERIAL_PARAMETERS[10] (m_T, intergranular-strain 90-degree multiplier) "
           "must be >= 0; got " << r_params[10] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[11] < 0.0)
        << "MATERIAL_PARAMETERS[11] (r_uc, intergranular-strain reference length) "
           "must be >= 0; got " << r_params[11] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[12] < 0.0)
        << "MATERIAL_PARAMETERS[12] (beta_r, intergranular-strain evolution exponent) "
           "must be >= 0; got " << r_params[12] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[13] < 0.0)
        << "MATERIAL_PARAMETERS[13] (chi, intergranular-strain interpolation exponent) "
           "must be >= 0; got " << r_params[13] << "." << std::endl;
    KRATOS_ERROR_IF(r_params[14] < 0.0)
        << "MATERIAL_PARAMETERS[14] (bulk_w, water bulk modulus) "
           "must be >= 0; got " << r_params[14] << "." << std::endl;

    // INITIAL_STRESS_VECTOR is recommended but not required: a missing or
    // zero initial stress is physically meaningful in some test setups
    // (e.g. the tensile_fallback case). Warn rather than error so the law
    // can still be used with sigma_old = 0 if the caller really wants that.
    if (rMaterialProperties.Has(INITIAL_STRESS_VECTOR)) {
        const Vector& s0 = rMaterialProperties[INITIAL_STRESS_VECTOR];
        KRATOS_ERROR_IF(s0.size() != 6)
            << "INITIAL_STRESS_VECTOR must have exactly 6 components in Kratos Voigt order "
               "[11,22,33,12,23,13]; got " << s0.size() << "." << std::endl;
        for (std::size_t i = 0; i < s0.size(); ++i) {
            KRATOS_ERROR_IF_NOT(std::isfinite(s0[i]))
                << "INITIAL_STRESS_VECTOR[" << i << "] is non-finite (got " << s0[i]
                << "). NaN/Inf cannot be a valid stress component." << std::endl;
        }
    } else {
        KRATOS_WARNING_FIRST_N("SandHypoplasticFortranDllLaw", 1)
            << "Properties[INITIAL_STRESS_VECTOR] not provided; sigma_old defaults to 0. "
               "Sand hypoplastic is stress-state dependent — set INITIAL_STRESS_VECTOR "
               "to the geostatic / hydrostatic initial Cauchy stress (tension-positive)." << std::endl;
    }

    return 0;
}

// ----- Initialization ------------------------------------------------------

void SandHypoplasticFortranDllLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    // Defense-in-depth: run the same validation Check() does, so callers that
    // skip Kratos's normal Check stage (e.g. standalone test drivers) cannot
    // sneak invalid material parameters into the UMAT call. ProcessInfo is
    // unused by our Check; pass a default-constructed one.
    {
        const ProcessInfo dummy_info;
        Check(rMaterialProperties, rElementGeometry, dummy_info);
    }

    HyperElastic3DLaw::InitializeMaterial(
        rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    if (mStressVectorCauchyFinalized.size() != 6)
        mStressVectorCauchyFinalized = ZeroVector(6);
    if (mStateVarsFinalized.size() != 15)
        mStateVarsFinalized = ZeroVector(15);
    if (mStressVectorCauchyTrial.size() != 6)
        mStressVectorCauchyTrial = ZeroVector(6);
    if (mStateVarsTrial.size() != 15)
        mStateVarsTrial = ZeroVector(15);

    // Seed initial Cauchy stress from Properties (Kratos Voigt order, tension+).
    // Without this seeding every particle starts at sigma = 0, which puts the
    // hypoplastic kernel into a numerically degenerate region near the origin
    // of the response envelope and (in tension+ convention) is also where the
    // Fortran inittension fallback can fire spuriously.
    if (rMaterialProperties.Has(INITIAL_STRESS_VECTOR)) {
        const Vector& s0 = rMaterialProperties[INITIAL_STRESS_VECTOR];
        if (s0.size() == 6) {
            for (int i = 0; i < 6; ++i) mStressVectorCauchyFinalized[i] = s0[i];
        }
    }

    // Initial void ratio e0:
    //   * The Fortran kernel itself initializes statev(7) on its FIRST UMAT
    //     call when statev(7) < 0.001, computing the pressure-dependent
    //     expression `props[15] * exp(-(3p/hs)^n)` if props[15] <= 10, or
    //     `props[15] - 10` if > 10 (literal). See 00_hypo_is_original.for
    //     line ~171-185 for the exact formula.
    //   * Path A used to pre-seed statev[6] from props[15] in this method,
    //     which inadvertently skipped the kernel's pressure-dependent
    //     branch — every particle started with a literal e0 (e.g. 0.65)
    //     instead of the calibrated pressure-corrected value (e.g. ~0.5975
    //     at p = 100 kPa for Hostun dense). Codex adversarial review
    //     2026-04-29 round 7 Finding 1.
    //   * Fix: leave statev[6] at zero (the constructor / ResetMaterial
    //     default). The kernel's init branch fires on the first Calculate
    //     and overwrites statev_abq[6] with the correct e0.
    //   * Verified by tests/path_a_results/test_void_ratio_init.py.

    EnsureUmatLoaded(rMaterialProperties);
}

void SandHypoplasticFortranDllLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    mStressVectorCauchyFinalized = ZeroVector(6);
    mStateVarsFinalized          = ZeroVector(15);
    mUmatResolved                = false;
    InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
}

// ----- Bridge layer (Step 4) ----------------------------------------------
//
// Conventions used throughout this section:
//   * Voigt orderings:
//       Kratos MPM stress/strain: [11, 22, 33, 12, 23, 13]   (idx 4 = yz, 5 = xz)
//       Abaqus 3D UMAT:           [11, 22, 33, 12, 13, 23]   (idx 4 = xz, 5 = yz)
//     -> swap idx 4 <-> 5 on every boundary (stress in/out, DSTRAN in,
//        statev intergranular components in/out, AND both rows AND cols
//        of DDSDDE).
//   * Engineering-shear convention: a Voigt strain vector stores
//        eps_voigt[3] = 2 * eps_tensor[0,1]  (and similar for [4],[5])
//     so the 3x3 <-> 6 conversions for STRAIN apply factors of 1/2 / 2;
//     STRESS conversions do not.
//   * Stored Cauchy stress is sigma in spatial (lab) frame. The bridge
//     pre-rotates it by R_inc (polar decomp of F_inc = F_new * F_old^-1)
//     into the corotational frame the UMAT expects.
//   * `mStateVarsTrial[0..5]` holds the intergranular strain in Abaqus
//     component order (matching the Fortran `define_h` convention) so the
//     UMAT call is just a memcpy + 4<->5 swap. We swap back into Kratos
//     order only when the value is exposed externally.

namespace
{

// Convert a 3x3 symmetric stress tensor to a 6-Voigt vector in Kratos
// MPM order [xx, yy, zz, xy, yz, xz] (no engineering-shear factor).
void StressTensorToKratosVoigt(const Matrix& T, Vector& v)
{
    v[0] = T(0, 0); v[1] = T(1, 1); v[2] = T(2, 2);
    v[3] = T(0, 1); v[4] = T(1, 2); v[5] = T(0, 2);
}

// Inverse of StressTensorToKratosVoigt.
void KratosVoigtToStressTensor(const Vector& v, Matrix& T)
{
    T(0, 0) = v[0]; T(1, 1) = v[1]; T(2, 2) = v[2];
    T(0, 1) = T(1, 0) = v[3];
    T(1, 2) = T(2, 1) = v[4];
    T(0, 2) = T(2, 0) = v[5];
}

// Convert a 3x3 symmetric STRAIN tensor to a 6-Voigt vector in Kratos
// MPM order with engineering-shear factor of 2 on off-diagonals.
void StrainTensorToKratosVoigt(const Matrix& T, Vector& v)
{
    v[0] = T(0, 0); v[1] = T(1, 1); v[2] = T(2, 2);
    v[3] = 2.0 * T(0, 1);
    v[4] = 2.0 * T(1, 2);
    v[5] = 2.0 * T(0, 2);
}

// Polar decomposition F = R * U via spectral decomposition of C = F^T * F.
// Computes R, U, and log(U) (left as a 3x3 symmetric tensor) in one pass —
// we already have the eigenpairs of C, so deriving U^{1/2}, U^{-1/2} and
// log(U) is just a swap of the diagonal weights.
//
// Returns true on convergence of the eigensolver. On false, R and logU are
// undefined; the caller should fall back to identity rotation + zero strain.
bool PolarDecomposeAndLog(const Matrix& F,
                          Matrix& R,
                          Matrix& logU)
{
    Matrix C(3, 3);
    noalias(C) = prod(trans(F), F);

    Matrix V(3, 3);          // GaussSeidelEigenSystem returns V with eigenvectors as columns:
    Matrix Lambda(3, 3);     // C = V * Lambda * V^T  (Lambda diagonal).
    const bool converged = MathUtils<double>::GaussSeidelEigenSystem(C, V, Lambda, 1.0e-18, 100);
    if (!converged) return false;

    // Build diagonal weights for sqrt(Lambda), 1/sqrt(Lambda), and 0.5*log(Lambda).
    // BDBtProductOperation(A, D, B) computes A = B * D * B^T, so we feed it the
    // weighted diagonal in D and V in B to get the corresponding tensor function.
    Matrix Dsqrt    = ZeroMatrix(3, 3);
    Matrix DsqrtInv = ZeroMatrix(3, 3);
    Matrix DlogU    = ZeroMatrix(3, 3);
    for (int i = 0; i < 3; ++i) {
        const double lam = Lambda(i, i);
        if (lam <= 0.0) return false;
        const double s = std::sqrt(lam);
        Dsqrt   (i, i) = s;
        DsqrtInv(i, i) = 1.0 / s;
        DlogU   (i, i) = 0.5 * std::log(lam);          // log(sqrt(lam)) = 0.5 * log(lam)
    }

    Matrix U_inv(3, 3);
    Matrix U(3, 3);
    MathUtils<double>::BDBtProductOperation(U,     Dsqrt,    V);      // U     = V * sqrt(Lambda) * V^T
    MathUtils<double>::BDBtProductOperation(U_inv, DsqrtInv, V);      // U^-1
    MathUtils<double>::BDBtProductOperation(logU,  DlogU,    V);      // log(U)

    // R = F * U^-1
    R.resize(3, 3, false);
    noalias(R) = prod(F, U_inv);
    return true;
}

// Apply R · T · R^T to a 6-vector that stores a 3x3 SYMMETRIC TENSOR in
// Kratos Voigt order [11, 22, 33, 12, 23, 13] (NO engineering-shear factor;
// off-diagonals are tensor components, not 2*tensor). Used for both stress
// and intergranular-strain pre-rotation in the Abaqus convention.
void RotateSymTensorVoigtKratos(const Matrix& R, double* v6)
{
    Matrix T(3, 3);
    T(0, 0) = v6[0]; T(1, 1) = v6[1]; T(2, 2) = v6[2];
    T(0, 1) = T(1, 0) = v6[3];
    T(1, 2) = T(2, 1) = v6[4];
    T(0, 2) = T(2, 0) = v6[5];

    Matrix RT(3, 3);
    noalias(RT) = prod(R, Matrix(prod(T, trans(R))));

    v6[0] = RT(0, 0); v6[1] = RT(1, 1); v6[2] = RT(2, 2);
    v6[3] = RT(0, 1);
    v6[4] = RT(1, 2);
    v6[5] = RT(0, 2);
}

// Linear-elastic Cauchy step that mirrors the Fortran calc_elasti_h fallback
// triggered by inittension != 0 (E = 100 kPa, nu = 0.48). Updates `sigma`
// in place by adding D : Delta_eps; D is the isotropic 6x6 tangent.
//
// Strain Voigt input is engineering shear; stress output is plain Voigt.
void ApplyElasticIncrement(Vector& sigma_kratos,
                           const Vector& dstran_kratos,
                           Matrix& C_kratos)
{
    constexpr double E  = 100.0;
    constexpr double nu = 0.48;
    const double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double mu  = E / (2.0 * (1.0 + nu));
    const double G2  = 2.0 * mu;

    C_kratos = ZeroMatrix(6, 6);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) C_kratos(i, j) = lam;
        C_kratos(i, i) += G2;
    }
    // Engineering-shear blocks: D_shear = mu (NOT 2*mu) because dstran has
    // the factor of 2 baked in already.
    C_kratos(3, 3) = mu;
    C_kratos(4, 4) = mu;
    C_kratos(5, 5) = mu;

    Vector dsigma = prod(C_kratos, dstran_kratos);
    sigma_kratos += dsigma;
}

} // unnamed namespace

void SandHypoplasticFortranDllLaw::CalculateMaterialResponseKirchhoff(Parameters& rValues)
{
    EnsureUmatLoaded(rValues.GetMaterialProperties());

    Flags& r_options = rValues.GetOptions();
    const Properties& r_props = rValues.GetMaterialProperties();
    const Vector& r_params = r_props[MATERIAL_PARAMETERS];

    // ---- 1. Pull F_{n+1} and recover F_n from inverse cache ---------------
    const Matrix& F_new = rValues.GetDeformationGradientF();
    const double  J_new = rValues.GetDeterminantF();

    // Guard against non-physical / numerically broken kinematics. NaN, inf,
    // and J <= 0 must hard-fail BEFORE polar decomposition. Polar-decomp
    // alone cannot catch an inverted F because F^T·F is positive-definite
    // for any nonsingular F (proper or improper); the eigensolver succeeds
    // and produces an improper-rotation R with det(R) = -1, which would
    // then be silently used to "co-rotate" stress. Codex adversarial review
    // 2026-04-29 round 9 Finding 1.
    KRATOS_ERROR_IF_NOT(std::isfinite(J_new))
        << "SandHypoplasticFortranDllLaw: det(F_new) is non-finite ("
        << J_new << "). Upstream element kinematics produced NaN/inf — the "
           "constitutive update cannot proceed." << std::endl;
    KRATOS_ERROR_IF(J_new <= 0.0)
        << "SandHypoplasticFortranDllLaw: det(F_new) = " << J_new
        << " is non-positive. The deformation gradient is inverted or singular; "
           "no physical Cauchy stress can be defined for this configuration. "
           "Reduce the time-step / strain-increment magnitude or investigate "
           "the upstream element." << std::endl;

    Matrix F_old(3, 3);
    if (mInverseDeformationGradientF0.size1() == 3) {
        double det_F_old_inv;
        Matrix tmp(3, 3);
        MathUtils<double>::InvertMatrix3(mInverseDeformationGradientF0, tmp, det_F_old_inv);
        noalias(F_old) = tmp;
    } else {
        noalias(F_old) = IdentityMatrix(3);
    }

    // ---- 2. Incremental F and polar decomposition --------------------------
    Matrix F_inc(3, 3);
    if (mInverseDeformationGradientF0.size1() == 3) {
        noalias(F_inc) = prod(F_new, mInverseDeformationGradientF0);
    } else {
        noalias(F_inc) = F_new;
    }

    // Same guard for the increment determinant (catches incremental inversion
    // even when both F_new and F_old are individually fine).
    const double J_inc = MathUtils<double>::Det3(F_inc);
    KRATOS_ERROR_IF_NOT(std::isfinite(J_inc))
        << "SandHypoplasticFortranDllLaw: det(F_inc) is non-finite ("
        << J_inc << "). Likely an inverse-of-singular-matrix in the upstream "
           "F-history cache." << std::endl;
    KRATOS_ERROR_IF(J_inc <= 0.0)
        << "SandHypoplasticFortranDllLaw: det(F_inc) = " << J_inc
        << " is non-positive. The incremental deformation gradient is inverted "
           "(F-history vs current F crossed an inversion boundary). Cannot "
           "proceed without producing an improper rotation in the polar "
           "decomposition." << std::endl;

    Matrix R_inc(3, 3);
    Matrix logU_inc(3, 3);
    const bool polar_ok = PolarDecomposeAndLog(F_inc, R_inc, logU_inc);
    KRATOS_ERROR_IF_NOT(polar_ok)
        << "SandHypoplasticFortranDllLaw: polar decomposition of F_inc failed "
           "(non-positive eigenvalue, NaN, or singular). The deformation "
           "gradient is degenerate and the constitutive update cannot proceed. "
           "Diagnostics: det(F_new) = " << J_new
        << ", det(F_inc) = " << MathUtils<double>::Det3(F_inc)
        << ", trace(F_inc^T F_inc) = "
        << (F_inc(0,0)*F_inc(0,0) + F_inc(1,0)*F_inc(1,0) + F_inc(2,0)*F_inc(2,0)
          + F_inc(0,1)*F_inc(0,1) + F_inc(1,1)*F_inc(1,1) + F_inc(2,1)*F_inc(2,1)
          + F_inc(0,2)*F_inc(0,2) + F_inc(1,2)*F_inc(1,2) + F_inc(2,2)*F_inc(2,2))
        << ". Reduce the time-step / strain-increment magnitude or investigate "
           "the upstream element kinematics." << std::endl;

    // ---- 3. DSTRAN = log(V_inc) in Kratos engineering-shear Voigt ----------
    // Both stress (sigma_old) and intergranular strain are pre-rotated by
    // R_inc into the new (rotated) lab frame before being passed to UMAT.
    // The strain increment must be in the SAME frame for the UMAT's rate-form
    // constitutive law to make sense. logU_inc is the right-stretch log
    // (in F_inc's "old" frame); to bring it into the new lab frame we apply
    // log(V) = R · log(U) · R^T (left-stretch log, by similarity).
    //
    // Codex adversarial review 2026-04-29 round 11 Finding 1: previously
    // logU_inc was passed as-is, leaving DSTRAN in the unrotated frame
    // while stress/state were rotated. Pure-stretch (R=I) and pure-rotation
    // (logU=0) tests both miss this — the bug only fires for combined R·U.
    // Verified by tests/path_a_results/test_combined_rotation_stretch.py.
    Matrix logV_inc(3, 3);
    noalias(logV_inc) = prod(R_inc, Matrix(prod(logU_inc, trans(R_inc))));
    Vector dstran_kratos(6);
    StrainTensorToKratosVoigt(logV_inc, dstran_kratos);

    // ---- 4. Pre-rotate previous Cauchy by R_inc -> corotational frame ------
    Matrix sigma_old_lab(3, 3);
    KratosVoigtToStressTensor(mStressVectorCauchyFinalized, sigma_old_lab);
    Matrix sigma_old_corot(3, 3);
    noalias(sigma_old_corot) = prod(R_inc, Matrix(prod(sigma_old_lab, trans(R_inc))));
    Vector sigma_corot_kratos(6);
    StressTensorToKratosVoigt(sigma_old_corot, sigma_corot_kratos);

    // ---- 5. Tensile-elastic fallback (mirrors inittension+calc_elasti_h) ---
    // Mirrors Fortran check_RKF_h verbatim. Inputs to the Fortran kernel are
    // tension-positive (Abaqus convention), but check_RKF_h works with the
    // shifted stress sig_star = sigma - p_t*I and then flips sign so its
    // internal `pmean` is COMPRESSION-positive of sig_star. Translation:
    //   pmean_compplus(sig_star) = -trace(sigma)/3 + p_t
    //   tmin_compplus(sig_star)  = p_t - max(principal_i of sigma)
    // Fallback if either drops to or below minstress = p_t/4. (For purely
    // compressive states with sigma in the -100 kPa range the conditions are
    // never close, so this branch only fires near genuinely tensile states.)
    bool tensile_fallback = false;
    Vector eig_principal(3);
    {
        Matrix eig_vec(3, 3), eig_val(3, 3);
        if (MathUtils<double>::GaussSeidelEigenSystem(sigma_old_corot, eig_vec, eig_val, 1.0e-18, 100)) {
            for (int i = 0; i < 3; ++i) eig_principal[i] = eig_val(i, i);
        } else {
            // Numerically degenerate symmetric 3x3 — punt to elastic.
            tensile_fallback = true;
            for (int i = 0; i < 3; ++i) eig_principal[i] = 0.0;
        }
    }
    const double p_t = (r_params.size() >= 2) ? r_params[1] : 0.0;
    const double minstress = p_t / 4.0;
    const double trace_sigma = eig_principal[0] + eig_principal[1] + eig_principal[2];
    const double pmean_compplus = -trace_sigma / 3.0 + p_t;
    const double max_principal = std::max({eig_principal[0],
                                           eig_principal[1],
                                           eig_principal[2]});
    const double tmin_compplus = p_t - max_principal;
    if (pmean_compplus <= minstress) tensile_fallback = true;
    if (tmin_compplus  <= minstress) tensile_fallback = true;

    Vector sigma_new_kratos(6);   // working Cauchy in Kratos Voigt order, corotational frame
    Matrix C_kratos(6, 6);

    if (tensile_fallback) {
        // Linear-elastic stand-in for the Fortran inittension+calc_elasti_h
        // branch. The kernel itself doesn't update statev under inittension,
        // so neither do we. BUT: the intergranular-strain slots [0..5] are a
        // SECOND-ORDER TENSOR in lab coordinates, and we still need to
        // co-rotate it by R_inc — otherwise a particle that enters fallback
        // mid-rotation keeps its IS history in stale axes (Codex review
        // 2026-04-29 round 6 Finding 2). Stress is already pre-rotated above
        // (sigma_old_corot), so we apply the same R_inc to the IS tensor for
        // consistency.
        sigma_new_kratos = sigma_corot_kratos;
        ApplyElasticIncrement(sigma_new_kratos, dstran_kratos, C_kratos);

        double statev_rot[15];
        for (int i = 0; i < 15; ++i) statev_rot[i] = mStateVarsFinalized[i];
        // Slots 0..5 are the IS tensor in Kratos Voigt order. Slots 6..14
        // (void ratio + postprocessing scalars) are NOT tensors — copy as-is.
        RotateSymTensorVoigtKratos(R_inc, statev_rot);
        for (int i = 0; i < 15; ++i) mStateVarsTrial[i] = statev_rot[i];
    } else {
        // ---- 6. Voigt-swap (Kratos -> Abaqus) on stress, dstran, statev ---
        double stress_abq[6];
        double dstran_abq[6];
        double stran_abq[6] = {0, 0, 0, 0, 0, 0};   // not used by hypoplastic kernel
        for (int i = 0; i < 6; ++i) {
            stress_abq[i] = sigma_corot_kratos[i];
            dstran_abq[i] = dstran_kratos[i];
        }
        std::swap(stress_abq[4], stress_abq[5]);
        std::swap(dstran_abq[4], dstran_abq[5]);

        double statev_abq[15];
        for (int i = 0; i < 15; ++i) statev_abq[i] = mStateVarsFinalized[i];

        // Intergranular strain (statev[0..5]) is a SECOND-ORDER TENSOR in lab
        // coordinates. To mirror the Abaqus pre-rotation we apply on stress,
        // we co-rotate IS by the same R_inc (R · IS · R^T) BEFORE passing to
        // UMAT. Without this, IS records the strain history in stale axes
        // and the model's hysteresis becomes non-objective under finite
        // rotation. Fixed in response to Codex adversarial review 2026-04-28
        // Finding 1 (the previous round only fixed stress co-rotation).
        // Verified by tests/path_a_results/test_is_rotation.py.
        RotateSymTensorVoigtKratos(R_inc, statev_abq);

        // Voigt-order swap idx 4<->5 on the (now rotated) intergranular
        // strain to convert Kratos order [del_11,del_22,del_33,del_12,del_23,del_13]
        // to Abaqus order [del_11,del_22,del_33,del_12,del_13,del_23].
        std::swap(statev_abq[4], statev_abq[5]);

        // ---- 7. Build the remaining UMAT arguments ------------------------
        char cmname[80];
        std::fill(cmname, cmname + 80, ' ');
        const char* tag = "SAND_HYPO";
        std::copy(tag, tag + std::strlen(tag), cmname);

        int ndi    = 3;
        int nshr   = 3;
        int ntens  = 6;
        int nstatv = 15;

        // Material parameters: copy into a plain double[] for safety
        // (Kratos Vector storage is contiguous in practice but the cast is
        // explicit here to avoid any &v.data()[0] surprises).
        const int nprops_int = static_cast<int>(r_params.size());
        std::vector<double> props_local(nprops_int);
        for (int i = 0; i < nprops_int; ++i) props_local[i] = r_params[i];
        int nprops = nprops_int;

        double sse = 0.0, spd = 0.0, scd = 0.0;
        double rpl = 0.0, drpldt = 0.0;
        double ddsddt[6] = {0, 0, 0, 0, 0, 0};
        double drplde[6] = {0, 0, 0, 0, 0, 0};
        double temp = 0.0, dtemp = 0.0;
        double predef = 0.0, dpred = 0.0;

        const ProcessInfo& r_info = rValues.GetProcessInfo();
        double dtime = r_info.Has(DELTA_TIME) ? r_info[DELTA_TIME] : 1.0;
        double time[2] = {0.0, 0.0};
        if (r_info.Has(TIME)) time[1] = r_info[TIME] - dtime;

        double coords[3] = {0, 0, 0};
        // DROT: rotation between current and previous co-rotational frames.
        // Path A handles the rotation in C++, so the UMAT receives identity here.
        double drot[9] = {1, 0, 0,  0, 1, 0,  0, 0, 1};
        double pnewdt = 1.0;
        double celent = 1.0;

        // Deformation gradient inputs in column-major (Fortran) layout.
        double dfgrd0[9];
        double dfgrd1[9];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dfgrd0[j * 3 + i] = F_old(i, j);
                dfgrd1[j * 3 + i] = F_new(i, j);
            }
        }

        int noel  = 0, npt = 0;
        int layer = 0, kspt = 0;
        // The local-fork Fortran (fortran/sand_hypo_is.for, see PATCHES.md)
        // disables the PLAXIS-mode trigger so we no longer need to dodge
        // (kstep=1, kinc=1). Pass through whatever Kratos provides.
        int kstep = r_info.Has(STEP) ? r_info[STEP] : 1;
        int kinc  = r_info.Has(NL_ITERATION_NUMBER) ? r_info[NL_ITERATION_NUMBER] : 1;

        double ddsdde_local[36] = {0};   // Fortran column-major, populated by UMAT

        // ---- 8. Call the Fortran UMAT --------------------------------------
        s_umat_ptr(stress_abq, statev_abq, ddsdde_local,
                   &sse, &spd, &scd,
                   &rpl, ddsddt, drplde, &drpldt,
                   stran_abq, dstran_abq, time, &dtime,
                   &temp, &dtemp, &predef, &dpred,
                   cmname, &ndi, &nshr, &ntens, &nstatv,
                   props_local.data(), &nprops,
                   coords, drot, &pnewdt, &celent,
                   dfgrd0, dfgrd1,
                   &noel, &npt, &layer, &kspt, &kstep, &kinc,
                   static_cast<std::size_t>(80));

        // ---- 9. Detect UMAT step rejection ---------------------------------
        // The local-fork kernel (fortran/sand_hypo_is.for, see PATCHES.md)
        // sets pnewdt = 0.25 when its RKF integrator gives up (error == 3)
        // and silently restores y_n. Without this signal the bridge would
        // commit the stale (unchanged) stress as a "successful" trial state,
        // corrupting the material history. KRATOS_ERROR here lets the calling
        // strategy catch + cut back; if the caller doesn't catch, the run
        // halts with a diagnostic instead of producing wrong physics.
        //
        // The Kratos ConstitutiveLaw API has no native pnewdt-equivalent
        // (unlike Abaqus), so a hard error is the right behavior — it makes
        // the failure loud and stops state corruption. A future enhancement
        // could replace this with a softer rValues-flag-based signal once
        // Kratos MPM grows a step-cut interface.
        KRATOS_ERROR_IF(pnewdt < 1.0)
            << "SandHypoplasticFortranDllLaw: Fortran UMAT rejected the step "
               "(pnewdt=" << pnewdt << "). RKF integration could not converge "
               "for the current strain increment. Reduce the time-step / "
               "strain-increment magnitude and retry."
            << " Diagnostic: |dstran|^2 = " << inner_prod(dstran_kratos, dstran_kratos)
            << ", trace(sigma_old) = "
            << (mStressVectorCauchyFinalized[0] + mStressVectorCauchyFinalized[1]
                + mStressVectorCauchyFinalized[2]) << "." << std::endl;

        // Finite-value hygiene on UMAT outputs (Codex round 11 Finding 2).
        // Even with pnewdt = 1.0, NaN/Inf can leak from internal degeneracies
        // (e.g. log of a near-zero void ratio, divide by p_t when p_t = 0).
        // Catch them BEFORE writing to mStressVectorCauchyTrial / mStateVarsTrial.
        for (int i = 0; i < 6; ++i) {
            KRATOS_ERROR_IF_NOT(std::isfinite(stress_abq[i]))
                << "SandHypoplasticFortranDllLaw: UMAT returned non-finite "
                   "stress[" << i << "] = " << stress_abq[i]
                << ". Trial state would be corrupted; rejecting step." << std::endl;
        }
        for (int i = 0; i < 15; ++i) {
            KRATOS_ERROR_IF_NOT(std::isfinite(statev_abq[i]))
                << "SandHypoplasticFortranDllLaw: UMAT returned non-finite "
                   "statev[" << i << "] = " << statev_abq[i] << "." << std::endl;
        }
        for (int i = 0; i < 36; ++i) {
            KRATOS_ERROR_IF_NOT(std::isfinite(ddsdde_local[i]))
                << "SandHypoplasticFortranDllLaw: UMAT returned non-finite "
                   "ddsdde[" << (i / 6) << "][" << (i % 6) << "] (column-major idx "
                << i << ") = " << ddsdde_local[i] << "." << std::endl;
        }

        // ---- 10. Decode UMAT output ---------------------------------------
        // 10a) stress: swap Abaqus -> Kratos (idx 4 <-> 5)
        std::swap(stress_abq[4], stress_abq[5]);
        for (int i = 0; i < 6; ++i) sigma_new_kratos[i] = stress_abq[i];

        // 9b) statev: swap intergranular strain back, then commit to trial.
        std::swap(statev_abq[4], statev_abq[5]);
        for (int i = 0; i < 15; ++i) mStateVarsTrial[i] = statev_abq[i];

        // 9c) DDSDDE: column-major -> row-major Matrix in Abaqus Voigt order,
        // then swap rows AND cols 4<->5 to land in Kratos Voigt order.
        Matrix C_abq(6, 6);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                C_abq(i, j) = ddsdde_local[j * 6 + i];

        auto k_swap = [](int idx) { return idx == 4 ? 5 : (idx == 5 ? 4 : idx); };
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                C_kratos(k_swap(i), k_swap(j)) = C_abq(i, j);
    }

    // ---- 10. Cauchy stress is already in the lab frame (Abaqus convention)
    // We pre-rotated sigma_old by R_inc (R · σ · R^T) so the UMAT's input was
    // already expressed in the new (n+1) lab orientation. The UMAT's output
    // therefore IS the new lab-frame Cauchy stress — no post-rotation needed.
    //
    // The earlier (R^T · σ · R) post-rotation was a bug: it was the inverse
    // of the pre-rotation, so the two cancelled and rigid-body rotation
    // produced no co-rotation of the stress (Codex adversarial review,
    // 2026-04-28, Finding 1). Convention now: pre-rotate IS R · σ · R^T,
    // post-rotate IS identity. Verified by test_rigid_rotation.py.
    Vector& sigma_lab_kratos = sigma_new_kratos;

    // Cauchy in lab frame is the trial state for Finalize.
    mStressVectorCauchyTrial = sigma_lab_kratos;

    // ---- 11. Cauchy -> Kirchhoff for output (Kratos MPM expects tau) -------
    Vector kirchhoff = J_new * sigma_lab_kratos;
    Matrix C_kirchhoff = J_new * C_kratos;

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        Vector& r_stress = rValues.GetStressVector();
        if (r_stress.size() != 6) r_stress.resize(6, false);
        noalias(r_stress) = kirchhoff;
    }

    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        Matrix& r_C = rValues.GetConstitutiveMatrix();
        if (r_C.size1() != 6 || r_C.size2() != 6) r_C.resize(6, 6, false);
        noalias(r_C) = C_kirchhoff;
    }
}

void SandHypoplasticFortranDllLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    // Compute the Kirchhoff response, then divide by J for the Cauchy slot.
    // CalculateMaterialResponseKirchhoff already KRATOS_ERRORs on
    // non-finite or non-positive det(F_new) (see Finding 1 guards above),
    // so by the time we get here J is guaranteed strictly positive and
    // finite. No silent J <= 0 fallback — invalid kinematics must surface.
    CalculateMaterialResponseKirchhoff(rValues);

    Flags& r_options = rValues.GetOptions();
    const double J = rValues.GetDeterminantF();
    const double inv_J = 1.0 / J;

    if (r_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        rValues.GetStressVector() *= inv_J;
    }
    if (r_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        rValues.GetConstitutiveMatrix() *= inv_J;
    }
}

void SandHypoplasticFortranDllLaw::CommitTrialAndRollFHistory(Parameters& rValues)
{
    // Commit the trial state produced by the latest CalculateMaterialResponse
    // (which both finalize overrides re-trigger at the finalize-time F so the
    // trial cache is guaranteed to match the F that the element commits).
    mStressVectorCauchyFinalized = mStressVectorCauchyTrial;
    mStateVarsFinalized          = mStateVarsTrial;

    const Matrix& F_new = rValues.GetDeformationGradientF();
    Matrix F_new_inv(3, 3);
    double det_F_new = 0.0;
    MathUtils<double>::InvertMatrix3(F_new, F_new_inv, det_F_new);
    mInverseDeformationGradientF0 = F_new_inv;
    mDeterminantF0 = det_F_new;
}

void SandHypoplasticFortranDllLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    // Codex round-12 [high]: MPMUpdatedLagrangian::FinalizeSolutionStep builds
    // a *fresh* GeneralVariables + ConstitutiveLaw::Parameters, calls
    // FinalizeMaterialResponse(Values, ...), and then commits Values.StressVector
    // straight into mMP.cauchy_stress_vector. The previous override only
    // committed our private trial cache and never wrote rValues.GetStressVector(),
    // so every material point's persisted stress (VTK output, restart) was
    // whatever InitializeGeneralVariables had left in the buffer. The 9 unit
    // tests passed only because they reused one Parameters instance across
    // Calculate and Finalize.
    //
    // Mirror the parent HyperElastic3DLaw three-stage pattern: set the
    // FINALIZE_MATERIAL_RESPONSE flag, re-run Calculate at the finalize-time F
    // (this populates rValues.GetStressVector() AND refreshes our trial cache
    // for the *exact* final F, even when the last NL-iteration Calculate was
    // at a slightly different trial F), then commit trial -> finalized and
    // roll the F-history. Calculate is deterministic in F + finalized state,
    // so re-calling is idempotent.
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseKirchhoff(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    CommitTrialAndRollFHistory(rValues);
}

void SandHypoplasticFortranDllLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // MPMUpdatedLagrangian invokes FinalizeMaterialResponse with
    // StressMeasure_Cauchy (mpm_updated_lagrangian.cpp:545,906), so this is
    // the override that actually drives MPM's persisted stress. Route through
    // the Cauchy variant of Calculate so rValues.GetStressVector() is the
    // Cauchy stress (not Kirchhoff).
    rValues.Set(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);
    this->CalculateMaterialResponseCauchy(rValues);
    rValues.Reset(ConstitutiveLaw::FINALIZE_MATERIAL_RESPONSE);

    CommitTrialAndRollFHistory(rValues);
}

// ----- DLL loader (copied from GMA SmallStrainUMAT3DLaw, renamed) ----------
//
// Per project_fortran_bridge_toolchain.md the gfortran-built DLL exports the
// symbol "umat_" (lowercase, trailing underscore) on BOTH Windows and Linux
// when built with default flags. GMA's loader hard-codes "umat" without
// underscore on Windows because their tooling differs; we override that for
// both branches because Path A specifically uses gfortran (PoC verified).

void SandHypoplasticFortranDllLaw::EnsureUmatLoaded(const Properties& rMaterialProperties)
{
    if (mUmatResolved && s_umat_ptr != nullptr) return;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPO_FORTRAN_DLL_PATH))
        << "SandHypoplasticFortranDllLaw: Properties[SAND_HYPO_FORTRAN_DLL_PATH] is not set. "
           "Set it to the absolute path of sand_hypo_is.dll in the materials.json (or "
           "via Properties.SetValue from Python before InitializeMaterial)." << std::endl;
    const std::string& dll_path = rMaterialProperties[SAND_HYPO_FORTRAN_DLL_PATH];

    // Finding 2 guard: the s_umat_ptr cache is process-wide, so loading a
    // SECOND distinct DLL would silently retarget every already-resolved
    // instance to the new kernel on its next call. Hard-error rather than
    // ship that footgun. See TECH_DEBT.md for the per-instance redesign
    // that would lift this restriction.
    KRATOS_ERROR_IF(s_umat_ptr != nullptr && !s_loaded_dll_path.empty()
                    && s_loaded_dll_path != dll_path)
        << "SandHypoplasticFortranDllLaw: a different sand-hypoplastic DLL is already "
           "loaded process-wide ('" << s_loaded_dll_path << "'); refusing to load '"
        << dll_path << "'. Loading multiple distinct sand-hypoplastic DLLs in one "
           "process is not supported by the current global-pointer cache. Use the "
           "same DLL for all materials, or restart the process." << std::endl;

    // Already loaded for the same path? Reuse.
    if (s_umat_ptr != nullptr && s_loaded_dll_path == dll_path) {
        mUmatResolved = true;
        return;
    }

    bool ok = false;
#ifdef KRATOS_COMPILED_IN_WINDOWS
    ok = LoadUmatWindows(dll_path);
#elif defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
    ok = LoadUmatLinux(dll_path);
#else
    KRATOS_ERROR << "SandHypoplasticFortranDllLaw: DLL loading is only supported "
                 << "on Windows and Linux. Path: " << dll_path << std::endl;
#endif

    KRATOS_ERROR_IF_NOT(ok)
        << "SandHypoplasticFortranDllLaw: failed to resolve umat_ from "
        << dll_path << std::endl;

    s_loaded_dll_path = dll_path;
    mUmatResolved     = true;
}

bool SandHypoplasticFortranDllLaw::LoadUmatWindows(const std::string& rDllPath)
{
#ifdef KRATOS_COMPILED_IN_WINDOWS
    // Try the path as given, then a .so->.dll fallback (Linux-style names in
    // a materials.json file shared between OSes).
    HMODULE h = ::LoadLibraryA(rDllPath.c_str());
    if (!h) {
        std::string alt = rDllPath;
        const auto pos = alt.find(".so");
        if (pos != std::string::npos) {
            alt.replace(pos, 3, ".dll");
            h = ::LoadLibraryA(alt.c_str());
        }
    }
    if (!h) {
        KRATOS_INFO("SandHypoplasticFortranDllLaw")
            << "LoadLibraryA failed for " << rDllPath << " (GetLastError=" << ::GetLastError() << ")" << std::endl;
        return false;
    }

    // gfortran's default name mangling — lowercase symbol with one trailing
    // underscore. Verified against poc_fortran_bridge/poc_umat.dll.
    auto fp = reinterpret_cast<SandHypoFortranUmatPtr>(
        ::GetProcAddress(h, "umat_"));
    if (!fp) {
        KRATOS_INFO("SandHypoplasticFortranDllLaw")
            << "GetProcAddress(\"umat_\") failed in " << rDllPath
            << " (GetLastError=" << ::GetLastError() << ")" << std::endl;
        return false;
    }

    s_dll_handle = static_cast<void*>(h);
    s_umat_ptr   = fp;
    return true;
#else
    (void)rDllPath;
    return false;
#endif
}

bool SandHypoplasticFortranDllLaw::LoadUmatLinux(const std::string& rDllPath)
{
#if defined(KRATOS_COMPILED_IN_LINUX) || defined(KRATOS_COMPILED_IN_OS)
    void* h = ::dlopen(rDllPath.c_str(), RTLD_LAZY);
    if (!h) {
        std::string alt = rDllPath;
        const auto pos = alt.find(".dll");
        if (pos != std::string::npos) {
            alt.replace(pos, 4, ".so");
            h = ::dlopen(alt.c_str(), RTLD_LAZY);
        }
    }
    if (!h) {
        KRATOS_INFO("SandHypoplasticFortranDllLaw")
            << "dlopen failed for " << rDllPath << ": " << ::dlerror() << std::endl;
        return false;
    }

    auto fp = reinterpret_cast<SandHypoFortranUmatPtr>(::dlsym(h, "umat_"));
    if (!fp) {
        KRATOS_INFO("SandHypoplasticFortranDllLaw")
            << "dlsym(\"umat_\") failed in " << rDllPath << ": " << ::dlerror() << std::endl;
        return false;
    }

    s_dll_handle = h;
    s_umat_ptr   = fp;
    return true;
#else
    (void)rDllPath;
    return false;
#endif
}

} // namespace Kratos

//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Hongcheng Jiao
//

// SandHypoplastic3DLaw (implementation).
// See sand_hypoplastic_3D_law.hpp for the kernel/law contract and the
// model references. Effective-Cauchy in/out at the kernel boundary;
// this file does total<->effective + slot 4<->5 swap + bulk_w tangent
// correction + Bauer void init + failure_class -> KRATOS_ERROR / pnewdt
// routing.
//
// References of the form .for:NNN point into the reference Fortran UMAT
// implementation of this model (umat_hypoplasticity by A. Niemunis,
// publicly distributed via soilmodels.info); they document algorithm
// provenance.
#include "sand_hypoplastic_3D_law.hpp"

#include <cmath>
#include <cstring>

#include "includes/checks.h"
#include "includes/mat_variables.h"
#include "utilities/math_utils.h"

// IS_EXPLICIT lives in this application's variables header (also pulled in
// transitively via sand_hypoplastic_3D_law.hpp). Kept as an explicit
// include: the explicit-mode branch keys off IS_EXPLICIT (see the
// GetLawFeatures docstring).
#include "mpm_application_variables.h"

namespace Kratos
{

namespace
{
constexpr std::size_t kVoigt = 6;
constexpr std::size_t kPropsCount = 16;

// Default kernel substep tolerance / limits. Match the Fortran umat caller
// defaults; exposing them through Properties is a possible future extension.
constexpr double kDefaultErrTol = 1.0e-3;
constexpr int    kDefaultMaxNInt = 10000;
constexpr double kDefaultDtMin   = 1.0e-17;

// Kratos uses [11,22,33,12,23,13]; the kernel (Abaqus convention) uses
// [11,22,33,12,13,23]. Slots 4 and 5 swap; the 0..3 range is identical.
inline std::size_t SwapIdx(std::size_t i)
{
    if (i == 4) return 5;
    if (i == 5) return 4;
    return i;
}

} // anonymous namespace

//********************************************************************
// Lifecycle / Clone
//********************************************************************

SandHypoplastic3DLaw::SandHypoplastic3DLaw()
    : ConstitutiveLaw(), mF0_3x3(IdentityMatrix(3, 3)) {}

SandHypoplastic3DLaw::SandHypoplastic3DLaw(const SandHypoplastic3DLaw& rOther)
    : ConstitutiveLaw(rOther),
      mSigEffAbq(rOther.mSigEffAbq),
      mStateAbq(rOther.mStateAbq),
      mPore(rOther.mPore),
      mDtsubPersistent(rOther.mDtsubPersistent),
      mInitialized(rOther.mInitialized),
      mFirstCallStateLoaded(rOther.mFirstCallStateLoaded),
      mF0_3x3(rOther.mF0_3x3),
      mEpsPrevExplicitKratos(rOther.mEpsPrevExplicitKratos)
{
}

ConstitutiveLaw::Pointer SandHypoplastic3DLaw::Clone() const
{
    return Kratos::make_shared<SandHypoplastic3DLaw>(*this);
}

//********************************************************************
// Features / metadata
//********************************************************************

void SandHypoplastic3DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Deformation_Gradient is the canonical strain measure advertised for
    // BOTH implicit and explicit MPM paths. MPMUpdatedLagrangian::Check
    // accepts Deformation_Gradient unconditionally for implicit and as
    // one of the two accepted measures for explicit. Explicit support is
    // handled by a commit-on-explicit branch at the end of
    // CalculateMaterialResponseCauchy (see the state-commit-policy
    // comment there).
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize = static_cast<SizeType>(kVoigt);
    rFeatures.mSpaceDimension = static_cast<SizeType>(WorkingSpaceDimension());
}

//********************************************************************
// Has / GetValue / SetValue
//********************************************************************

bool SandHypoplastic3DLaw::Has(const Variable<double>& rThisVariable)
{
    return false;
}

bool SandHypoplastic3DLaw::Has(const Variable<Vector>& rThisVariable)
{
    // Expose INTERNAL_VARIABLES for inspection and post-processing.
    // 15-slot layout (mirrors the Fortran UMAT statev layout):
    //   [0..5] intergranular strain delta, Kratos Voigt order (slot 4<->5
    //          swap applied from kernel-internal Abaqus order)
    //   [6]   void ratio
    //   [7]   -pore (Fortran statev(8) sign-flipped convention)
    //   [8]   dtsub_persistent
    //   [9..14] reserved / zero
    return rThisVariable == INTERNAL_VARIABLES;
}

double& SandHypoplastic3DLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    rValue = 0.0;
    return rValue;
}

Vector& SandHypoplastic3DLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        constexpr std::size_t kStatevSize = 15;
        if (rValue.size() != kStatevSize) rValue.resize(kStatevSize, false);
        // 0..5: intergranular strain in Kratos Voigt order. SwapIdx(i)
        // maps the Abaqus-internal indices to Kratos-external indices
        // (slots 4 and 5 swap; 0..3 unchanged). Convention check: a
        // pure 1-3 shear (Kratos slot 5) loads slot 5 of the output.
        for (std::size_t i = 0; i < 6; ++i) rValue[i] = mStateAbq[SwapIdx(i)];
        rValue[6]  = mStateAbq[6];     // void ratio
        rValue[7]  = -mPore;            // Fortran statev(8) sign flip
        rValue[8]  = mDtsubPersistent;
        for (std::size_t i = 9; i < kStatevSize; ++i) rValue[i] = 0.0;
        return rValue;
    }
    return rValue;
}

void SandHypoplastic3DLaw::SetValue(const Variable<double>& /*rThisVariable*/,
                                     const double& /*rValue*/,
                                     const ProcessInfo& /*rCurrentProcessInfo*/)
{
}

void SandHypoplastic3DLaw::SetValue(const Variable<Vector>& /*rThisVariable*/,
                                     const Vector& /*rValue*/,
                                     const ProcessInfo& /*rCurrentProcessInfo*/)
{
}

//********************************************************************
// InitializeMaterial
//********************************************************************

void SandHypoplastic3DLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                               const GeometryType& /*rElementGeometry*/,
                                               const Vector& /*rShapeFunctionsValues*/)
{
    // Validate the props vector once at material init so the user gets a
    // fast KRATOS_ERROR rather than a silent first-step contract violation.
    // Always run -- props validation is cheap and helps surface bad JSON
    // during model setup even on a deserialized law that already has
    // committed history (mInitialized == true).
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PARAMETERS))
        << "SandHypoplastic3DLaw requires SAND_HYPOPLASTIC_PARAMETERS in Properties.";

    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PARAMETERS];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SAND_HYPOPLASTIC_PARAMETERS must have 16 entries; got " << props_vec.size();

    // Fortran's check_parms is IEEE-blind; finiteness is enforced here so
    // bad JSON inputs surface during model setup instead of mid-solve.
    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SAND_HYPOPLASTIC_PARAMETERS[" << i << "] is non-finite (" << props_vec[i] << ").";
    }

    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SAND_HYPOPLASTIC_PARAMETERS rejected by check_parms: failed='"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    // Idempotent reset. InitializeMaterial is called once per material
    // point at simulation start, AND may be called again after
    // deserialization (restart / partition migration) depending on caller
    // lifecycle. The committed history members (mSigEffAbq, mStateAbq,
    // mPore, mDtsubPersistent, mFirstCallStateLoaded, mF0_3x3) are
    // serialized via save/load; clobbering them here would silently
    // discard restored history. Reset only when this is the FIRST init
    // -- detected via mInitialized still false. Kratos-registered
    // prototypes start with
    // mInitialized == false (default constructor), and Clone() copies the
    // prototype's default state to each material point, so the per-MP first
    // init still resets as expected.
    //
    // INITIAL_STRESS_VECTOR + Bauer void initialization is still deferred to
    // the first `CalculateMaterialResponseCauchy` call so all step-0 inputs
    // (including the element's effective Almansi strain) are in scope at
    // the same time.
    if (!mInitialized) {
        mSigEffAbq.fill(0.0);
        mStateAbq.fill(0.0);
        mPore = 0.0;
        mDtsubPersistent = 0.0;
        mF0_3x3 = IdentityMatrix(3, 3);
        mEpsPrevExplicitKratos.fill(0.0);
        mFirstCallStateLoaded = false;
        mInitialized = true;
    }
}

//********************************************************************
// ResetMaterial
//********************************************************************

void SandHypoplastic3DLaw::ResetMaterial(const Properties& rMaterialProperties,
                                          const GeometryType& /*rElementGeometry*/,
                                          const Vector& /*rShapeFunctionsValues*/)
{
    // Re-validate props on every reset -- catches the case where Reset is
    // called with a modified material (e.g. user adjusted parameters and
    // wants to reset history) and the new values are bad.
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PARAMETERS))
        << "SandHypoplastic3DLaw::ResetMaterial requires SAND_HYPOPLASTIC_PARAMETERS in Properties.";
    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PARAMETERS];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SAND_HYPOPLASTIC_PARAMETERS must have 16 entries; got " << props_vec.size();
    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SAND_HYPOPLASTIC_PARAMETERS[" << i << "] is non-finite at ResetMaterial.";
    }
    SandHypoCpp::Props16 buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SAND_HYPOPLASTIC_PARAMETERS rejected by check_parms in ResetMaterial: '"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    // Clear all committed history. mInitialized stays true -- after
    // ResetMaterial the law IS validly initialized, just with default
    // state. mFirstCallStateLoaded resets to false so the next Calculate
    // re-reads INITIAL_STRESS_VECTOR and re-runs Bauer void init.
    mSigEffAbq.fill(0.0);
    mStateAbq.fill(0.0);
    mPore = 0.0;
    mDtsubPersistent = 0.0;
    mF0_3x3 = IdentityMatrix(3, 3);
    mEpsPrevExplicitKratos.fill(0.0);
    mFirstCallStateLoaded = false;
    mInitialized = true;
}

//********************************************************************
// CalculateMaterialResponseCauchy
//********************************************************************

void SandHypoplastic3DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_ERROR_IF_NOT(mInitialized)
        << "SandHypoplastic3DLaw::CalculateMaterialResponseCauchy called before InitializeMaterial.";

    const Properties& props = rValues.GetMaterialProperties();
    const Vector& props_vec = props[SAND_HYPOPLASTIC_PARAMETERS];
    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    if (check.error != 0) {
        KRATOS_ERROR << "SandHypoplastic3DLaw: check_parms failed at solve time: '"
                     << (check.failed ? check.failed : "?") << "' bad_value="
                     << check.bad_value;
    }
    const double bulk_w = check.parms[14];

    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double dtime = process_info[DELTA_TIME];

    // Detect IS_EXPLICIT for branching.
    // - In implicit, the law defers state commit to `FinalizeMaterialResponseCauchy`
    //   (called once at end of converged step).
    // - In explicit, the element calls Calculate ONCE per step via the
    //   `CALCULATE_EXPLICIT_MP_STRESS` hook (see MPMUpdatedLagrangian::
    //   CalculateOnIntegrationPoints) and never calls Finalize, so the law
    //   must commit internal state at end of Calculate when is_explicit.
    //   The element's `FinalizeStepVariables` handles the
    //   `mMP.cauchy_stress_vector` writeback separately.
    const bool is_explicit = (process_info.Has(IS_EXPLICIT))
                                 ? process_info.GetValue(IS_EXPLICIT)
                                 : false;

    // Option flags drive what the caller wants this Calculate to publish.
    // The law must NOT touch the stress / tangent output buffers when
    // the corresponding compute flag is unset -- some callers omit
    // SetStressVector / SetConstitutiveMatrix entirely, and GetStressVector
    // / GetConstitutiveMatrix on an unset pointer KRATOS_DEBUG_ERRORs.
    const Flags& options = rValues.GetOptions();
    const bool compute_stress  = options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool compute_tangent = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // Read F here for both paths -- the implicit branch needs it
    // for polar decomposition + Almansi-from-F write-back; the explicit
    // branch only uses it as a defense-in-depth kinematic gate / restart
    // hand-off (explicit kinematics come from rValues.StrainVector,
    // not F).
    const Matrix& F_total = rValues.GetDeformationGradientF();
    KRATOS_ERROR_IF(F_total.size1() != 3 || F_total.size2() != 3)
        << "DeformationGradientF must be 3x3; got " << F_total.size1() << "x" << F_total.size2();

    // Physical-kinematics gate. det(F) <= 0 means
    // self-inverted / collapsed material -- not a valid input to a continuum
    // constitutive update. Catch BEFORE the 3x3 inversion in
    // CalculateAlmansiFromF / before forming F_rel for polar decomposition,
    // because both inversion paths silently produce NaN/Inf on a singular F
    // and the polar decomp can return an improper rotation (det R = -1) for
    // orientation-reversing F (C = F^T F is SPD regardless of sign(det F),
    // so the eigenvalue-positivity check inside PolarDecompFAndLogU does
    // NOT block det(F) < 0). On the explicit path F is typically I or
    // I+symmetric_strain_increment (see mpm_explicit_utilities.cpp);
    // keeping the gate is cheap and rejects upstream contract breakage.
    const double det_F_total = MathUtils<double>::Det(F_total);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_total) && det_F_total > 0.0)
        << "SandHypoplastic3DLaw: det(F_total) = " << det_F_total
        << " is non-finite or non-positive (singular / orientation-reversing kinematics).";

    // ----- First-call state load (one-shot) ----------------------------
    // Order matters. Steps mirror Fortran umat .for:165-209.
    if (!mFirstCallStateLoaded) {
        std::array<double, kVoigt> sig_total_kratos{};
        if (props.Has(INITIAL_STRESS_VECTOR)) {
            const Vector& s0 = props[INITIAL_STRESS_VECTOR];
            KRATOS_ERROR_IF(s0.size() != kVoigt)
                << "INITIAL_STRESS_VECTOR must have 6 entries; got " << s0.size();
            for (std::size_t i = 0; i < kVoigt; ++i) sig_total_kratos[i] = s0[i];
        }

        // Slot 4<->5 swap: Kratos -> Abaqus.
        std::array<double, kVoigt> sig_total_abq{};
        for (std::size_t i = 0; i < kVoigt; ++i) sig_total_abq[i] = sig_total_kratos[SwapIdx(i)];

        const double init_void = InitVoidRatioBauer(props_vec, sig_total_abq);
        mStateAbq.fill(0.0);
        mStateAbq[6] = init_void;

        mPore = 0.0; // pore pressure is not yet exposed via a Kratos variable.

        // total -> effective (in-place on the normal slots)
        std::array<double, kVoigt> sig_eff_abq = sig_total_abq;
        TotalToEffectiveAbq(sig_eff_abq, mPore);
        mSigEffAbq = sig_eff_abq;

        // mF0_3x3 is initialized to Identity at construction time; the
        // first step's F_rel therefore equals the current total F.

        mFirstCallStateLoaded = true;
    }

    // ----- Kernel-input pack -- branches on IS_EXPLICIT ----------------
    //
    // IMPLICIT path: F_rel = F_total · F0^-1, polar
    // decomp F_rel = R · U. R pre-rotates stored σ_n and δ_n into the
    // end-of-step spatial frame (Hughes-Winget), log(V) is the
    // corotational strain increment. Kratos MPM does NO host-side stress
    // rotation on the implicit path (see MPMUpdatedLagrangian);
    // objectivity is owned by the law. Almansi-from-F is written back to
    // rValues.StrainVector for MPM's `mMP.almansi_strain_vector`
    // post-processing persistence.
    //
    // EXPLICIT path: the element's
    // CalculateExplicitKinematics populates rValues.StrainVector with the
    // Jaumann-corrected cumulative Almansi
    // (see MPMExplicitUtilities::CalculateExplicitKinematics) -- THAT is
    // the only strain signal carrying rotation in explicit MPM, because
    // the F handed in via SetDeformationGradientF is either I or
    // I+symmetric_strain_increment and carries no spin.
    // The kernel-driving Δε is therefore the difference between the
    // current rValues.StrainVector and the persisted mEpsPrevExplicitKratos
    // -- NO polar decomposition, NO R-rotation of σ_n / δ_n (the Jaumann
    // history already accounts for spin in the strain stream). The strain
    // write-back is skipped so the next step's CalculateExplicitKinematics
    // can read the value it stored, apply its own Jaumann correction, and
    // produce the next strain consistently.

    std::array<double, kVoigt> sig_n_abq{};
    std::array<double, SandHypoCpp::kStateDim> q_n{};
    std::array<double, kVoigt> deps_abq{};

    if (!is_explicit) {
        // ----- IMPLICIT branch -----------------------------------------
        std::array<double, kVoigt> strain_now_kratos{};
        CalculateAlmansiFromF(F_total, strain_now_kratos);

        Matrix F0_inv(3, 3);
        double det_F0 = 0.0;
        MathUtils<double>::InvertMatrix(mF0_3x3, F0_inv, det_F0);
        KRATOS_ERROR_IF_NOT(std::isfinite(det_F0) && det_F0 > 0.0)
            << "SandHypoplastic3DLaw: persisted mF0_3x3 has det = " << det_F0
            << " (non-finite or non-positive). Internal invariant violation -- previous "
            << "step committed a corrupt F_total.";
        Matrix F_rel(3, 3);
        noalias(F_rel) = prod(F_total, F0_inv);

        // Second kinematic gate: det(F_rel) > 0 (orientation-preserving relative step).
        // Note that C = F_rel^T F_rel is SPD for any non-singular F_rel, so
        // PolarDecompFAndLogU alone cannot reject this case -- it would return
        // an improper rotation (det R = -1).
        const double det_F_rel = MathUtils<double>::Det(F_rel);
        KRATOS_ERROR_IF_NOT(std::isfinite(det_F_rel) && det_F_rel > 0.0)
            << "SandHypoplastic3DLaw: det(F_rel) = " << det_F_rel
            << " is non-finite or non-positive (singular / orientation-reversing relative step).";

        Matrix R(3, 3);
        PolarDecompFAndLogU(F_rel, R, deps_abq);

        // Pre-rotate stored stress and intergranular strain by R.
        RotateStressVoigtAbq(mSigEffAbq, R, sig_n_abq);

        std::array<double, kVoigt> delta_n_abq{};
        for (std::size_t i = 0; i < kVoigt; ++i) delta_n_abq[i] = mStateAbq[i];
        std::array<double, kVoigt> delta_n_abq_rotated{};
        RotateStrainVoigtAbq(delta_n_abq, R, delta_n_abq_rotated);
        for (std::size_t i = 0; i < kVoigt; ++i) q_n[i] = delta_n_abq_rotated[i];
        q_n[6] = mStateAbq[6]; // void ratio scalar -- no rotation

        // Write Almansi back so MPM's FinalizeStepVariables persists it into
        // `mMP.almansi_strain_vector` for post-processing.
        if (rValues.IsSetStrainVector()) {
            Vector& strain_out_kratos = rValues.GetStrainVector();
            if (strain_out_kratos.size() != kVoigt) strain_out_kratos.resize(kVoigt, false);
            for (std::size_t i = 0; i < kVoigt; ++i) strain_out_kratos[i] = strain_now_kratos[i];
        }
    } else {
        // ----- EXPLICIT branch -----------------------------------------
        //
        // Consume the element-populated Jaumann-corrected cumulative
        // Almansi from rValues.GetStrainVector(). Differentiate against
        // mEpsPrevExplicitKratos to get the kernel-driving Δε. No
        // polar decomp, no R-rotation of σ_n/δ_n.
        KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector())
            << "SandHypoplastic3DLaw: explicit branch requires the element to "
               "set the StrainVector before calling Calculate. MPM normally "
               "populates rValues.GetStrainVector() from "
               "mMP.almansi_strain_vector via SetStrainVector. If you are "
               "calling Calculate manually, attach a Vector(6) via "
               "Parameters::SetStrainVector first.";
        const Vector& strain_now_vec = rValues.GetStrainVector();
        KRATOS_ERROR_IF(strain_now_vec.size() != kVoigt)
            << "SandHypoplastic3DLaw: explicit branch expects StrainVector of "
               "size 6 (Kratos Voigt); got " << strain_now_vec.size() << ".";

        std::array<double, kVoigt> strain_now_kratos{};
        for (std::size_t i = 0; i < kVoigt; ++i) {
            const double v = strain_now_vec[i];
            KRATOS_ERROR_IF_NOT(std::isfinite(v) && std::abs(v) <= 1.0e30)
                << "SandHypoplastic3DLaw: explicit StrainVector[" << i << "] = "
                << v << " is non-finite or > 1e30. CalculateExplicitKinematics "
                   "produced a corrupt strain.";
            strain_now_kratos[i] = v;
        }

        // Δε in Kratos Voigt order (engineering shear). Off-diagonals are
        // doubled the same way on both terms so the difference preserves
        // engineering-shear convention -- the kernel and the Fortran umat
        // both expect engineering shear (kernel norm_D guard mirrors
        // .for:201-209 `umatisnan_h(norm_D)`).
        std::array<double, kVoigt> deps_kratos{};
        for (std::size_t i = 0; i < kVoigt; ++i) {
            deps_kratos[i] = strain_now_kratos[i] - mEpsPrevExplicitKratos[i];
            KRATOS_ERROR_IF_NOT(std::isfinite(deps_kratos[i])
                                && std::abs(deps_kratos[i]) <= 1.0e30)
                << "SandHypoplastic3DLaw: explicit deps_kratos[" << i << "] = "
                << deps_kratos[i] << " is non-finite or > 1e30. Likely a "
                   "discontinuity in the Jaumann strain stream "
                   "(mEpsPrevExplicitKratos out of sync after restart, or "
                   "element strain history reset mid-simulation).";
        }

        // Slot 4<->5 swap: Kratos -> Abaqus (the Voigt convention the kernel
        // expects; same as the implicit path's log(V) output).
        for (std::size_t i = 0; i < kVoigt; ++i) deps_abq[i] = deps_kratos[SwapIdx(i)];

        // No rotation of σ_n / δ_n -- the Jaumann history baked into Δε
        // already accounts for spin.
        sig_n_abq = mSigEffAbq;
        for (std::size_t i = 0; i < kVoigt; ++i) q_n[i] = mStateAbq[i];
        q_n[6] = mStateAbq[6];

        // Do NOT write back to rValues.StrainVector in explicit. The next
        // step's CalculateExplicitKinematics READS this slot to apply the
        // Jaumann spin correction (see MPMExplicitUtilities);
        // overwriting it would corrupt the strain stream.
    }

    const double dtsub_in = mDtsubPersistent; // first call: 0; integrate_step normalizes -> dtime.

    // ----- Call kernel -------------------------------------------------
    const auto step = SandHypoCpp::integrate_step(sig_n_abq.data(),
                                                  q_n.data(),
                                                  deps_abq.data(),
                                                  check.parms,
                                                  dtime,
                                                  dtsub_in,
                                                  kDefaultErrTol,
                                                  kDefaultMaxNInt,
                                                  kDefaultDtMin);

    // ----- Failure routing ---------------------------------------------
    using SandHypoCpp::FailureClass;
    if (step.failure_class == FailureClass::ContractViolation) {
        KRATOS_ERROR << "SandHypoplastic3DLaw: kernel ContractViolation (error=" << step.error
                     << "); inputs / persisted state corrupt. Not retryable. dtime=" << dtime
                     << ", dtsub_in=" << dtsub_in << ".";
    }
    if (step.failure_class == FailureClass::Fatal) {
        KRATOS_ERROR << "SandHypoplastic3DLaw: kernel Fatal (error=" << step.error
                     << ", parms-only failure inside get_tan / get_F_sig_q).";
    }
    if (step.failure_class == FailureClass::RkfReject) {
        // Caller-driven cutback. Absent a wired Kratos pnewdt channel,
        // the law raises a KRATOS_ERROR whose message includes the
        // literal "pnewdt=0.25" substring, mirroring the Abaqus UMAT
        // cutback convention so external drivers can detect a retryable
        // RKF reject and reduce the step. Wiring a real pnewdt-style
        // cutback channel is a possible future improvement.
        KRATOS_ERROR << "SandHypoplastic3DLaw: kernel RkfReject (error=" << step.error
                     << "); pnewdt=0.25 outer cutback required. dtime=" << dtime
                     << " too large for current state.";
    }

    // ----- Pore evolution: pore_next = pore - bulk_w * tr(deps) -------
    //
    // Mirrors Fortran solout_h .for:1856 in-place rebind. `bulk_w * tr(deps)`
    // uses the same deps that drove the kernel (pre-swap doesn't matter for
    // the trace -- normal slots [0..2] are the same in either order).
    //
    // The pore update happens AFTER integrate_step
    // returns, so kernel guards don't cover it. A typo-scale but finite
    // bulk_w (e.g. 1e29 instead of 2.2e9) times a small volumetric strain
    // can produce a non-finite or overflow pore_next, which then corrupts
    // the total-stress conversion and persists across steps. Verify
    // pore_next finite and bounded before propagating.
    const double tr_deps = deps_abq[0] + deps_abq[1] + deps_abq[2];
    const double pore_next = mPore - bulk_w * tr_deps;
    KRATOS_ERROR_IF_NOT(std::isfinite(pore_next) && std::abs(pore_next) <= 1.0e30)
        << "SandHypoplastic3DLaw: pore_next = " << pore_next
        << " is non-finite or > 1e30 (overflow). bulk_w=" << bulk_w
        << ", tr_deps=" << tr_deps << ", mPore=" << mPore << ".";

    // ----- Convert effective -> total on the published output ----------
    //
    // Note: pore_next, NOT mPore -- the bulk_w*tr(deps) correction must
    // propagate through here. Mirrors solout_h .for:1856
    // followed by .for:1862.
    std::array<double, kVoigt> sig_total_next_abq;
    std::memcpy(sig_total_next_abq.data(), step.sig, sizeof(sig_total_next_abq));
    EffectiveToTotalAbq(sig_total_next_abq, pore_next);

    // Verify the post-conversion total stress is finite and bounded
    // before publishing.
    for (std::size_t i = 0; i < kVoigt; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(sig_total_next_abq[i])
                            && std::abs(sig_total_next_abq[i]) <= 1.0e30)
            << "SandHypoplastic3DLaw: sig_total_next[" << i << "] = "
            << sig_total_next_abq[i] << " is non-finite or > 1e30.";
    }

    // ----- Abaqus -> Kratos swap, option-guarded -----------------------
    if (compute_stress) {
        Vector& stress_out_kratos = rValues.GetStressVector();
        if (stress_out_kratos.size() != kVoigt) stress_out_kratos.resize(kVoigt, false);
        SwapVectorAbqToKratos(sig_total_next_abq, stress_out_kratos);
    }

    if (compute_tangent) {
        Matrix& tangent_out_kratos = rValues.GetConstitutiveMatrix();
        if (tangent_out_kratos.size1() != kVoigt || tangent_out_kratos.size2() != kVoigt) {
            tangent_out_kratos.resize(kVoigt, kVoigt, false);
        }
        // Tangent: kernel returns 6x6 in Abaqus order. Swap rows AND columns.
        std::array<std::array<double, kVoigt>, kVoigt> D_abq{};
        for (std::size_t i = 0; i < kVoigt; ++i)
            for (std::size_t j = 0; j < kVoigt; ++j)
                D_abq[i][j] = step.D[i][j];
        SwapTangentAbqToKratos(D_abq, tangent_out_kratos);
        AddBulkWNormalBlock(tangent_out_kratos, bulk_w);
    }


    // ----- State commit policy ----------------------------------------
    //
    // IMPLICIT path: persisted state is NOT mutated here. The element may
    // call Calculate multiple times during NR iterations with different
    // trial F; stress / tangent are published but the commit is deferred to
    // FinalizeMaterialResponseCauchy (called once at end of converged step).
    //
    // EXPLICIT path: the element calls Calculate ONCE per step via the
    // `CALCULATE_EXPLICIT_MP_STRESS` hook and NEVER calls Finalize
    // (MPMUpdatedLagrangian::FinalizeSolutionStep errors out
    // for explicit). So the law MUST commit internal state here, at the
    // end of Calculate, otherwise mSigEffAbq / mStateAbq / mPore /
    // mDtsubPersistent / mEpsPrevExplicitKratos / mF0_3x3 stay at initial
    // values forever and the law is effectively stateless across
    // explicit steps. The element's `FinalizeStepVariables` (called
    // immediately after Calculate returns) handles `mMP.cauchy_stress_vector`
    // writeback separately.
    //
    // mEpsPrevExplicitKratos commits to the current step's Jaumann-corrected
    // strain so the next step's Δε is well-formed.
    //
    // mF0_3x3 is committed to F_total even though the explicit branch
    // does not use it for kinematics -- a serialized state round-tripped
    // from explicit to implicit (restart) would otherwise see an
    // out-of-date F0 on the first implicit step. Cost is one Matrix copy
    // per step; the alternative (leaving F0 stale) is a silent restart
    // hazard.
    if (is_explicit) {
        std::memcpy(mSigEffAbq.data(), step.sig, sizeof(double) * kVoigt);
        for (std::size_t i = 0; i < SandHypoCpp::kStateDim; ++i) mStateAbq[i] = step.q[i];
        mPore = pore_next;
        mDtsubPersistent = step.dtsub_next;
        mF0_3x3 = F_total;

        // Persist the begin-of-(next-)step baseline for the
        // explicit Δε differencing scheme. Re-read from rValues since the
        // local copy `strain_now_kratos` lives inside the !is_explicit
        // scope only.
        const Vector& strain_now_vec = rValues.GetStrainVector();
        for (std::size_t i = 0; i < kVoigt; ++i) mEpsPrevExplicitKratos[i] = strain_now_vec[i];
    }
}

//********************************************************************
// FinalizeMaterialResponseCauchy
//********************************************************************

void SandHypoplastic3DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // Re-run the integration with the same inputs to recover scratch state,
    // then commit to the persisted members. This deliberately mirrors the
    // Calculate path -- if anything diverges, the kernel's own determinism
    // (no RNG / no internal time-dependent caching) guarantees the second
    // run produces the same StepResult.
    KRATOS_ERROR_IF_NOT(mInitialized)
        << "SandHypoplastic3DLaw::FinalizeMaterialResponseCauchy called before InitializeMaterial.";

    const Properties& props = rValues.GetMaterialProperties();
    const Vector& props_vec = props[SAND_HYPOPLASTIC_PARAMETERS];
    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SandHypoplastic3DLaw::Finalize: check_parms failed: '"
        << (check.failed ? check.failed : "?") << "'.";
    const double bulk_w = check.parms[14];

    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double dtime = process_info[DELTA_TIME];

    // Same Almansi-from-F path as Calculate -- the element's
    // CalculateKinematics is called again on the Finalize side and again
    // doesn't populate rVariables.StrainVector. Compute Almansi from F
    // here too.
    const Matrix& F_total = rValues.GetDeformationGradientF();
    KRATOS_ERROR_IF(F_total.size1() != 3 || F_total.size2() != 3)
        << "DeformationGradientF must be 3x3 in Finalize; got " << F_total.size1() << "x" << F_total.size2();

    // Same kinematic gate as Calculate.
    const double det_F_total = MathUtils<double>::Det(F_total);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_total) && det_F_total > 0.0)
        << "SandHypoplastic3DLaw::Finalize: det(F_total) = " << det_F_total
        << " is non-finite or non-positive.";

    std::array<double, kVoigt> strain_now_kratos{};
    CalculateAlmansiFromF(F_total, strain_now_kratos);

    // First-call init MUST also happen at finalize time if Calculate*
    // was never invoked first (defensive; the standard Kratos call order
    // is Calculate* then Finalize*).
    if (!mFirstCallStateLoaded) {
        std::array<double, kVoigt> sig_total_kratos{};
        if (props.Has(INITIAL_STRESS_VECTOR)) {
            const Vector& s0 = props[INITIAL_STRESS_VECTOR];
            KRATOS_ERROR_IF(s0.size() != kVoigt)
                << "INITIAL_STRESS_VECTOR must have 6 entries; got " << s0.size();
            for (std::size_t i = 0; i < kVoigt; ++i) sig_total_kratos[i] = s0[i];
        }
        std::array<double, kVoigt> sig_total_abq{};
        for (std::size_t i = 0; i < kVoigt; ++i) sig_total_abq[i] = sig_total_kratos[SwapIdx(i)];
        const double init_void = InitVoidRatioBauer(props_vec, sig_total_abq);
        mStateAbq.fill(0.0);
        mStateAbq[6] = init_void;
        mPore = 0.0;
        std::array<double, kVoigt> sig_eff_abq = sig_total_abq;
        TotalToEffectiveAbq(sig_eff_abq, mPore);
        mSigEffAbq = sig_eff_abq;
        // mF0_3x3 already Identity from construction.
        mFirstCallStateLoaded = true;
    }

    // ----- Same corotational kinematics as Calculate. The kernel is
    // deterministic and Calculate did not mutate persisted state, so this
    // re-derives R and log(U) with identical inputs and identical results.
    Matrix F0_inv(3, 3);
    double det_F0 = 0.0;
    MathUtils<double>::InvertMatrix(mF0_3x3, F0_inv, det_F0);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F0) && det_F0 > 0.0)
        << "SandHypoplastic3DLaw::Finalize: persisted mF0_3x3 has det = " << det_F0
        << " (invariant violation).";
    Matrix F_rel(3, 3);
    noalias(F_rel) = prod(F_total, F0_inv);
    const double det_F_rel = MathUtils<double>::Det(F_rel);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_rel) && det_F_rel > 0.0)
        << "SandHypoplastic3DLaw::Finalize: det(F_rel) = " << det_F_rel
        << " (singular / orientation-reversing).";

    Matrix R(3, 3);
    std::array<double, kVoigt> deps_abq{};
    PolarDecompFAndLogU(F_rel, R, deps_abq);

    std::array<double, kVoigt> sig_n_abq_rotated{};
    RotateStressVoigtAbq(mSigEffAbq, R, sig_n_abq_rotated);

    std::array<double, kVoigt> delta_n_abq{};
    for (std::size_t i = 0; i < kVoigt; ++i) delta_n_abq[i] = mStateAbq[i];
    std::array<double, kVoigt> delta_n_abq_rotated{};
    RotateStrainVoigtAbq(delta_n_abq, R, delta_n_abq_rotated);

    std::array<double, SandHypoCpp::kStateDim> q_n_rotated{};
    for (std::size_t i = 0; i < kVoigt; ++i) q_n_rotated[i] = delta_n_abq_rotated[i];
    q_n_rotated[6] = mStateAbq[6];

    // Same strain write-back as Calculate: persist Almansi so MPM's
    // FinalizeStepVariables gets a meaningful value into
    // mMP.almansi_strain_vector for post-processing.
    if (rValues.IsSetStrainVector()) {
        Vector& strain_out_kratos = rValues.GetStrainVector();
        if (strain_out_kratos.size() != kVoigt) strain_out_kratos.resize(kVoigt, false);
        for (std::size_t i = 0; i < kVoigt; ++i) strain_out_kratos[i] = strain_now_kratos[i];
    }

    std::array<double, kVoigt> sig_n_abq = sig_n_abq_rotated;
    std::array<double, SandHypoCpp::kStateDim> q_n = q_n_rotated;

    const auto step = SandHypoCpp::integrate_step(sig_n_abq.data(),
                                                  q_n.data(),
                                                  deps_abq.data(),
                                                  check.parms,
                                                  dtime,
                                                  mDtsubPersistent,
                                                  kDefaultErrTol,
                                                  kDefaultMaxNInt,
                                                  kDefaultDtMin);

    using SandHypoCpp::FailureClass;
    KRATOS_ERROR_IF(step.failure_class != FailureClass::Ok)
        << "SandHypoplastic3DLaw::Finalize: kernel returned non-Ok failure_class. "
        << "Calculate* should have caught this earlier.";

    const double tr_deps = deps_abq[0] + deps_abq[1] + deps_abq[2];
    const double pore_next = mPore - bulk_w * tr_deps;
    KRATOS_ERROR_IF_NOT(std::isfinite(pore_next) && std::abs(pore_next) <= 1.0e30)
        << "SandHypoplastic3DLaw::Finalize: pore_next = " << pore_next
        << " is non-finite or > 1e30 (F17 guard). bulk_w=" << bulk_w
        << ", tr_deps=" << tr_deps << ", mPore=" << mPore << ".";

    // Publish stress + tangent to rValues so MPM persists them. MPM uses a
    // FRESH Parameters object in FinalizeSolutionStep
    // (see MPMUpdatedLagrangian): the stress vector that was published in
    // Calculate is GONE. After this Finalize returns, the element copies
    // rVariables.StressVector back to mMP.cauchy_stress_vector.
    // Skipping the publish would leave the material point with
    // zero (or stale) stress while the law's private state advances.
    //
    // Both writes are option-guarded. MPM Finalize
    // sets COMPUTE_STRESS only -- the tangent guard below evaluates to
    // false on the MPM path and the unrequested tangent buffer is left
    // untouched. Non-MPM callers that omit SetStressVector entirely also
    // pass through cleanly.
    const Flags& fin_options = rValues.GetOptions();
    const bool fin_compute_stress  = fin_options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool fin_compute_tangent = fin_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    std::array<double, kVoigt> sig_total_next_abq;
    std::memcpy(sig_total_next_abq.data(), step.sig, sizeof(sig_total_next_abq));
    EffectiveToTotalAbq(sig_total_next_abq, pore_next);

    // Same finite/bounded guard on the post-conversion total stress as in Calculate.
    for (std::size_t i = 0; i < kVoigt; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(sig_total_next_abq[i])
                            && std::abs(sig_total_next_abq[i]) <= 1.0e30)
            << "SandHypoplastic3DLaw::Finalize: sig_total_next[" << i << "] = "
            << sig_total_next_abq[i] << " is non-finite or > 1e30.";
    }

    if (fin_compute_stress) {
        Vector& stress_out_kratos = rValues.GetStressVector();
        if (stress_out_kratos.size() != kVoigt) stress_out_kratos.resize(kVoigt, false);
        SwapVectorAbqToKratos(sig_total_next_abq, stress_out_kratos);
    }

    if (fin_compute_tangent) {
        Matrix& tangent_out_kratos = rValues.GetConstitutiveMatrix();
        if (tangent_out_kratos.size1() != kVoigt || tangent_out_kratos.size2() != kVoigt) {
            tangent_out_kratos.resize(kVoigt, kVoigt, false);
        }
        std::array<std::array<double, kVoigt>, kVoigt> D_abq{};
        for (std::size_t i = 0; i < kVoigt; ++i)
            for (std::size_t j = 0; j < kVoigt; ++j)
                D_abq[i][j] = step.D[i][j];
        SwapTangentAbqToKratos(D_abq, tangent_out_kratos);
        AddBulkWNormalBlock(tangent_out_kratos, bulk_w);
    }

    // Commit persisted state. The kernel returned σ_n+1 and δ_n+1 in the
    // end-of-step (= current/t_n+1) spatial frame, which is what MPM and
    // the next-step caller expect -- no post-rotation needed.
    std::memcpy(mSigEffAbq.data(), step.sig, sizeof(double) * kVoigt);
    for (std::size_t i = 0; i < SandHypoCpp::kStateDim; ++i) mStateAbq[i] = step.q[i];
    mPore = pore_next;
    mDtsubPersistent = step.dtsub_next;
    // Commit F_total so next step's F_rel = F_next * F_total^-1.
    mF0_3x3 = F_total;
}

//********************************************************************
// Check
//********************************************************************

int SandHypoplastic3DLaw::Check(const Properties& rMaterialProperties,
                                 const GeometryType& /*rElementGeometry*/,
                                 const ProcessInfo& /*rCurrentProcessInfo*/) const
{
    // Explicit MPM setups are deliberately accepted at Check time.
    // Explicit MPM's call sequence (see MPMUpdatedLagrangian::
    // CalculateOnIntegrationPoints):
    //
    //     CalculateOnIntegrationPoints(CALCULATE_EXPLICIT_MP_STRESS, ...)
    //       -> CalculateExplicitStresses           (calls this law's Calculate)
    //       -> FinalizeStepVariables               (writes mMP.cauchy_stress_vector
    //                                                from rVariables.StressVector)
    //
    // commits stress to the MP via the element, so `FinalizeMaterialResponseCauchy`
    // never being called in explicit (FinalizeSolutionStep errors out
    // for explicit) is intentional and safe -- the law's
    // own internal state (mSigEffAbq / mStateAbq / mPore / mDtsubPersistent /
    // mF0_3x3) commits at the END of `CalculateMaterialResponseCauchy`'s
    // explicit branch instead.

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PARAMETERS))
        << "SandHypoplastic3DLaw::Check: SAND_HYPOPLASTIC_PARAMETERS missing in Properties.";
    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PARAMETERS];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SandHypoplastic3DLaw::Check: SAND_HYPOPLASTIC_PARAMETERS must have 16 entries.";

    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SandHypoplastic3DLaw::Check: props[" << i << "] is non-finite.";
    }

    SandHypoCpp::Props16 buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SandHypoplastic3DLaw::Check: check_parms rejected '"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    return 0;
}

//********************************************************************
// Voigt swap helpers
//********************************************************************

void SandHypoplastic3DLaw::SwapVectorKratosToAbq(const Vector& kratos6,
                                                  std::array<double, 6>& abq6)
{
    for (std::size_t i = 0; i < 6; ++i) abq6[i] = kratos6[SwapIdx(i)];
}

void SandHypoplastic3DLaw::SwapVectorAbqToKratos(const std::array<double, 6>& abq6,
                                                  Vector& kratos6)
{
    if (kratos6.size() != 6) kratos6.resize(6, false);
    for (std::size_t i = 0; i < 6; ++i) kratos6[i] = abq6[SwapIdx(i)];
}

void SandHypoplastic3DLaw::SwapTangentAbqToKratos(const std::array<std::array<double, 6>, 6>& abq,
                                                   Matrix& kratos6x6)
{
    if (kratos6x6.size1() != 6 || kratos6x6.size2() != 6) kratos6x6.resize(6, 6, false);
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j)
            kratos6x6(i, j) = abq[SwapIdx(i)][SwapIdx(j)];
}

//********************************************************************
// total <-> effective (Abaqus order; normal slots 0..2 only)
//********************************************************************

void SandHypoplastic3DLaw::TotalToEffectiveAbq(std::array<double, 6>& sigAbq, double pore)
{
    sigAbq[0] += pore;
    sigAbq[1] += pore;
    sigAbq[2] += pore;
}

void SandHypoplastic3DLaw::EffectiveToTotalAbq(std::array<double, 6>& sigAbq, double pore)
{
    sigAbq[0] -= pore;
    sigAbq[1] -= pore;
    sigAbq[2] -= pore;
}

//********************************************************************
// Total Almansi strain from deformation gradient
//********************************************************************

void SandHypoplastic3DLaw::CalculateAlmansiFromF(const Matrix& F,
                                                  std::array<double, 6>& e_kratos)
{
    Matrix b(3, 3);
    noalias(b) = prod(F, trans(F));
    Matrix b_inv(3, 3);
    double det_b = 0.0;
    MathUtils<double>::InvertMatrix(b, b_inv, det_b);
    e_kratos[0] = 0.5 * (1.0 - b_inv(0, 0));
    e_kratos[1] = 0.5 * (1.0 - b_inv(1, 1));
    e_kratos[2] = 0.5 * (1.0 - b_inv(2, 2));
    // Kratos Voigt order [11, 22, 33, 12, 23, 13]
    e_kratos[3] = -b_inv(0, 1); // gamma_12 = 2 * e_12
    e_kratos[4] = -b_inv(1, 2); // gamma_23 = 2 * e_23 -- Kratos slot 4
    e_kratos[5] = -b_inv(0, 2); // gamma_13 = 2 * e_13 -- Kratos slot 5
}

//********************************************************************
// Polar decomposition + log(U)
//********************************************************************

void SandHypoplastic3DLaw::PolarDecompFAndLogU(const Matrix& F_rel,
                                                Matrix& R_out,
                                                std::array<double, 6>& logU_voigt_abq)
{
    // C = F^T F  (right Cauchy-Green; symmetric positive definite)
    Matrix C(3, 3);
    noalias(C) = prod(trans(F_rel), F_rel);

    // Spectral decomposition C = B * Λ * B^T (Kratos convention: B columns are
    // eigenvectors, Λ diagonal). GaussSeidelEigenSystem is iterative Jacobi.
    Matrix B(3, 3), Lambda(3, 3);
    const bool converged = MathUtils<double>::GaussSeidelEigenSystem(C, B, Lambda);
    KRATOS_ERROR_IF_NOT(converged)
        << "SandHypoplastic3DLaw::PolarDecompFAndLogU: GaussSeidelEigenSystem failed to converge on C = F^T F.";

    // Build √Λ, 1/√Λ, log(√Λ) on the diagonal. Eigenvalues must be > 0 for a
    // valid polar decomposition (F_rel must be invertible / orientation-
    // preserving). Guard non-positivity defensively.
    Matrix SqrtLam(3, 3, 0.0), InvSqrtLam(3, 3, 0.0), LogSqrtLam(3, 3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        const double l = Lambda(i, i);
        KRATOS_ERROR_IF(l <= 0.0)
            << "SandHypoplastic3DLaw::PolarDecompFAndLogU: non-positive eigenvalue of C ("
            << l << "); F_rel is singular or has det <= 0.";
        const double sl = std::sqrt(l);
        SqrtLam(i, i)    = sl;
        InvSqrtLam(i, i) = 1.0 / sl;
        LogSqrtLam(i, i) = std::log(sl);
    }

    // U^-1   = B * (1/√Λ) * B^T   (just used to recover R)
    // log(U) = B * log(√Λ) * B^T  (Hencky strain of the stretch, in the
    //          BEGIN-of-step / right-stretch basis)
    Matrix Uinv(3, 3);
    noalias(Uinv)  = prod(B, Matrix(prod(InvSqrtLam, trans(B))));
    Matrix LogU(3, 3);
    noalias(LogU)  = prod(B, Matrix(prod(LogSqrtLam, trans(B))));

    // R = F_rel * U^-1
    R_out = prod(F_rel, Uinv);

    // Convert log(U) (begin-of-step basis) -> log(V) = R * log(U) * R^T
    // (end-of-step / current basis). Without this, the law would feed
    // the kernel a strain increment in the OLD axes while the stress and
    // intergranular strain that the same law pre-rotates by R land in
    // the NEW axes -- combined rotation+stretch would silently corrupt
    // the constitutive update.
    Matrix RLogU(3, 3);
    noalias(RLogU) = prod(R_out, LogU);
    Matrix LogV(3, 3);
    noalias(LogV)  = prod(RLogU, trans(R_out));

    // log(V) -> Voigt-Abaqus strain with engineering shear convention.
    // Abaqus order: [11, 22, 33, 12, 13, 23].
    logU_voigt_abq[0] = LogV(0, 0);
    logU_voigt_abq[1] = LogV(1, 1);
    logU_voigt_abq[2] = LogV(2, 2);
    logU_voigt_abq[3] = 2.0 * LogV(0, 1); // gamma_12
    logU_voigt_abq[4] = 2.0 * LogV(0, 2); // gamma_13  (Abaqus slot 4)
    logU_voigt_abq[5] = 2.0 * LogV(1, 2); // gamma_23  (Abaqus slot 5)
}

//********************************************************************
// Voigt-Abaqus tensor rotation: σ' = R · σ · R^T (stress, tensorial)
//********************************************************************

void SandHypoplastic3DLaw::RotateStressVoigtAbq(const std::array<double, 6>& sig_in,
                                                 const Matrix& R,
                                                 std::array<double, 6>& sig_out)
{
    Matrix M(3, 3);
    M(0, 0) = sig_in[0]; M(1, 1) = sig_in[1]; M(2, 2) = sig_in[2];
    M(0, 1) = M(1, 0) = sig_in[3];          // 12
    M(0, 2) = M(2, 0) = sig_in[4];          // 13 (Abaqus slot 4)
    M(1, 2) = M(2, 1) = sig_in[5];          // 23 (Abaqus slot 5)
    Matrix RM(3, 3);
    noalias(RM) = prod(R, M);
    Matrix Rot(3, 3);
    noalias(Rot) = prod(RM, trans(R));
    sig_out[0] = Rot(0, 0); sig_out[1] = Rot(1, 1); sig_out[2] = Rot(2, 2);
    sig_out[3] = Rot(0, 1);
    sig_out[4] = Rot(0, 2);
    sig_out[5] = Rot(1, 2);
}

//********************************************************************
// Voigt-Abaqus tensor rotation: ε' = R · ε · R^T (strain, engineering shear)
//********************************************************************

void SandHypoplastic3DLaw::RotateStrainVoigtAbq(const std::array<double, 6>& eps_in,
                                                 const Matrix& R,
                                                 std::array<double, 6>& eps_out)
{
    Matrix M(3, 3);
    M(0, 0) = eps_in[0]; M(1, 1) = eps_in[1]; M(2, 2) = eps_in[2];
    M(0, 1) = M(1, 0) = 0.5 * eps_in[3];    // γ_12 / 2 = ε_12
    M(0, 2) = M(2, 0) = 0.5 * eps_in[4];    // (Abaqus slot 4)
    M(1, 2) = M(2, 1) = 0.5 * eps_in[5];    // (Abaqus slot 5)
    Matrix RM(3, 3);
    noalias(RM) = prod(R, M);
    Matrix Rot(3, 3);
    noalias(Rot) = prod(RM, trans(R));
    eps_out[0] = Rot(0, 0); eps_out[1] = Rot(1, 1); eps_out[2] = Rot(2, 2);
    eps_out[3] = 2.0 * Rot(0, 1);
    eps_out[4] = 2.0 * Rot(0, 2);
    eps_out[5] = 2.0 * Rot(1, 2);
}

//********************************************************************
// Bauer void-ratio init (Fortran umat .for:171-185)
//********************************************************************

double SandHypoplastic3DLaw::InitVoidRatioBauer(const Vector& props16,
                                                 const std::array<double, 6>& sigTotalAbq)
{
    // ameanstress is COMPRESSION-positive, computed from TOTAL stress
    // (which is tension-positive).
    const double ameanstress = -(sigTotalAbq[0] + sigTotalAbq[1] + sigTotalAbq[2]) / 3.0;
    const double e_0_encoded = props16[15]; // props(16) = parms[15]
    if (e_0_encoded > 10.0) {
        return e_0_encoded - 10.0;          // raw encoding override
    }
    if (ameanstress < 0.001) {
        return e_0_encoded;                  // near-zero stress: raw e_0
    }
    const double hs = props16[2];            // props(3)
    const double en = props16[3];            // props(4)
    return e_0_encoded * std::exp(-std::pow(3.0 * ameanstress / hs, en));
}

//********************************************************************
// bulk_w * 1\otimes1 normal-normal block (Kratos coords)
//********************************************************************

void SandHypoplastic3DLaw::AddBulkWNormalBlock(Matrix& C_kratos, double bulk_w)
{
    if (bulk_w == 0.0) return;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            C_kratos(i, j) += bulk_w;
}

} // namespace Kratos

// Path B - SandHypoplasticCppLaw bridge (implementation).
// See sand_hypoplastic_cpp_law.h and the path-B-cpp-port memory for the
// kernel/bridge contract. Effective-Cauchy in/out at the kernel boundary;
// this file does total<->effective + slot 4<->5 swap + bulk_w tangent
// correction + Bauer void init + failure_class -> KRATOS_ERROR / pnewdt
// routing.
#include "sand_hypoplastic_cpp_law.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <string>

#include "includes/checks.h"
#include "includes/mat_variables.h"
#include "utilities/math_utils.h"

// MPM IS_EXPLICIT lives in MPMApplication's variables header. Pulled in so
// the bridge can fail closed on explicit-mode contracts the round-1 design
// does not yet honor (see GetLawFeatures docstring + Codex review 2026-05-11
// finding F3).
//
// SandHypoplasticApplication already depends on MPMApplication through its
// CMakeLists (`kratos_add_dependency(.../MPMApplication)` +
// `target_link_libraries(... KratosMPMCore)`), so this `#include` does NOT
// introduce a NEW link-time dependency -- only a header-level coupling on
// the IS_EXPLICIT symbol. The round-1 bridge accepts that coupling because
// the MPMApplication-side IS_EXPLICIT is the canonical signal for the
// failure mode it guards (explicit MPM skipping FinalizeSolutionStep). If
// Path B is ever needed in a SandHypoplasticApplication build that drops
// MPMApplication, replace this include with a runtime
// `KratosComponents<Variable<bool>>::Has("IS_EXPLICIT")` lookup. Codex
// 2026-05-11 F6 (medium): documented coupling, no new build dependency.
#include "mpm_application_variables.h"

namespace Kratos
{

namespace
{
constexpr std::size_t kVoigt = 6;
constexpr std::size_t kPropsCount = 16;

// Default kernel substep tolerance / limits. Match the Fortran umat caller
// defaults; surfacing them in a future Properties read is a follow-up.
constexpr double kDefaultErrTol = 1.0e-3;
constexpr int    kDefaultMaxNInt = 10000;
constexpr double kDefaultDtMin   = 1.0e-17;

// Kratos uses [11,22,33,12,23,13]; Path B kernel (Abaqus) uses
// [11,22,33,12,13,23]. Slots 4 and 5 swap; the 0..3 range is identical.
inline std::size_t SwapIdx(std::size_t i)
{
    if (i == 4) return 5;
    if (i == 5) return 4;
    return i;
}

// ====================================================================
// Per-MP debug-dump (post-round-4 diagnostic, 2026-05-12)
// ====================================================================
//
// When environment variables `PATH_B_DEBUG_DUMP_STEP_MIN` and
// `PATH_B_DEBUG_DUMP_STEP_MAX` are both set (inclusive range, integer
// step indices), every Calculate call within that step range writes a
// JSON-lines record to `PATH_B_DEBUG_DUMP_FILE` (default
// `path_b_debug_dump.jsonl` in the current working directory).
//
// Designed for the use case "Kratos MPM crashes natively at step N
// AFTER the bridge returns; what changed between step N-1 (success)
// and step N (death) in the per-MP kernel input/output?" — set
// MIN=N-1 MAX=N, run the failing sim, then diff the resulting JSONL.
//
// Zero overhead when the env vars are unset (single static-bool
// check). Thread-safe via a global mutex (cheap because the dump
// path is only hot during debug sessions). Best run with
// OMP_NUM_THREADS=1 to keep per-MP call ordering deterministic across
// repeated runs.
//
// The record schema is documented in
// memory/reference_kratos_mpm_3d_singularity.md.

struct DebugDumper
{
    bool enabled = false;
    int step_min = -1;
    int step_max = -1;
    std::ofstream stream;
    std::mutex mutex;
    int current_step = -1;
    int call_count_in_step = 0;

    static DebugDumper& Instance()
    {
        static DebugDumper inst;
        return inst;
    }

    void EnsureInit()
    {
        std::lock_guard<std::mutex> g(mutex);
        if (enabled || step_min >= 0) return; // already initialized
        const char* v_min = std::getenv("PATH_B_DEBUG_DUMP_STEP_MIN");
        const char* v_max = std::getenv("PATH_B_DEBUG_DUMP_STEP_MAX");
        std::cerr << "[PATH_B_DEBUG] EnsureInit: PATH_B_DEBUG_DUMP_STEP_MIN="
                  << (v_min ? v_min : "(unset)")
                  << " PATH_B_DEBUG_DUMP_STEP_MAX="
                  << (v_max ? v_max : "(unset)") << std::endl;
        if (!v_min || !v_max || !*v_min || !*v_max) {
            step_min = -2; // sentinel "checked, not requested"
            return;
        }
        step_min = std::atoi(v_min);
        step_max = std::atoi(v_max);
        if (step_min < 0 || step_max < step_min) {
            std::cerr << "[PATH_B_DEBUG] bad range: " << step_min << "-" << step_max << std::endl;
            step_min = -2;
            return;
        }
        const char* v_file = std::getenv("PATH_B_DEBUG_DUMP_FILE");
        std::string path = (v_file && *v_file) ? v_file : "path_b_debug_dump.jsonl";
        stream.open(path, std::ios::out | std::ios::trunc);
        if (!stream.is_open()) {
            std::cerr << "[PATH_B_DEBUG] FAILED to open dump file: " << path << std::endl;
            step_min = -2;
            return;
        }
        stream << std::setprecision(15);
        enabled = true;
        std::cerr << "[PATH_B_DEBUG] init OK: dumping steps " << step_min
                  << ".." << step_max << " -> " << path << std::endl;
    }

    bool ShouldDumpForStep(int step)
    {
        EnsureInit();
        return enabled && step >= step_min && step <= step_max;
    }

    // Bump per-step call counter. Caller must hold the dumper mutex.
    int NextCallIdLocked(int step)
    {
        if (step != current_step) {
            current_step = step;
            call_count_in_step = 0;
        }
        return call_count_in_step++;
    }
};

// Helpers to format arrays / matrices as JSON arrays.
inline void DumpArrayJson(std::ostream& s, const double* p, std::size_t n)
{
    s << "[";
    for (std::size_t i = 0; i < n; ++i) {
        if (i) s << ",";
        s << p[i];
    }
    s << "]";
}
inline void DumpVoigtJson(std::ostream& s, const std::array<double, 6>& v)
{
    DumpArrayJson(s, v.data(), 6);
}
inline void DumpMatrix3x3Json(std::ostream& s, const Matrix& M)
{
    s << "[";
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            if (i || j) s << ",";
            s << M(i, j);
        }
    }
    s << "]";
}
inline void DumpMatrix6x6Json(std::ostream& s, const Matrix& M)
{
    s << "[";
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            if (i || j) s << ",";
            s << M(i, j);
        }
    }
    s << "]";
}

} // anonymous namespace

//********************************************************************
// Lifecycle / Clone
//********************************************************************

SandHypoplasticCppLaw::SandHypoplasticCppLaw()
    : ConstitutiveLaw(), mF0_3x3(IdentityMatrix(3, 3)) {}

SandHypoplasticCppLaw::SandHypoplasticCppLaw(const SandHypoplasticCppLaw& rOther)
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

ConstitutiveLaw::Pointer SandHypoplasticCppLaw::Clone() const
{
    return Kratos::make_shared<SandHypoplasticCppLaw>(*this);
}

//********************************************************************
// Features / metadata
//********************************************************************

void SandHypoplasticCppLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Deformation_Gradient is the canonical strain measure consumed by
    // the bridge in BOTH implicit and explicit MPM paths (via polar
    // decomposition + log(V) from F). MPMUpdatedLagrangian::Check
    // (mpm_updated_lagrangian.cpp:1728-1743) accepts Deformation_Gradient
    // unconditionally for implicit and as one of the two accepted measures
    // for explicit. Round 2 (2026-05-11) added explicit support by removing
    // the bridge's own IS_EXPLICIT Check-time rejection and adding a
    // commit-on-explicit branch at the end of CalculateMaterialResponseCauchy.
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize = static_cast<SizeType>(kVoigt);
    rFeatures.mSpaceDimension = static_cast<SizeType>(WorkingSpaceDimension());
}

//********************************************************************
// Has / GetValue / SetValue
//********************************************************************

bool SandHypoplasticCppLaw::Has(const Variable<double>& rThisVariable)
{
    return false;
}

bool SandHypoplasticCppLaw::Has(const Variable<Vector>& rThisVariable)
{
    // Expose INTERNAL_VARIABLES for cross-path CSV comparison harness
    // (Path B smoke-test against Path A baselines). 15-slot layout:
    //   [0..5] intergranular strain delta, Kratos Voigt order (slot 4<->5
    //          swap applied from kernel-internal Abaqus order)
    //   [6]   void ratio
    //   [7]   -pore (Fortran statev(8) sign-flipped convention)
    //   [8]   dtsub_persistent
    //   [9..14] reserved / zero
    return rThisVariable == INTERNAL_VARIABLES;
}

double& SandHypoplasticCppLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    rValue = 0.0;
    return rValue;
}

Vector& SandHypoplasticCppLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == INTERNAL_VARIABLES) {
        constexpr std::size_t kStatevSize = 15;
        if (rValue.size() != kStatevSize) rValue.resize(kStatevSize, false);
        // 0..5: del (IS) in Kratos Voigt order. SwapIdx(i) maps the
        // Abaqus-internal indices to Kratos-external indices (slots 4 and
        // 5 swap; 0..3 unchanged). Convention check at end-of-bridge: a
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

void SandHypoplasticCppLaw::SetValue(const Variable<double>& /*rThisVariable*/,
                                     const double& /*rValue*/,
                                     const ProcessInfo& /*rCurrentProcessInfo*/)
{
}

void SandHypoplasticCppLaw::SetValue(const Variable<Vector>& /*rThisVariable*/,
                                     const Vector& /*rValue*/,
                                     const ProcessInfo& /*rCurrentProcessInfo*/)
{
}

//********************************************************************
// InitializeMaterial
//********************************************************************

void SandHypoplasticCppLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                               const GeometryType& /*rElementGeometry*/,
                                               const Vector& /*rShapeFunctionsValues*/)
{
    // Validate the props vector once at material init so the user gets a
    // fast KRATOS_ERROR rather than a silent first-step contract violation.
    // Always run -- props validation is cheap and helps surface bad JSON
    // during model setup even on a deserialized law that already has
    // committed history (mInitialized == true).
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PROPS_16))
        << "SandHypoplasticCppLaw requires SAND_HYPOPLASTIC_PROPS_16 in Properties.";

    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PROPS_16];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SAND_HYPOPLASTIC_PROPS_16 must have 16 entries; got " << props_vec.size();

    // Fortran's check_parms is IEEE-blind; we enforce finiteness here so
    // bad JSON inputs surface during model setup instead of mid-solve.
    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SAND_HYPOPLASTIC_PROPS_16[" << i << "] is non-finite (" << props_vec[i] << ").";
    }

    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SAND_HYPOPLASTIC_PROPS_16 rejected by check_parms: failed='"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    // Idempotent reset (Codex 2026-05-11 F8 high finding). InitializeMaterial
    // is called once per material point at simulation start, AND may be
    // called again after deserialization (restart / partition migration)
    // depending on caller lifecycle. The committed history members
    // (mSigEffAbq, mStateAbq, mPore, mDtsubPersistent, mFirstCallStateLoaded,
    // mF0_3x3) are now serialized via save/load; clobbering them
    // here would silently discard restored history and the F5 fix would be
    // dead weight. Reset only when this is the FIRST init -- detected via
    // mInitialized still false. Kratos-registered prototypes start with
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

void SandHypoplasticCppLaw::ResetMaterial(const Properties& rMaterialProperties,
                                          const GeometryType& /*rElementGeometry*/,
                                          const Vector& /*rShapeFunctionsValues*/)
{
    // Re-validate props on every reset -- catches the case where Reset is
    // called with a modified material (e.g. user adjusted parameters and
    // wants to reset history) and the new values are bad.
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PROPS_16))
        << "SandHypoplasticCppLaw::ResetMaterial requires SAND_HYPOPLASTIC_PROPS_16 in Properties.";
    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PROPS_16];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SAND_HYPOPLASTIC_PROPS_16 must have 16 entries; got " << props_vec.size();
    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SAND_HYPOPLASTIC_PROPS_16[" << i << "] is non-finite at ResetMaterial.";
    }
    SandHypoCpp::Props16 buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SAND_HYPOPLASTIC_PROPS_16 rejected by check_parms in ResetMaterial: '"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    // Clear all committed bridge history. mInitialized stays true -- after
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

void SandHypoplasticCppLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    KRATOS_ERROR_IF_NOT(mInitialized)
        << "SandHypoplasticCppLaw::CalculateMaterialResponseCauchy called before InitializeMaterial.";

    const Properties& props = rValues.GetMaterialProperties();
    const Vector& props_vec = props[SAND_HYPOPLASTIC_PROPS_16];
    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    if (check.error != 0) {
        KRATOS_ERROR << "SandHypoplasticCppLaw: check_parms failed at solve time: '"
                     << (check.failed ? check.failed : "?") << "' bad_value="
                     << check.bad_value;
    }
    const double bulk_w = check.parms[14];

    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double dtime = process_info[DELTA_TIME];

    // Round 2 (explicit MPM support): detect IS_EXPLICIT for branching.
    // - In implicit, the bridge defers state commit to `FinalizeMaterialResponseCauchy`
    //   (called once at end of converged step).
    // - In explicit, the element calls Calculate ONCE per step via
    //   `CALCULATE_EXPLICIT_MP_STRESS` hook (.cpp:1496-1502) and never calls
    //   Finalize, so the bridge must commit internal state at end of Calculate
    //   when is_explicit. The element's `FinalizeStepVariables` (line 1501)
    //   handles the `mMP.cauchy_stress_vector` writeback separately.
    const bool is_explicit = (process_info.Has(IS_EXPLICIT))
                                 ? process_info.GetValue(IS_EXPLICIT)
                                 : false;

    // Option flags drive what the caller wants this Calculate to publish.
    // The bridge must NOT touch the stress / tangent output buffers when
    // the corresponding compute flag is unset -- some callers omit
    // SetStressVector / SetConstitutiveMatrix entirely, and GetStressVector
    // / GetConstitutiveMatrix on an unset pointer KRATOS_DEBUG_ERRORs
    // (Codex 2026-05-11 F9 medium finding).
    const Flags& options = rValues.GetOptions();
    const bool compute_stress  = options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool compute_tangent = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    // F11 fix: read F here for both paths -- the implicit branch needs it
    // for polar decomposition + Almansi-from-F write-back; the explicit
    // branch only uses it as a defense-in-depth kinematic-gate / restart
    // hand-off (round-3 explicit kinematics come from rValues.StrainVector,
    // not F -- see PROJ-TD-2 closure in TECH_DEBT.md).
    const Matrix& F_total = rValues.GetDeformationGradientF();
    KRATOS_ERROR_IF(F_total.size1() != 3 || F_total.size2() != 3)
        << "DeformationGradientF must be 3x3; got " << F_total.size1() << "x" << F_total.size2();

    // F14 (Codex 2026-05-11): physical-kinematics gate. det(F) <= 0 means
    // self-inverted / collapsed material -- not a valid input to a continuum
    // constitutive update. Catch BEFORE the 3x3 inversion in
    // CalculateAlmansiFromF / before forming F_rel for polar decomposition,
    // because both inversion paths silently produce NaN/Inf on a singular F
    // and the polar decomp can return an improper rotation (det R = -1) for
    // orientation-reversing F (C = F^T F is SPD regardless of sign(det F),
    // so the eigenvalue-positivity check inside PolarDecompFAndLogU does
    // NOT block det(F) < 0). On the explicit path F is typically I or
    // I+symmetric_strain_increment (mpm_explicit_utilities.cpp:346-347);
    // keeping the gate is cheap and rejects upstream contract breakage.
    const double det_F_total = MathUtils<double>::Det(F_total);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_total) && det_F_total > 0.0)
        << "SandHypoplasticCppLaw: det(F_total) = " << det_F_total
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

        mPore = 0.0; // pore is not yet exposed via a Kratos variable; round-2 hook.

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
    // IMPLICIT path (Codex 2026-05-11 F13): F_rel = F_total · F0^-1, polar
    // decomp F_rel = R · U. R pre-rotates stored σ_n and δ_n into the
    // end-of-step spatial frame (Hughes-Winget), log(V) is the
    // corotational strain increment. Kratos MPM does NO host-side stress
    // rotation on the implicit path (mpm_updated_lagrangian.cpp:568,926);
    // objectivity is owned by the law. Almansi-from-F is written back to
    // rValues.StrainVector for MPM's `mMP.almansi_strain_vector`
    // post-processing persistence (.cpp:927).
    //
    // EXPLICIT path (round 3, 2026-05-12, PROJ-TD-2 fix): the element's
    // CalculateExplicitKinematics populates rValues.StrainVector with the
    // Jaumann-corrected cumulative Almansi
    // (mpm_explicit_utilities.cpp:314-322 -> .cpp:565-569) -- THAT is the
    // only strain signal carrying rotation in explicit MPM, because the
    // F handed in via SetDeformationGradientF (.cpp:255) is either I or
    // I+symmetric_strain_increment (.cpp:346-347) and carries no spin.
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
        // ----- IMPLICIT branch (Round 1+2 logic, unchanged) ------------
        std::array<double, kVoigt> strain_now_kratos{};
        CalculateAlmansiFromF(F_total, strain_now_kratos);

        Matrix F0_inv(3, 3);
        double det_F0 = 0.0;
        MathUtils<double>::InvertMatrix(mF0_3x3, F0_inv, det_F0);
        KRATOS_ERROR_IF_NOT(std::isfinite(det_F0) && det_F0 > 0.0)
            << "SandHypoplasticCppLaw: persisted mF0_3x3 has det = " << det_F0
            << " (non-finite or non-positive). Internal invariant violation -- previous "
            << "step committed a corrupt F_total.";
        Matrix F_rel(3, 3);
        noalias(F_rel) = prod(F_total, F0_inv);

        // F14 guard #2: det(F_rel) > 0 (orientation-preserving relative step).
        // Note that C = F_rel^T F_rel is SPD for any non-singular F_rel, so
        // PolarDecompFAndLogU alone cannot reject this case -- it would return
        // an improper rotation (det R = -1).
        const double det_F_rel = MathUtils<double>::Det(F_rel);
        KRATOS_ERROR_IF_NOT(std::isfinite(det_F_rel) && det_F_rel > 0.0)
            << "SandHypoplasticCppLaw: det(F_rel) = " << det_F_rel
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
        // ----- EXPLICIT branch (Round 3, 2026-05-12, PROJ-TD-2 fix) ----
        //
        // Consume the element-populated Jaumann-corrected cumulative
        // Almansi from rValues.GetStrainVector(). Differentiate against
        // mEpsPrevExplicitKratos to get the kernel-driving Δε. No
        // polar decomp, no R-rotation of σ_n/δ_n.
        KRATOS_ERROR_IF_NOT(rValues.IsSetStrainVector())
            << "SandHypoplasticCppLaw: explicit branch requires the element to "
               "set the StrainVector before calling Calculate. MPM normally "
               "populates rValues.GetStrainVector() from "
               "mMP.almansi_strain_vector via SetStrainVector. If you are "
               "calling Calculate manually, attach a Vector(6) via "
               "Parameters::SetStrainVector first.";
        const Vector& strain_now_vec = rValues.GetStrainVector();
        KRATOS_ERROR_IF(strain_now_vec.size() != kVoigt)
            << "SandHypoplasticCppLaw: explicit branch expects StrainVector of "
               "size 6 (Kratos Voigt); got " << strain_now_vec.size() << ".";

        std::array<double, kVoigt> strain_now_kratos{};
        for (std::size_t i = 0; i < kVoigt; ++i) {
            const double v = strain_now_vec[i];
            KRATOS_ERROR_IF_NOT(std::isfinite(v) && std::abs(v) <= 1.0e30)
                << "SandHypoplasticCppLaw: explicit StrainVector[" << i << "] = "
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
                << "SandHypoplasticCppLaw: explicit deps_kratos[" << i << "] = "
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
        // Jaumann spin correction to (mpm_explicit_utilities.cpp:319-322);
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
        KRATOS_ERROR << "SandHypoplasticCppLaw: kernel ContractViolation (error=" << step.error
                     << "); inputs / persisted state corrupt. Not retryable. dtime=" << dtime
                     << ", dtsub_in=" << dtsub_in << ".";
    }
    if (step.failure_class == FailureClass::Fatal) {
        KRATOS_ERROR << "SandHypoplasticCppLaw: kernel Fatal (error=" << step.error
                     << ", parms-only failure inside get_tan / get_F_sig_q).";
    }
    if (step.failure_class == FailureClass::RkfReject) {
        // Caller-driven cutback. Absent a wired Kratos pnewdt channel,
        // the bridge raises a KRATOS_ERROR that includes the literal
        // "pnewdt=0.25" magic substring used by the reference-case driver
        // (`tests/path_a_results/run_reference_cases.py`) to detect a
        // retryable RKF reject and halve d_eps. Round-2: implement a
        // real `CONSTITUTIVE_LAW_PNEWDT`-style cutback. The magic-substring
        // bridge is mirroring Path A's `pnewdt=0.25` rejection convention
        // (Codex 2026-05-11 F16 / round 8 deferred).
        KRATOS_ERROR << "SandHypoplasticCppLaw: kernel RkfReject (error=" << step.error
                     << "); pnewdt=0.25 outer cutback required. dtime=" << dtime
                     << " too large for current state.";
    }

    // ----- Pore evolution: pore_next = pore - bulk_w * tr(deps) -------
    //
    // Mirrors Fortran solout_h .for:1856 in-place rebind. `bulk_w * tr(deps)`
    // uses the same deps that drove the kernel (pre-swap doesn't matter for
    // the trace -- normal slots [0..2] are the same in either order).
    //
    // F17 (Codex 2026-05-11): the pore update happens AFTER integrate_step
    // returns, so kernel guards don't cover it. A typo-scale but finite
    // bulk_w (e.g. 1e29 instead of 2.2e9) times a small volumetric strain
    // can produce a non-finite or overflow pore_next, which then corrupts
    // the total-stress conversion and persists across steps. Verify
    // pore_next finite and bounded before propagating.
    const double tr_deps = deps_abq[0] + deps_abq[1] + deps_abq[2];
    const double pore_next = mPore - bulk_w * tr_deps;
    KRATOS_ERROR_IF_NOT(std::isfinite(pore_next) && std::abs(pore_next) <= 1.0e30)
        << "SandHypoplasticCppLaw: pore_next = " << pore_next
        << " is non-finite or > 1e30 (overflow). bulk_w=" << bulk_w
        << ", tr_deps=" << tr_deps << ", mPore=" << mPore << ".";

    // ----- Convert effective -> total on the bridge output -------------
    //
    // Note: pore_next, NOT mPore -- bulk_w*tr(deps) correction propagates
    // through here (memo Codex 2026-05-02 P2 follow-up). solout_h .for:1856
    // followed by .for:1862.
    std::array<double, kVoigt> sig_total_next_abq;
    std::memcpy(sig_total_next_abq.data(), step.sig, sizeof(sig_total_next_abq));
    EffectiveToTotalAbq(sig_total_next_abq, pore_next);

    // F17 #2: verify the post-conversion total stress is finite and bounded
    // before publishing.
    for (std::size_t i = 0; i < kVoigt; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(sig_total_next_abq[i])
                            && std::abs(sig_total_next_abq[i]) <= 1.0e30)
            << "SandHypoplasticCppLaw: sig_total_next[" << i << "] = "
            << sig_total_next_abq[i] << " is non-finite or > 1e30.";
    }

    // ----- Abaqus -> Kratos swap, option-guarded (F9) ------------------
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

    // ----- Per-MP debug dump (post-round-4 diagnostic) ----------------
    // No-op when env vars PATH_B_DEBUG_DUMP_STEP_MIN / _MAX are unset.
    // When set, writes one JSON line per Calculate call within the
    // configured step range to PATH_B_DEBUG_DUMP_FILE. Captures the
    // bridge-side input and output as they would reach Kratos, so the
    // user can diff successful-step records against the step where
    // Kratos crashes after the bridge returns.
    {
        const int step_index = (process_info.Has(STEP)) ? process_info.GetValue(STEP) : -1;
        auto& dumper = DebugDumper::Instance();
        if (dumper.ShouldDumpForStep(step_index)) {
            std::lock_guard<std::mutex> g(dumper.mutex);
            const int call_id = dumper.NextCallIdLocked(step_index);
            std::ostringstream out;
            out << std::setprecision(15);
            out << "{"
                << "\"step\":" << step_index
                << ",\"call\":" << call_id
                << ",\"dtime\":" << dtime
                << ",\"is_explicit\":" << (is_explicit ? "true" : "false")
                << ",\"first_call_state_loaded\":" << (mFirstCallStateLoaded ? "true" : "false")
                << ",\"det_F\":" << det_F_total
                << ",\"F\":";
            DumpMatrix3x3Json(out, F_total);
            out << ",\"deps_abq\":";
            DumpVoigtJson(out, deps_abq);
            out << ",\"sig_n_abq_in\":";
            DumpVoigtJson(out, sig_n_abq);
            out << ",\"q_n_in\":";
            DumpArrayJson(out, q_n.data(), SandHypoCpp::kStateDim);
            out << ",\"pore_in\":" << mPore
                << ",\"dtsub_in\":" << dtsub_in
                << ",\"failure_class\":" << static_cast<int>(step.failure_class)
                << ",\"error\":" << step.error
                << ",\"nfev\":" << step.nfev
                << ",\"used_elastic\":" << (step.used_elastic ? "true" : "false")
                << ",\"dtsub_next\":" << step.dtsub_next
                << ",\"pore_next\":" << pore_next
                << ",\"sig_eff_out\":";
            DumpArrayJson(out, step.sig, 6);
            out << ",\"q_out\":";
            DumpArrayJson(out, step.q, SandHypoCpp::kStateDim);
            out << ",\"sig_total_next_abq\":";
            DumpVoigtJson(out, sig_total_next_abq);
            // D in Abaqus order (post-swap to Kratos happens only in
            // the Kratos output buffer above; the raw kernel D is more
            // useful for diff analysis since it's the law's intrinsic
            // output before any Voigt-swap or bulk_w correction).
            out << ",\"D_abq\":[";
            for (std::size_t i = 0; i < 6; ++i)
                for (std::size_t j = 0; j < 6; ++j) {
                    if (i || j) out << ",";
                    out << step.D[i][j];
                }
            out << "]}";
            dumper.stream << out.str() << "\n";
            dumper.stream.flush();
        }
    }

    // ----- State commit policy ----------------------------------------
    //
    // IMPLICIT path: persisted state is NOT mutated here. The element may
    // call Calculate multiple times during NR iterations with different
    // trial F; we publish stress / tangent but defer the commit to
    // FinalizeMaterialResponseCauchy (called once at end of converged step).
    //
    // EXPLICIT path: the element calls Calculate ONCE per step via the
    // `CALCULATE_EXPLICIT_MP_STRESS` hook and NEVER calls Finalize
    // (mpm_updated_lagrangian.cpp:882-883 errors out FinalizeSolutionStep
    // for explicit). So the bridge MUST commit internal state here, at the
    // end of Calculate, otherwise mSigEffAbq / mStateAbq / mPore /
    // mDtsubPersistent / mEpsPrevExplicitKratos / mF0_3x3 stay at initial
    // values forever and the bridge is effectively stateless across
    // explicit steps. The element's `FinalizeStepVariables` (line 1501,
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

        // Round 3: persist the begin-of-(next-)step baseline for the
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

void SandHypoplasticCppLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // Re-run the integration with the same inputs to recover scratch state,
    // then commit to the persisted members. This deliberately mirrors the
    // Calculate path -- if anything diverges, the kernel's own determinism
    // (no RNG / no internal time-dependent caching) guarantees the second
    // run produces the same StepResult.
    KRATOS_ERROR_IF_NOT(mInitialized)
        << "SandHypoplasticCppLaw::FinalizeMaterialResponseCauchy called before InitializeMaterial.";

    const Properties& props = rValues.GetMaterialProperties();
    const Vector& props_vec = props[SAND_HYPOPLASTIC_PROPS_16];
    SandHypoCpp::Props16 props_buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) props_buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(props_buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SandHypoplasticCppLaw::Finalize: check_parms failed: '"
        << (check.failed ? check.failed : "?") << "'.";
    const double bulk_w = check.parms[14];

    const ProcessInfo& process_info = rValues.GetProcessInfo();
    const double dtime = process_info[DELTA_TIME];

    // F11 fix: same Almansi-from-F path as Calculate -- the element's
    // CalculateKinematics is called again on the Finalize side and again
    // doesn't populate rVariables.StrainVector. Compute Almansi from F
    // here too.
    const Matrix& F_total = rValues.GetDeformationGradientF();
    KRATOS_ERROR_IF(F_total.size1() != 3 || F_total.size2() != 3)
        << "DeformationGradientF must be 3x3 in Finalize; got " << F_total.size1() << "x" << F_total.size2();

    // F14 (Codex 2026-05-11): same kinematic gate as Calculate.
    const double det_F_total = MathUtils<double>::Det(F_total);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_total) && det_F_total > 0.0)
        << "SandHypoplasticCppLaw::Finalize: det(F_total) = " << det_F_total
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

    // ----- Same corotational kinematics as Calculate (F13). The kernel is
    // deterministic and Calculate did not mutate persisted state, so this
    // re-derives R and log(U) with identical inputs and identical results.
    Matrix F0_inv(3, 3);
    double det_F0 = 0.0;
    MathUtils<double>::InvertMatrix(mF0_3x3, F0_inv, det_F0);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F0) && det_F0 > 0.0)
        << "SandHypoplasticCppLaw::Finalize: persisted mF0_3x3 has det = " << det_F0
        << " (invariant violation).";
    Matrix F_rel(3, 3);
    noalias(F_rel) = prod(F_total, F0_inv);
    const double det_F_rel = MathUtils<double>::Det(F_rel);
    KRATOS_ERROR_IF_NOT(std::isfinite(det_F_rel) && det_F_rel > 0.0)
        << "SandHypoplasticCppLaw::Finalize: det(F_rel) = " << det_F_rel
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
    // FinalizeStepVariables (.cpp:927) gets a meaningful value into
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
        << "SandHypoplasticCppLaw::Finalize: kernel returned non-Ok failure_class. "
        << "Calculate* should have caught this earlier.";

    const double tr_deps = deps_abq[0] + deps_abq[1] + deps_abq[2];
    const double pore_next = mPore - bulk_w * tr_deps;
    KRATOS_ERROR_IF_NOT(std::isfinite(pore_next) && std::abs(pore_next) <= 1.0e30)
        << "SandHypoplasticCppLaw::Finalize: pore_next = " << pore_next
        << " is non-finite or > 1e30 (F17 guard). bulk_w=" << bulk_w
        << ", tr_deps=" << tr_deps << ", mPore=" << mPore << ".";

    // Publish stress + tangent to rValues so MPM persists them. MPM uses a
    // FRESH Parameters object in FinalizeSolutionStep (mpm_updated_lagrangian
    // .cpp:890-909): the stress vector that was published in Calculate is
    // GONE. After this Finalize returns, the element copies
    // rVariables.StressVector back to mMP.cauchy_stress_vector
    // (.cpp:927). Skipping the publish would leave the material point with
    // zero (or stale) stress while the bridge's private state advances --
    // Codex 2026-05-11 finding F2.
    //
    // Both writes are option-guarded (Codex 2026-05-11 F9). MPM Finalize
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

    // F17 guard on the post-conversion total stress (same as Calculate).
    for (std::size_t i = 0; i < kVoigt; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(sig_total_next_abq[i])
                            && std::abs(sig_total_next_abq[i]) <= 1.0e30)
            << "SandHypoplasticCppLaw::Finalize: sig_total_next[" << i << "] = "
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

int SandHypoplasticCppLaw::Check(const Properties& rMaterialProperties,
                                 const GeometryType& /*rElementGeometry*/,
                                 const ProcessInfo& /*rCurrentProcessInfo*/) const
{
    // Round 2 (explicit MPM support, 2026-05-11): the IS_EXPLICIT Check-time
    // rejection from round 1 (Codex F10) is REMOVED. Explicit MPM's call
    // sequence (`mpm_updated_lagrangian.cpp:1496-1502`):
    //
    //     CalculateOnIntegrationPoints(CALCULATE_EXPLICIT_MP_STRESS, ...)
    //       -> CalculateExplicitStresses           (calls our Calculate)
    //       -> FinalizeStepVariables               (writes mMP.cauchy_stress_vector
    //                                                from rVariables.StressVector)
    //
    // commits stress to the MP via the element, so our `FinalizeMaterialResponseCauchy`
    // never being called in explicit (`.cpp:882-883` errors out
    // FinalizeSolutionStep) is intentional and safe -- the bridge's
    // own internal state (mSigEffAbq / mStateAbq / mPore / mDtsubPersistent /
    // mF0_3x3) commits at the END of `CalculateMaterialResponseCauchy`'s
    // explicit branch instead. Implicit behavior is unchanged.

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SAND_HYPOPLASTIC_PROPS_16))
        << "SandHypoplasticCppLaw::Check: SAND_HYPOPLASTIC_PROPS_16 missing in Properties.";
    const Vector& props_vec = rMaterialProperties[SAND_HYPOPLASTIC_PROPS_16];
    KRATOS_ERROR_IF(props_vec.size() != kPropsCount)
        << "SandHypoplasticCppLaw::Check: SAND_HYPOPLASTIC_PROPS_16 must have 16 entries.";

    for (std::size_t i = 0; i < kPropsCount; ++i) {
        KRATOS_ERROR_IF_NOT(std::isfinite(props_vec[i]))
            << "SandHypoplasticCppLaw::Check: props[" << i << "] is non-finite.";
    }

    SandHypoCpp::Props16 buf;
    for (std::size_t i = 0; i < kPropsCount; ++i) buf[i] = props_vec[i];
    const auto check = SandHypoCpp::check_parms(buf);
    KRATOS_ERROR_IF(check.error != 0)
        << "SandHypoplasticCppLaw::Check: check_parms rejected '"
        << (check.failed ? check.failed : "?") << "' bad_value=" << check.bad_value;

    return 0;
}

//********************************************************************
// Voigt swap helpers
//********************************************************************

void SandHypoplasticCppLaw::SwapVectorKratosToAbq(const Vector& kratos6,
                                                  std::array<double, 6>& abq6)
{
    for (std::size_t i = 0; i < 6; ++i) abq6[i] = kratos6[SwapIdx(i)];
}

void SandHypoplasticCppLaw::SwapVectorAbqToKratos(const std::array<double, 6>& abq6,
                                                  Vector& kratos6)
{
    if (kratos6.size() != 6) kratos6.resize(6, false);
    for (std::size_t i = 0; i < 6; ++i) kratos6[i] = abq6[SwapIdx(i)];
}

void SandHypoplasticCppLaw::SwapTangentAbqToKratos(const std::array<std::array<double, 6>, 6>& abq,
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

void SandHypoplasticCppLaw::TotalToEffectiveAbq(std::array<double, 6>& sigAbq, double pore)
{
    sigAbq[0] += pore;
    sigAbq[1] += pore;
    sigAbq[2] += pore;
}

void SandHypoplasticCppLaw::EffectiveToTotalAbq(std::array<double, 6>& sigAbq, double pore)
{
    sigAbq[0] -= pore;
    sigAbq[1] -= pore;
    sigAbq[2] -= pore;
}

//********************************************************************
// Total Almansi strain from deformation gradient
//********************************************************************

void SandHypoplasticCppLaw::CalculateAlmansiFromF(const Matrix& F,
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

void SandHypoplasticCppLaw::PolarDecompFAndLogU(const Matrix& F_rel,
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
        << "SandHypoplasticCppLaw::PolarDecompFAndLogU: GaussSeidelEigenSystem failed to converge on C = F^T F.";

    // Build √Λ, 1/√Λ, log(√Λ) on the diagonal. Eigenvalues must be > 0 for a
    // valid polar decomposition (F_rel must be invertible / orientation-
    // preserving). Guard non-positivity defensively.
    Matrix SqrtLam(3, 3, 0.0), InvSqrtLam(3, 3, 0.0), LogSqrtLam(3, 3, 0.0);
    for (std::size_t i = 0; i < 3; ++i) {
        const double l = Lambda(i, i);
        KRATOS_ERROR_IF(l <= 0.0)
            << "SandHypoplasticCppLaw::PolarDecompFAndLogU: non-positive eigenvalue of C ("
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
    // (end-of-step / current basis). Without this, the bridge would feed
    // the kernel a strain increment in the OLD axes while the stress and
    // intergranular strain that the same bridge pre-rotates by R land in
    // the NEW axes -- combined rotation+stretch silently corrupts the
    // constitutive update (Codex 2026-05-11 F15 high finding). Codex
    // verified this against Path A's Hughes-Winget comments.
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

void SandHypoplasticCppLaw::RotateStressVoigtAbq(const std::array<double, 6>& sig_in,
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

void SandHypoplasticCppLaw::RotateStrainVoigtAbq(const std::array<double, 6>& eps_in,
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

double SandHypoplasticCppLaw::InitVoidRatioBauer(const Vector& props16,
                                                 const std::array<double, 6>& sigTotalAbq)
{
    // ameanstress is COMPRESSION-positive, computed from TOTAL stress
    // (tension-positive); see memory bridge contract pore equations.
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

void SandHypoplasticCppLaw::AddBulkWNormalBlock(Matrix& C_kratos, double bulk_w)
{
    if (bulk_w == 0.0) return;
    for (std::size_t i = 0; i < 3; ++i)
        for (std::size_t j = 0; j < 3; ++j)
            C_kratos(i, j) += bulk_w;
}

} // namespace Kratos

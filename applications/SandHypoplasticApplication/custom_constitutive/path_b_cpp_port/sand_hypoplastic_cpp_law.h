// Path B - SandHypoplasticCppLaw bridge
// Pure C++ Kratos ConstitutiveLaw wrapper around the Path B kernel
// (sand_hypoplastic_kernel.h). Effective-Cauchy in / out at the kernel
// boundary; this wrapper handles Voigt slot 4<->5 swap, total<->effective
// pore-pressure conversion, void-ratio first-call init (Bauer 1996),
// bulk_w*1\otimes1 tangent correction, and failure_class -> KRATOS_ERROR /
// pnewdt routing.
//
// Path B isolation rule: this file does not include Path A or Path C
// headers. The kernel is the only sand-hypoplastic translation unit it
// depends on; the rest is Kratos core / SandHypoplasticApplication-shared.
#pragma once

#include <array>

#include "includes/constitutive_law.h"
#include "sand_hypoplastic_application_variables.h"

#include "sand_hypoplastic_kernel.h"

namespace Kratos
{

class KRATOS_API(SAND_HYPOPLASTIC_APPLICATION) SandHypoplasticCppLaw
    : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SandHypoplasticCppLaw);

    SandHypoplasticCppLaw();
    SandHypoplasticCppLaw(const SandHypoplasticCppLaw& rOther);
    ~SandHypoplasticCppLaw() override = default;

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override { return 3; }
    SizeType GetStrainSize() const override { return 6; }

    // `Deformation_Gradient` matches the strain-measure gate in
    // MPMUpdatedLagrangian::Check (mpm_updated_lagrangian.cpp:1728-1743):
    // implicit runs require Deformation_Gradient; explicit runs accept
    // EITHER Deformation_Gradient OR Velocity_Gradient. The bridge consumes
    // F directly (via polar decomposition + log(V)) in BOTH paths -- the
    // element-supplied `rVariables.StrainVector` is bypassed entirely for
    // the kernel-driving Δε. Round 2 (2026-05-11) added explicit MPM
    // support: Check no longer rejects `IS_EXPLICIT`, and
    // `CalculateMaterialResponseCauchy` commits internal state at the end
    // of the call when explicit (because the element never calls Finalize
    // in explicit -- .cpp:882-883).
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Deformation_Gradient; }
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }

    void GetLawFeatures(Features& rFeatures) override;

    bool Has(const Variable<double>& rThisVariable) override;
    bool Has(const Variable<Vector>& rThisVariable) override;

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override;
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override;

    void SetValue(const Variable<double>& rThisVariable,
                  const double& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;
    void SetValue(const Variable<Vector>& rThisVariable,
                  const Vector& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    // Override ResetMaterial -- the base ConstitutiveLaw::ResetMaterial
    // throws KRATOS_ERROR ("Calling virtual function for ResetMaterial"),
    // so any caller (MPM ResetConstitutiveLaw, restart workflows, user
    // utilities) that hits the unoverridden default would crash. Clears
    // all committed bridge history and leaves the law in a state
    // equivalent to a freshly-initialized one. Codex 2026-05-11 F12.
    void ResetMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues) override;

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) const override;

private:
    // Persisted committed state. All vectors stored in Abaqus Voigt order
    // [11,22,33,12,13,23] -- the kernel-native order. Kratos exposure
    // applies the slot 4<->5 swap on the boundary only.
    std::array<double, 6> mSigEffAbq{};      // effective Cauchy, tension+
    std::array<double, 7> mStateAbq{};       // del[0..5] + void
    double mPore = 0.0;                       // compression-positive
    double mDtsubPersistent = 0.0;            // umat-style persisted dtsub
    bool mInitialized = false;                // InitializeMaterial done
    bool mFirstCallStateLoaded = false;       // INITIAL_STRESS_VECTOR + Bauer void done

    // Hughes-Winget corotational bookkeeping (Codex 2026-05-11 F13 fix).
    // mF0_3x3 is the total deformation gradient committed at the END of
    // the previous step. At first call it must equal Identity(3,3). Each
    // step computes the relative deformation gradient
    // `F_rel = F_total_now * F0_3x3^-1`, polar-decomposes
    // `F_rel = R * U`, uses `R` to rotate the stored stress and
    // intergranular strain into the current/end-of-step spatial frame,
    // and uses `log(U)` as the corotational strain increment Δε fed to
    // the kernel. Kratos MPM does NO host-side stress rotation -- the
    // contract puts objectivity inside the law (verified against
    // mpm_updated_lagrangian.cpp:568,926: stress passed in/out unchanged).
    Matrix mF0_3x3;

    // Kratos<->Abaqus Voigt slot 4<->5 swap helpers. Vectors swap idx 4 and
    // 5; the 6x6 tangent swaps both rows 4<->5 and columns 4<->5.
    static void SwapVectorKratosToAbq(const Vector& kratos6,
                                      std::array<double, 6>& abq6);
    static void SwapVectorAbqToKratos(const std::array<double, 6>& abq6,
                                      Vector& kratos6);
    static void SwapTangentAbqToKratos(const std::array<std::array<double, 6>, 6>& abq6x6,
                                       Matrix& kratos6x6);

    // Total <-> effective stress conversion on the three normal slots only
    // (idx 0,1,2 in Abaqus order). Shear slots stay invariant under pore.
    // Convention: sig_eff = sig_total + pore * I_normal (sig tension+, pore
    // compression+); the inverse subtracts.
    static void TotalToEffectiveAbq(std::array<double, 6>& sigAbq, double pore);
    static void EffectiveToTotalAbq(std::array<double, 6>& sigAbq, double pore);

    // Bauer 1996 void-ratio init from total-stress trace. Mirrors Fortran
    // umat .for:171-185 logic on `props(16)`. Returns the initialized
    // void ratio.
    static double InitVoidRatioBauer(const Vector& props16,
                                     const std::array<double, 6>& sigTotalAbq);

    // Compute spatial Almansi strain from a total deformation gradient
    // (MPM passes FT = F * F0 via SetDeformationGradientF, .cpp:255).
    // Returns Kratos-Voigt order [11, 22, 33, 12, 23, 13] with engineering
    // shear convention on off-diagonals (gamma_ij = 2 * e_ij). This is the
    // canonical MPM strain measure -- the element's CalculateKinematics
    // does NOT populate rVariables.StrainVector itself, so the bridge owns
    // computing it (Codex 2026-05-11 F11 high finding).
    //
    // Used ONLY for the strain write-back into rValues.GetStrainVector()
    // so MPM's FinalizeStepVariables persists a meaningful Almansi into
    // `mMP.almansi_strain_vector` for post-processing. The kernel-driving
    // strain increment uses log(U) from polar decomposition instead -- see
    // PolarDecompFAndLogU below.
    //
    // Formula: e = 0.5 * (I - b^-1), b = F * F^T.
    static void CalculateAlmansiFromF(const Matrix& F,
                                      std::array<double, 6>& e_kratos);

    // Polar decomposition F = R * U with R orthogonal, U symmetric positive
    // definite. Implemented via C = F^T F -> eigendecomposition
    // C = B * Λ * B^T, then U = B √Λ B^T, U^-1 = B (1/√Λ) B^T, R = F U^-1.
    // Returns the Hencky strain of the relative stretch in the END-of-step
    // (current / left-stretch) basis: `log(V) = R · log(U) · R^T`, as a
    // Voigt-strain in Abaqus order (engineering shear convention,
    // off-diagonals doubled). Returning log(V) instead of log(U) keeps the
    // strain increment in the SAME basis as the pre-rotated stress and
    // intergranular state fed to the kernel (Codex 2026-05-11 F15 fix:
    // mixing log(U) with R-rotated σ/δ silently corrupts the constitutive
    // update under combined rotation+stretch).
    static void PolarDecompFAndLogU(const Matrix& F_rel,
                                    Matrix& R_out,
                                    std::array<double, 6>& logV_voigt_abq);

    // Rotate a 6-Voigt symmetric STRESS tensor (Abaqus order) by R:
    // σ' = R · σ · R^T. Stress off-diagonals are tensorial -- no engineering
    // factor.
    static void RotateStressVoigtAbq(const std::array<double, 6>& sig_in,
                                     const Matrix& R,
                                     std::array<double, 6>& sig_out);

    // Rotate a 6-Voigt symmetric STRAIN tensor (Abaqus order, engineering
    // shear convention) by R: ε' = R · ε · R^T. Off-diagonals are doubled
    // on the way in (engineering -> tensorial) and doubled again on the
    // way out (tensorial -> engineering).
    static void RotateStrainVoigtAbq(const std::array<double, 6>& eps_in,
                                     const Matrix& R,
                                     std::array<double, 6>& eps_out);

    // bulk_w * 1\otimes1 added to the normal-normal 3x3 sub-block of the
    // Kratos tangent. Operates AFTER the slot 4<->5 swap, in Kratos
    // coordinates -- normal slots [0..2] are the same in both orderings.
    static void AddBulkWNormalBlock(Matrix& C_kratos, double bulk_w);

    // Map kernel `failure_class` to Kratos behavior: Ok -> commit,
    // RkfReject -> set pnewdt for outer cutback (caller-driven), and
    // ContractViolation / Fatal -> KRATOS_ERROR. Inlined in the .cpp.

    // Disambiguated copies of the 15-slot Path-B-side statev for
    // SetValue/GetValue access (kernel statev only carries 7 slots; the
    // bridge owns pore + dtsub_persistent + the spare slots).
    // Empty for Round 1; placeholder hook if external statev IO is wired
    // later.

    friend class Serializer;
    // Serialize every committed history member. Restart, partition
    // migration, and any Kratos Serializer round trip must preserve the
    // bridge's path-dependent state -- it cannot be reconstructed from
    // materials JSON + INITIAL_STRESS_VECTOR alone once loading has
    // begun, because that only captures the initial-configuration
    // snapshot (Codex 2026-05-11 F5 high finding). Order is fixed by the
    // tag strings -- do not reorder save/load entries asymmetrically.
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("mSigEffAbq",            mSigEffAbq);
        rSerializer.save("mStateAbq",             mStateAbq);
        rSerializer.save("mPore",                 mPore);
        rSerializer.save("mDtsubPersistent",      mDtsubPersistent);
        rSerializer.save("mInitialized",          mInitialized);
        rSerializer.save("mFirstCallStateLoaded", mFirstCallStateLoaded);
        rSerializer.save("mF0_3x3",               mF0_3x3);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("mSigEffAbq",            mSigEffAbq);
        rSerializer.load("mStateAbq",             mStateAbq);
        rSerializer.load("mPore",                 mPore);
        rSerializer.load("mDtsubPersistent",      mDtsubPersistent);
        rSerializer.load("mInitialized",          mInitialized);
        rSerializer.load("mFirstCallStateLoaded", mFirstCallStateLoaded);
        rSerializer.load("mF0_3x3",               mF0_3x3);
    }
};

} // namespace Kratos

// Path B - SandHypoplasticCppLaw2DPlaneStrain (round 4, 2026-05-12)
//
// 2D plane-strain interface for Kratos MPM 2D, wrapping the 3D
// SandHypoplasticCppLaw bridge. Composition pattern: holds an internal
// 3D law instance and delegates Calculate / Finalize / state to it,
// padding 2D inputs to 3D and reducing 3D outputs to 2D at the Kratos
// boundary.
//
// Motivation (granular_flow_3D_pathb validation findings, 2026-05-12):
// Kratos MPM 3D8N + sand hypoplastic (both Path A and Path B) has an
// upstream Kratos C++ crash around sim time t ≈ 0.04-0.06s in gravity-
// driven column collapse, while Kratos MPM 2D + Mohr-Coulomb plane-
// strain runs cleanly to t=0.1s+. To run a real Bui-style column
// collapse with the sand hypoplastic kernel, we therefore expose the
// kernel through the mature Kratos MPM 2D plane-strain element pipeline
// instead of 3D8N.
//
// Plane-strain mapping:
//   2D Voigt order (Kratos):  [εxx, εyy, γxy]                      (StrainSize = 3)
//   3D Voigt order (Kratos):  [εxx, εyy, εzz, γxy, γyz, γxz]       (StrainSize = 6)
//   2D index [0,1,2]  ↔  3D index [0,1,3]
//   Plane strain enforces 3D index [2] = 0 (εzz), [4] = 0 (γyz), [5] = 0 (γxz).
//   3D F is built as
//       F_3D = [[F00, F01, 0],
//               [F10, F11, 0],
//               [0,    0,  1]]
//   from the 2x2 plane-strain F.
//
// Implementation note: the kernel + 3D bridge are reused unchanged.
// The 2D wrapper only marshals values at the Kratos `Parameters`
// interface. State (sig_eff, q, pore, dtsub, F0, mEpsPrevExplicitKratos)
// lives inside the internal 3D law instance; we do not duplicate it.
#pragma once

#include <array>

#include "includes/constitutive_law.h"
#include "sand_hypoplastic_application_variables.h"

#include "sand_hypoplastic_cpp_law.h"   // we compose this

namespace Kratos
{

class KRATOS_API(SAND_HYPOPLASTIC_APPLICATION) SandHypoplasticCppLaw2DPlaneStrain
    : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SandHypoplasticCppLaw2DPlaneStrain);

    SandHypoplasticCppLaw2DPlaneStrain();
    SandHypoplasticCppLaw2DPlaneStrain(const SandHypoplasticCppLaw2DPlaneStrain& rOther);
    ~SandHypoplasticCppLaw2DPlaneStrain() override = default;

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override { return 2; }
    SizeType GetStrainSize() const override { return 3; }

    // Same kinematic contract as the 3D bridge (Deformation_Gradient
    // measure, Cauchy stress), just one dimension smaller. Kratos MPM 2D
    // plane-strain element passes a 2x2 F and 3-component Voigt strain.
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Deformation_Gradient; }
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }

    void GetLawFeatures(Features& rFeatures) override;

    // Has / GetValue / SetValue: delegate to the inner 3D law. INTERNAL_VARIABLES
    // exposes the same 15-slot statev as the 3D bridge.
    bool Has(const Variable<double>& rThisVariable) override { return m3DLaw.Has(rThisVariable); }
    bool Has(const Variable<Vector>& rThisVariable) override { return m3DLaw.Has(rThisVariable); }

    double& GetValue(const Variable<double>& rThisVariable, double& rValue) override
    {
        return m3DLaw.GetValue(rThisVariable, rValue);
    }
    Vector& GetValue(const Variable<Vector>& rThisVariable, Vector& rValue) override
    {
        return m3DLaw.GetValue(rThisVariable, rValue);
    }

    void SetValue(const Variable<double>& rThisVariable,
                  const double& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override
    {
        m3DLaw.SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
    void SetValue(const Variable<Vector>& rThisVariable,
                  const Vector& rValue,
                  const ProcessInfo& rCurrentProcessInfo) override
    {
        m3DLaw.SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override
    {
        m3DLaw.InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
    }

    void ResetMaterial(const Properties& rMaterialProperties,
                       const GeometryType& rElementGeometry,
                       const Vector& rShapeFunctionsValues) override
    {
        m3DLaw.ResetMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
    }

    void CalculateMaterialResponseCauchy(Parameters& rValues) override;
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    int Check(const Properties& rMaterialProperties,
              const GeometryType& rElementGeometry,
              const ProcessInfo& rCurrentProcessInfo) const override
    {
        return m3DLaw.Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    }

private:
    // The 3D bridge holds all the constitutive + corotational +
    // failure-class machinery. We just adapt its boundary.
    SandHypoplasticCppLaw m3DLaw;

    // Run the inner 3D law with padded inputs / unpadded outputs.
    // Shared by CalculateMaterialResponseCauchy and FinalizeMaterialResponseCauchy
    // (both delegate to the 3D law with the same marshalling). `is_finalize`
    // selects which 3D entry-point to call.
    void RunInner3D(Parameters& rValues_2D, bool is_finalize);

    // 2D Voigt index -> 3D Voigt index (Kratos order).
    //   2D [0,1,2] = [εxx, εyy, γxy]
    //   3D [0,1,2,3,4,5] = [εxx, εyy, εzz, γxy, γyz, γxz]
    // So 2D[0]→3D[0], 2D[1]→3D[1], 2D[2]→3D[3].
    static constexpr std::array<std::size_t, 3> k2Dto3D = {0, 1, 3};

    friend class Serializer;
    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("m3DLaw", m3DLaw);
    }
    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("m3DLaw", m3DLaw);
    }
};

} // namespace Kratos

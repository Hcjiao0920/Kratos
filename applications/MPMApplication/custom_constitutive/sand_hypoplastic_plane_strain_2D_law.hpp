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

// SandHypoplasticPlaneStrain2DLaw
//
// 2D plane-strain interface for Kratos MPM 2D, wrapping the 3D
// SandHypoplastic3DLaw (hypoplastic sand model after von Wolffersdorff
// 1996, with intergranular strain extension after Niemunis & Herle
// 1997; see the 3D law for the full model documentation). Composition
// pattern: holds an internal 3D law instance and delegates Calculate /
// Finalize / state access to it, padding 2D inputs to 3D and reducing
// 3D outputs to 2D at the Kratos boundary.
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
// Implementation note: the constitutive kernel and the 3D law are
// reused unchanged; this wrapper only marshals values at the Kratos
// `Parameters` interface. State (sig_eff, q, pore, dtsub, F0,
// mEpsPrevExplicitKratos) lives inside the internal 3D law instance
// and is not duplicated here.

#if !defined (KRATOS_SAND_HYPOPLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define       KRATOS_SAND_HYPOPLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED

#include <array>

#include "includes/constitutive_law.h"
#include "mpm_application_variables.h"

#include "sand_hypoplastic_3D_law.hpp"   // composed 3D law

namespace Kratos
{

class KRATOS_API(MPM_APPLICATION) SandHypoplasticPlaneStrain2DLaw
    : public ConstitutiveLaw
{
public:
    using BaseType = ConstitutiveLaw;
    using SizeType = std::size_t;

    KRATOS_CLASS_POINTER_DEFINITION(SandHypoplasticPlaneStrain2DLaw);

    SandHypoplasticPlaneStrain2DLaw();
    SandHypoplasticPlaneStrain2DLaw(const SandHypoplasticPlaneStrain2DLaw& rOther);
    ~SandHypoplasticPlaneStrain2DLaw() override = default;

    ConstitutiveLaw::Pointer Clone() const override;

    SizeType WorkingSpaceDimension() override { return 2; }
    SizeType GetStrainSize() const override { return 3; }

    // Same kinematic contract as the 3D law (Deformation_Gradient
    // measure, Cauchy stress), just one dimension smaller. Kratos MPM 2D
    // plane-strain element passes a 2x2 F and 3-component Voigt strain.
    StrainMeasure GetStrainMeasure() override { return StrainMeasure_Deformation_Gradient; }
    StressMeasure GetStressMeasure() override { return StressMeasure_Cauchy; }

    void GetLawFeatures(Features& rFeatures) override;

    // Has / GetValue / SetValue: delegate to the inner 3D law. INTERNAL_VARIABLES
    // exposes the same 15-slot state-variable vector as the 3D law.
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
    // The inner 3D law holds all the constitutive, corotational and
    // failure-handling machinery; this wrapper only adapts its boundary.
    SandHypoplastic3DLaw m3DLaw;

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

#endif // KRATOS_SAND_HYPOPLASTIC_PLANE_STRAIN_2D_LAW_H_INCLUDED

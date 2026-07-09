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

// SandHypoplasticPlaneStrain2DLaw (implementation).
// See header for design rationale. Composition wrapper that pads 2D
// plane-strain inputs to 3D, calls the inner 3D law, reduces 3D outputs
// to 2D.
#include "sand_hypoplastic_plane_strain_2D_law.hpp"

#include <cmath>

#include "includes/checks.h"
#include "includes/mat_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{

constexpr std::array<std::size_t, 3> SandHypoplasticPlaneStrain2DLaw::k2Dto3D;

//********************************************************************
// Lifecycle / Clone
//********************************************************************

SandHypoplasticPlaneStrain2DLaw::SandHypoplasticPlaneStrain2DLaw()
    : ConstitutiveLaw(), m3DLaw() {}

SandHypoplasticPlaneStrain2DLaw::SandHypoplasticPlaneStrain2DLaw(const SandHypoplasticPlaneStrain2DLaw& rOther)
    : ConstitutiveLaw(rOther), m3DLaw(rOther.m3DLaw) {}

ConstitutiveLaw::Pointer SandHypoplasticPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<SandHypoplasticPlaneStrain2DLaw>(*this);
}

//********************************************************************
// Features / metadata
//********************************************************************

void SandHypoplasticPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(FINITE_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    rFeatures.mStrainSize     = static_cast<SizeType>(3);
    rFeatures.mSpaceDimension = static_cast<SizeType>(WorkingSpaceDimension());
}

//********************************************************************
// Calculate / Finalize (delegate to inner 3D law with pad / unpad)
//********************************************************************

void SandHypoplasticPlaneStrain2DLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    RunInner3D(rValues, /*is_finalize=*/false);
}

void SandHypoplasticPlaneStrain2DLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    RunInner3D(rValues, /*is_finalize=*/true);
}

void SandHypoplasticPlaneStrain2DLaw::RunInner3D(Parameters& rValues_2D, bool is_finalize)
{
    constexpr std::size_t kVoigt2D = 3;
    constexpr std::size_t kVoigt3D = 6;

    // ----- 1. Pad F: 2x2 plane-strain -> 3x3 with F33=1 ---------------
    const Matrix& F2D = rValues_2D.GetDeformationGradientF();
    KRATOS_ERROR_IF(F2D.size1() != 2 || F2D.size2() != 2)
        << "SandHypoplasticPlaneStrain2DLaw: expected 2x2 F from Kratos MPM 2D, got "
        << F2D.size1() << "x" << F2D.size2();

    Matrix F3D = ZeroMatrix(3, 3);
    F3D(0, 0) = F2D(0, 0); F3D(0, 1) = F2D(0, 1);
    F3D(1, 0) = F2D(1, 0); F3D(1, 1) = F2D(1, 1);
    F3D(2, 2) = 1.0; // plane strain: no z deformation

    const double detF = rValues_2D.GetDeterminantF(); // = det(F2D) = det(F3D) under plane strain

    // ----- 2. Pad strain: 3-vector -> 6-vector with zeros in zz/yz/xz -
    const Vector& strain2D = rValues_2D.GetStrainVector();
    KRATOS_ERROR_IF(strain2D.size() != kVoigt2D)
        << "SandHypoplasticPlaneStrain2DLaw: expected 3-component 2D Voigt strain, got "
        << strain2D.size();
    Vector strain3D = ZeroVector(kVoigt3D);
    for (std::size_t i = 0; i < kVoigt2D; ++i) strain3D[k2Dto3D[i]] = strain2D[i];
    // strain3D[2] = 0  (εzz, plane-strain)
    // strain3D[4] = 0  (γyz)
    // strain3D[5] = 0  (γxz)

    // ----- 3. Scratch 6-vector stress + 6x6 tangent -------------------
    Vector stress3D = ZeroVector(kVoigt3D);
    Matrix C3D      = ZeroMatrix(kVoigt3D, kVoigt3D);

    // ----- 4. Build inner 3D Parameters object ------------------------
    ConstitutiveLaw::Parameters params_3D;
    params_3D.SetMaterialProperties(rValues_2D.GetMaterialProperties());
    params_3D.SetProcessInfo(rValues_2D.GetProcessInfo());
    params_3D.SetDeformationGradientF(F3D);
    params_3D.SetDeterminantF(detF);
    params_3D.SetStrainVector(strain3D);
    params_3D.SetStressVector(stress3D);
    params_3D.SetConstitutiveMatrix(C3D);
    params_3D.SetOptions(rValues_2D.GetOptions());

    // Shape functions / DN_DX: forward if set on the caller's params.
    // Kratos MPM elements typically don't read these inside the
    // constitutive law (the inner 3D law ignores them), but pass-through
    // is free and makes the inner Parameters object structurally
    // equivalent.
    if (rValues_2D.IsSetShapeFunctionsValues()) {
        params_3D.SetShapeFunctionsValues(rValues_2D.GetShapeFunctionsValues());
    }
    if (rValues_2D.IsSetShapeFunctionsDerivatives()) {
        params_3D.SetShapeFunctionsDerivatives(rValues_2D.GetShapeFunctionsDerivatives());
    }

    // ----- 5. Call inner 3D law --------------------------------------
    if (is_finalize) {
        m3DLaw.FinalizeMaterialResponseCauchy(params_3D);
    } else {
        m3DLaw.CalculateMaterialResponseCauchy(params_3D);
    }

    // ----- 6. Unpad outputs back to 2D --------------------------------
    const Flags& options = rValues_2D.GetOptions();
    const bool compute_stress  = options.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool compute_tangent = options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    if (compute_stress && rValues_2D.IsSetStressVector()) {
        Vector& stress_out_2D = rValues_2D.GetStressVector();
        if (stress_out_2D.size() != kVoigt2D) stress_out_2D.resize(kVoigt2D, false);
        // Extract σxx, σyy, σxy from 3D σ. σzz, σyz, σxz are computed by
        // the 3D kernel under plane-strain constraint but not exposed to
        // the 2D Kratos element (it has no slot for them).
        for (std::size_t i = 0; i < kVoigt2D; ++i) stress_out_2D[i] = stress3D[k2Dto3D[i]];
    }

    if (compute_tangent && rValues_2D.IsSetConstitutiveMatrix()) {
        Matrix& C_out_2D = rValues_2D.GetConstitutiveMatrix();
        if (C_out_2D.size1() != kVoigt2D || C_out_2D.size2() != kVoigt2D) {
            C_out_2D.resize(kVoigt2D, kVoigt2D, false);
        }
        // Plane-strain reduction of 6x6 tangent: pick rows and columns
        // corresponding to 2D Voigt slots ([0, 1, 3] of 3D order). This
        // is the standard plane-strain Schur reduction when εzz=0 is
        // an enforced kinematic constraint and σzz is whatever the law
        // produces.
        for (std::size_t i = 0; i < kVoigt2D; ++i)
            for (std::size_t j = 0; j < kVoigt2D; ++j)
                C_out_2D(i, j) = C3D(k2Dto3D[i], k2Dto3D[j]);
    }

    // Strain write-back. The 3D law writes Almansi-from-F-3D into
    // strain3D in its IMPLICIT path (for MPM's mMP.almansi_strain_vector
    // persistence) and skips the write-back in its EXPLICIT path. Under
    // plane strain the 3D Almansi has zero zz/yz/xz channels by
    // construction (F_3D has F33=1, F13=F23=F31=F32=0), so just
    // truncate to 3 components for the 2D caller. If the 3D path
    // didn't touch strain3D (explicit), this re-copies the same value
    // that was padded in -- a no-op for the caller's state.
    if (rValues_2D.IsSetStrainVector()) {
        Vector& strain_out_2D = rValues_2D.GetStrainVector();
        if (strain_out_2D.size() != kVoigt2D) strain_out_2D.resize(kVoigt2D, false);
        for (std::size_t i = 0; i < kVoigt2D; ++i) strain_out_2D[i] = strain3D[k2Dto3D[i]];
    }
}

} // namespace Kratos

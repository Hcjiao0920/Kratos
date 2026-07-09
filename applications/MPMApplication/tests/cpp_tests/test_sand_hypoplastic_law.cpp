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

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "testing/testing.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

// Application includes
#include "mpm_application_variables.h"

// Material laws
#include "custom_constitutive/sand_hypoplastic_3D_law.hpp"
#include "custom_constitutive/sand_hypoplastic_plane_strain_2D_law.hpp"

namespace Kratos
{
namespace Testing
{
    typedef Node NodeType;

    namespace
    {
        // Hostun dense sand calibration (Herle & Gudehus 1999), the common
        // reference parameter set in the hypoplasticity literature. Order:
        // [phi_deg, p_t, h_s, n, e_d0, e_c0, e_i0, alpha, beta,
        //  m_R, m_T, r_uc, beta_r, chi, bulk_w, e_0]
        Vector HostunDenseParameters()
        {
            Vector props(16);
            props[0]  = 34.0;    // phi (deg)
            props[1]  = 1.0;     // p_t (kPa)
            props[2]  = 1.0e6;   // h_s (kPa)
            props[3]  = 0.29;    // n
            props[4]  = 0.61;    // e_d0
            props[5]  = 0.96;    // e_c0
            props[6]  = 1.09;    // e_i0
            props[7]  = 0.13;    // alpha
            props[8]  = 2.00;    // beta
            props[9]  = 5.0;     // m_R
            props[10] = 2.0;     // m_T
            props[11] = 1.0e-4;  // r_uc
            props[12] = 0.5;     // beta_r
            props[13] = 6.0;     // chi
            props[14] = 0.0;     // bulk_w (dry)
            props[15] = 0.65;    // e_0 (dense)
            return props;
        }

        void SetIsotropicInitialStress(Properties& rProps, const double MeanStress)
        {
            Vector initial_stress = ZeroVector(6);
            initial_stress[0] = MeanStress;
            initial_stress[1] = MeanStress;
            initial_stress[2] = MeanStress;
            rProps.SetValue(INITIAL_STRESS_VECTOR, initial_stress);
        }

        // One 3D constitutive call (Calculate + Finalize) driven by a total
        // deformation gradient, mirroring the implicit MPM element protocol.
        void RunOneStep3D(SandHypoplastic3DLaw& rLaw,
                          Properties& rProps,
                          ModelPart& rModelPart,
                          const Matrix& rF,
                          Vector& rStress,
                          Matrix& rTangent)
        {
            ConstitutiveLaw::Parameters cl_parameters;
            Vector strain_vector = ZeroVector(6);
            rStress = ZeroVector(6);
            rTangent = ZeroMatrix(6, 6);

            Flags& r_options = cl_parameters.GetOptions();
            r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

            cl_parameters.SetProcessInfo(rModelPart.GetProcessInfo());
            cl_parameters.SetMaterialProperties(rProps);
            cl_parameters.SetStrainVector(strain_vector);
            cl_parameters.SetStressVector(rStress);
            cl_parameters.SetConstitutiveMatrix(rTangent);
            cl_parameters.SetDeformationGradientF(rF);
            cl_parameters.SetDeterminantF(MathUtils<double>::Det(rF));

            rLaw.CalculateMaterialResponseCauchy(cl_parameters);
            rLaw.FinalizeMaterialResponseCauchy(cl_parameters);
        }
    }

    // Isotropic compression step from p = 100 kPa: checks the barotropic
    // stiffness of the model (reference values from the validated
    // implementation, cross-checked at machine precision against the
    // reference Fortran UMAT).
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawSandHypoplasticIsotropicCompression, KratosMPMFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        r_model_part.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

        Properties material_properties;
        material_properties.SetValue(SAND_HYPOPLASTIC_PARAMETERS, HostunDenseParameters());
        SetIsotropicInitialStress(material_properties, -100.0);

        SandHypoplastic3DLaw cl;
        Vector dummy;
        Geometry<NodeType> dummy_geometry;
        cl.InitializeMaterial(material_properties, dummy_geometry, dummy);

        // Volumetric compression, eps_v = -1e-4 in one step.
        const double e = -1.0e-4 / 3.0;
        Matrix F = IdentityMatrix(3);
        F(0, 0) = 1.0 + e;
        F(1, 1) = 1.0 + e;
        F(2, 2) = 1.0 + e;

        Vector stress;
        Matrix tangent;
        RunOneStep3D(cl, material_properties, r_model_part, F, stress, tangent);

        const double tol_stress  = 5.0e-2;  // kPa; ~4e-4 relative
        const double tol_tangent = 5.0;     // kPa; ~2e-5 relative
        KRATOS_EXPECT_NEAR(stress[0], -120.95504747489647, tol_stress);
        KRATOS_EXPECT_NEAR(stress[1], -120.95504747489647, tol_stress);
        KRATOS_EXPECT_NEAR(stress[2], -120.95504747489647, tol_stress);
        KRATOS_EXPECT_NEAR(stress[3], 0.0, 1.0e-6);
        KRATOS_EXPECT_NEAR(stress[4], 0.0, 1.0e-6);
        KRATOS_EXPECT_NEAR(stress[5], 0.0, 1.0e-6);
        KRATOS_EXPECT_NEAR(tangent(0, 0), 312714.66968642984, tol_tangent);
        KRATOS_EXPECT_NEAR(tangent(0, 1), 138387.60405120772, tol_tangent);
    }

    // Simple shear in the 1-3 plane: discriminates the Kratos
    // [11,22,33,12,23,13] vs Abaqus [11,22,33,12,13,23] Voigt orders. A
    // correct implementation loads sigma_13 (Kratos slot 5) and leaves
    // sigma_23 (Kratos slot 4) at zero.
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawSandHypoplasticSimpleShearVoigtOrder, KratosMPMFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        r_model_part.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

        Properties material_properties;
        material_properties.SetValue(SAND_HYPOPLASTIC_PARAMETERS, HostunDenseParameters());
        SetIsotropicInitialStress(material_properties, -100.0);

        SandHypoplastic3DLaw cl;
        Vector dummy;
        Geometry<NodeType> dummy_geometry;
        cl.InitializeMaterial(material_properties, dummy_geometry, dummy);

        Matrix F = IdentityMatrix(3);
        F(0, 2) = 2.0e-4;  // simple shear gamma_13

        Vector stress;
        Matrix tangent;
        RunOneStep3D(cl, material_properties, r_model_part, F, stress, tangent);

        const double tol_stress = 5.0e-2;
        KRATOS_EXPECT_NEAR(stress[0], -102.22774978721847, tol_stress);
        KRATOS_EXPECT_NEAR(stress[1], -102.22941662438211, tol_stress);
        KRATOS_EXPECT_NEAR(stress[2], -102.23108346155178, tol_stress);
        KRATOS_EXPECT_NEAR(stress[3], 0.0, 1.0e-6);                       // sigma_12
        KRATOS_EXPECT_NEAR(stress[4], 0.0, 1.0e-6);                       // sigma_23
        KRATOS_EXPECT_NEAR(stress[5], 16.668371658019943, tol_stress);    // sigma_13
    }

    // Free-surface / tensile regime: from a zero-stress state, a small
    // tensile step must route through the elastic fallback and return an
    // exactly zero stress (low-mean-pressure regularization), never NaN or
    // a spurious deviator.
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawSandHypoplasticTensileFallback, KratosMPMFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        r_model_part.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

        Properties material_properties;
        material_properties.SetValue(SAND_HYPOPLASTIC_PARAMETERS, HostunDenseParameters());
        SetIsotropicInitialStress(material_properties, 0.0);

        SandHypoplastic3DLaw cl;
        Vector dummy;
        Geometry<NodeType> dummy_geometry;
        cl.InitializeMaterial(material_properties, dummy_geometry, dummy);

        Matrix F = IdentityMatrix(3);
        F(0, 0) = 1.0 + 1.0e-5;
        F(1, 1) = 1.0 + 1.0e-5;
        F(2, 2) = 1.0 + 1.0e-5;

        Vector stress;
        Matrix tangent;
        RunOneStep3D(cl, material_properties, r_model_part, F, stress, tangent);

        for (std::size_t i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(stress[i], 0.0, 1.0e-12);
        }
        // The fallback still provides a finite, positive-definite-ish tangent.
        KRATOS_EXPECT_NEAR(tangent(0, 0), 878.37837837837765, 1.0e-3);
    }

    // The plane-strain wrapper must reproduce the 3D law exactly on an
    // in-plane path (same-run comparison; platform independent).
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawSandHypoplasticPlaneStrainMatches3D, KratosMPMFastSuite)
    {
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");
        r_model_part.GetProcessInfo().SetValue(DELTA_TIME, 1.0);

        Properties material_properties;
        material_properties.SetValue(SAND_HYPOPLASTIC_PARAMETERS, HostunDenseParameters());
        SetIsotropicInitialStress(material_properties, -100.0);

        Vector dummy;
        Geometry<NodeType> dummy_geometry;

        // 3D reference: F = diag(1 - 1e-4, 1, 1).
        SandHypoplastic3DLaw cl3d;
        cl3d.InitializeMaterial(material_properties, dummy_geometry, dummy);
        Matrix F3 = IdentityMatrix(3);
        F3(0, 0) = 1.0 - 1.0e-4;
        Vector stress3;
        Matrix tangent3;
        RunOneStep3D(cl3d, material_properties, r_model_part, F3, stress3, tangent3);

        // 2D plane-strain twin: F = [[1 - 1e-4, 0], [0, 1]].
        SandHypoplasticPlaneStrain2DLaw cl2d;
        cl2d.InitializeMaterial(material_properties, dummy_geometry, dummy);

        ConstitutiveLaw::Parameters cl_parameters;
        Vector strain2 = ZeroVector(3);
        Vector stress2 = ZeroVector(3);
        Matrix tangent2 = ZeroMatrix(3, 3);

        Flags& r_options = cl_parameters.GetOptions();
        r_options.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
        r_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);

        cl_parameters.SetProcessInfo(r_model_part.GetProcessInfo());
        cl_parameters.SetMaterialProperties(material_properties);
        cl_parameters.SetStrainVector(strain2);
        cl_parameters.SetStressVector(stress2);
        cl_parameters.SetConstitutiveMatrix(tangent2);

        Matrix F2 = IdentityMatrix(2);
        F2(0, 0) = 1.0 - 1.0e-4;
        cl_parameters.SetDeformationGradientF(F2);
        cl_parameters.SetDeterminantF(1.0 - 1.0e-4);

        cl2d.CalculateMaterialResponseCauchy(cl_parameters);
        cl2d.FinalizeMaterialResponseCauchy(cl_parameters);

        // 2D Voigt [xx, yy, xy] maps to 3D Voigt slots [0, 1, 3].
        KRATOS_EXPECT_NEAR(stress2[0], stress3[0], 1.0e-12);
        KRATOS_EXPECT_NEAR(stress2[1], stress3[1], 1.0e-12);
        KRATOS_EXPECT_NEAR(stress2[2], stress3[3], 1.0e-12);
    }

    // Material input validation happens once at InitializeMaterial: a
    // missing or wrongly sized parameter vector must fail loudly there.
    KRATOS_TEST_CASE_IN_SUITE(MPMConstitutiveLawSandHypoplasticRejectsBadParameters, KratosMPMFastSuite)
    {
        SandHypoplastic3DLaw cl;
        Vector dummy;
        Geometry<NodeType> dummy_geometry;

        Properties no_params;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            cl.InitializeMaterial(no_params, dummy_geometry, dummy),
            "SAND_HYPOPLASTIC_PARAMETERS");

        Properties short_params;
        Vector too_short(4);
        too_short[0] = 34.0; too_short[1] = 1.0; too_short[2] = 1.0e6; too_short[3] = 0.29;
        short_params.SetValue(SAND_HYPOPLASTIC_PARAMETERS, too_short);
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
            cl.InitializeMaterial(short_params, dummy_geometry, dummy),
            "16 entries");
    }

} // namespace Testing
} // namespace Kratos

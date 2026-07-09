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

// Component-level tests of the standalone sand-hypoplasticity kernel
// (parameter validation, RKF-substepped integration driver). The full
// Kratos-facing behaviour is covered by test_sand_hypoplastic_law.cpp;
// these cases pin the kernel's own control-flow contracts, following the
// component-test precedent of the Mohr-Coulomb / Cam-Clay flow rules.

// System includes
#include <cmath>

// External includes

// Project includes
#include "testing/testing.h"

// Application includes
#include "custom_constitutive/sand_hypoplastic_kernel.hpp"

namespace Kratos
{
namespace Testing
{
    namespace
    {
        // Hostun dense sand (Herle & Gudehus 1999). Fortran props(1..16) order.
        SandHypoCpp::Props16 HostunDenseProps()
        {
            return {34.0, 1.0, 1.0e6, 0.29, 0.61, 0.96, 1.09, 0.13, 2.00,
                    5.0, 2.0, 1.0e-4, 0.5, 6.0, 0.0, 0.65};
        }

        // A well-posed compressive state: isotropic effective stress of
        // -100 kPa (tension-positive, Voigt-Abaqus) with zero intergranular
        // strain and a dense void ratio.
        void CompressiveState(double* sig, double* q)
        {
            for (int i = 0; i < 6; ++i) sig[i] = (i < 3) ? -100.0 : 0.0;
            for (int i = 0; i < 6; ++i) q[i] = 0.0;
            q[6] = 0.65;  // void ratio
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMSandHypoplasticKernelCheckParms, KratosMPMFastSuite)
    {
        // The reference calibration validates cleanly.
        auto ok = SandHypoCpp::check_parms(HostunDenseProps());
        KRATOS_EXPECT_EQ(ok.error, 0);
        KRATOS_EXPECT_TRUE(ok.failed == nullptr);
        // phi is converted to radians in the kernel-internal parms.
        const double pi = 3.14159265358979323846;
        KRATOS_EXPECT_NEAR(ok.parms[0], 34.0 * pi / 180.0, 1.0e-12);

        // Non-positive friction angle is rejected and named.
        auto bad_phi_props = HostunDenseProps();
        bad_phi_props[0] = 0.0;
        auto bad_phi = SandHypoCpp::check_parms(bad_phi_props);
        KRATOS_EXPECT_TRUE(bad_phi.error != 0);
        KRATOS_EXPECT_TRUE(bad_phi.failed != nullptr);

        // Intergranular strain enabled (m_R > 0.5) with a non-positive
        // reference length r_uc would divide by zero inside the model;
        // the validator rejects the combination up front.
        auto bad_ruc_props = HostunDenseProps();
        bad_ruc_props[11] = 0.0;
        auto bad_ruc = SandHypoCpp::check_parms(bad_ruc_props);
        KRATOS_EXPECT_TRUE(bad_ruc.error != 0);
        KRATOS_EXPECT_TRUE(bad_ruc.failed != nullptr);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMSandHypoplasticKernelZeroStrainIncrementIsIdentity, KratosMPMFastSuite)
    {
        auto checked = SandHypoCpp::check_parms(HostunDenseProps());
        KRATOS_EXPECT_EQ(checked.error, 0);

        double sig[6], q[7];
        CompressiveState(sig, q);
        const double deps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        const auto r = SandHypoCpp::integrate_step(sig, q, deps, checked.parms,
                                                   /*dtime=*/1.0, /*dtsub_in=*/0.0);
        KRATOS_EXPECT_EQ(r.error, 0);
        for (int i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(r.sig[i], sig[i], 1.0e-10);
        }
        for (int i = 0; i < 7; ++i) {
            KRATOS_EXPECT_NEAR(r.q[i], q[i], 1.0e-12);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMSandHypoplasticKernelCompressionPhysics, KratosMPMFastSuite)
    {
        auto checked = SandHypoCpp::check_parms(HostunDenseProps());
        KRATOS_EXPECT_EQ(checked.error, 0);

        double sig[6], q[7];
        CompressiveState(sig, q);
        // Volumetric compression, eps_v = -1e-4.
        const double e = -1.0e-4 / 3.0;
        const double deps[6] = {e, e, e, 0.0, 0.0, 0.0};

        const auto r = SandHypoCpp::integrate_step(sig, q, deps, checked.parms, 1.0, 0.0);
        KRATOS_EXPECT_EQ(r.error, 0);
        KRATOS_EXPECT_FALSE(r.used_elastic);

        // Mean pressure grows (more compressive normals), stress stays
        // isotropic, void ratio decreases per e_dot = (1 + e) * tr(D).
        for (int i = 0; i < 3; ++i) {
            KRATOS_EXPECT_TRUE(r.sig[i] < -100.0);
        }
        KRATOS_EXPECT_NEAR(r.sig[0], r.sig[1], 1.0e-8);
        KRATOS_EXPECT_NEAR(r.sig[1], r.sig[2], 1.0e-8);
        for (int i = 3; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(r.sig[i], 0.0, 1.0e-10);
        }
        KRATOS_EXPECT_TRUE(r.q[6] < 0.65);
        KRATOS_EXPECT_NEAR(r.q[6], 0.65 + (1.0 + 0.65) * (-1.0e-4), 1.0e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(MPMSandHypoplasticKernelTensileStateRoutesToElasticFallback, KratosMPMFastSuite)
    {
        auto checked = SandHypoCpp::check_parms(HostunDenseProps());
        KRATOS_EXPECT_EQ(checked.error, 0);

        // Physically untenable for a cohesionless sand: tensile principal
        // stresses. The admissibility gate must route the step to the
        // linear-elastic fallback instead of running the main model on it.
        double sig[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.65};
        const double e = 1.0e-6;
        const double deps[6] = {e, e, e, 0.0, 0.0, 0.0};

        const auto r = SandHypoCpp::integrate_step(sig, q, deps, checked.parms, 1.0, 0.0);
        KRATOS_EXPECT_EQ(r.error, 0);
        KRATOS_EXPECT_TRUE(r.used_elastic);
        for (int i = 0; i < 6; ++i) {
            KRATOS_EXPECT_TRUE(std::isfinite(r.sig[i]));
        }
        // The intergranular strain state is frozen on the fallback path.
        for (int i = 0; i < 6; ++i) {
            KRATOS_EXPECT_NEAR(r.q[i], q[i], 1.0e-12);
        }
    }

} // namespace Testing
} // namespace Kratos

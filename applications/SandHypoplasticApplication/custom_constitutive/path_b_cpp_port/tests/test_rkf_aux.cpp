// Standalone unit-test exerciser for the Path B kernel RKF-auxiliary
// routines `calc_elasti` and `norm_res`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// Coverage:
//
//   calc_elasti:
//     (E1) Hydrostatic strain at textbook E=100, nu=0.48 (the
//          inittension fallback constants from .for:329 area).
//          Hand-pinned d_sig and bulk-modulus K = E/(3(1-2nu)).
//     (E2) Pure shear strain (engineering form gamma_12). Tests that
//          the half-factor in II[3][3] correctly produces D[3][3] = G
//          and d_sig[3] = G * gamma.
//     (E3) D structural verification: lambda + 2 mu on normal diagonal,
//          lambda on normal off-diagonal, G on shear diagonal, zero
//          for normal-shear cross terms, full Voigt symmetry.
//     (E4) Defensive deviation: youngel = 0 -> D and d_sig both zero.
//     (E5) Defensive deviation: youngel = -100 (Fortran umat's test=3
//          sentinel) -> D and d_sig both zero. The Fortran original
//          would propagate uninitialized DDtan here.
//
//   norm_res:
//     (N1) Identical y_til == y_hat -> norm_R = 0, no NaN.
//     (N2) Pure stress residual: hand-pinned norm_R from a single
//          slot bump, against the relative-error formula.
//     (N3) Pure void residual: norm_R = del_void / void_hat.
//     (N4) NaN sentinel: NaN in y_hat[0] (a stress slot) -> NaN
//          propagates into norm_sig -> sentinel fires; norm_R = 1e20.
//     (N5) Inf sentinel: +Inf in y_hat[12] (void) -> norm_R = 1e20
//          and nan_detected = true.
//     (N6) Overflow sentinel: 1e31 in y_hat[12] -> norm_R = 1e20.
//     (N7) Zero-stress edge: y_hat sig identically zero -> norm_sig = 0,
//          err[0..5] = 0 (no division by zero, no NaN poisoning).

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <limits>

using Kratos::SandHypoCpp::ElasticIncrementResult;
using Kratos::SandHypoCpp::NormResResult;
using Kratos::SandHypoCpp::calc_elasti;
using Kratos::SandHypoCpp::kYDim;
using Kratos::SandHypoCpp::norm_res;

namespace {

int g_pass = 0;
int g_fail = 0;

void expect_eq_bool(bool actual, bool expected, const char* label) {
    if (actual == expected) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr, "[FAIL] %s: expected %d, got %d\n",
                     label, static_cast<int>(expected), static_cast<int>(actual));
    }
}

void expect_near(double actual, double expected, double tol,
                 const char* label) {
    const double diff = std::abs(actual - expected);
    if (diff <= tol) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr,
                     "[FAIL] %s: expected %.17g, got %.17g (|diff|=%.3g)\n",
                     label, expected, actual, diff);
    }
}

}  // namespace

int main() {
    // ====================================================================
    // calc_elasti
    // ====================================================================

    // (E1) Hydrostatic deps with textbook inittension constants.
    //      E = 100, nu = 0.48.
    //        lambda = E*nu / ((1+nu)(1-2nu))
    //               = 100 * 0.48 / (1.48 * 0.04) = 48 / 0.0592 = 810.81081...
    //        mu     = E / (2(1+nu)) = 100 / 2.96 = 33.78378...
    //        K      = lambda + (2/3)*mu = 810.81 + 22.52 = 833.33...
    //               (verify: K = E/(3(1-2nu)) = 100/(3*0.04) = 833.33 OK)
    //      For deps = (-1e-3, -1e-3, -1e-3, 0, 0, 0):
    //        d_sig[0..2] = (3 lambda + 2 mu) * (-1e-3) = 2500 * (-1e-3) = -2.5
    //        d_sig[3..5] = 0
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        const ElasticIncrementResult r = calc_elasti(deps, 100.0, 0.48);

        expect_near(r.d_sig[0], -2.5, 1e-12, "E1 d_sig[0] == -2.5 (hydrostatic)");
        expect_near(r.d_sig[1], -2.5, 1e-12, "E1 d_sig[1] == -2.5");
        expect_near(r.d_sig[2], -2.5, 1e-12, "E1 d_sig[2] == -2.5");
        expect_near(r.d_sig[3], 0.0,  0.0,   "E1 d_sig[3] == 0 (no shear deps)");
        expect_near(r.d_sig[4], 0.0,  0.0,   "E1 d_sig[4] == 0");
        expect_near(r.d_sig[5], 0.0,  0.0,   "E1 d_sig[5] == 0");

        // Verify D[0][0] = lambda + 2 mu = 100 * 0.52 / (1.48 * 0.04) = 878.378...
        const double E      = 100.0;
        const double nu     = 0.48;
        const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double mu     = E / (2.0 * (1.0 + nu));
        expect_near(r.D[0][0], lambda + 2.0 * mu, 1e-10,
                    "E1 D[0][0] == lambda + 2*mu");
        expect_near(r.D[0][1], lambda, 1e-10,
                    "E1 D[0][1] == lambda");
        expect_near(r.D[3][3], mu, 1e-10,
                    "E1 D[3][3] == mu (= G)");
    }

    // (E2) Pure shear deps. Engineering shear gamma_12 = 1e-3.
    //      d_sig[3] = G * gamma. G = E/(2(1+nu)) = 33.78378...
    //      For gamma = 1e-3: d_sig[3] = 33.78378e-3.
    {
        const double deps[6] = {0.0, 0.0, 0.0, 1.0e-3, 0.0, 0.0};
        const ElasticIncrementResult r = calc_elasti(deps, 100.0, 0.48);

        const double G = 100.0 / (2.0 * 1.48);
        expect_near(r.d_sig[3], G * 1.0e-3, 1e-14,
                    "E2 d_sig[3] == G * gamma_12");
        expect_near(r.d_sig[0], 0.0, 0.0, "E2 d_sig[0] == 0");
        expect_near(r.d_sig[1], 0.0, 0.0, "E2 d_sig[1] == 0");
        expect_near(r.d_sig[2], 0.0, 0.0, "E2 d_sig[2] == 0");
        expect_near(r.d_sig[4], 0.0, 0.0, "E2 d_sig[4] == 0");
        expect_near(r.d_sig[5], 0.0, 0.0, "E2 d_sig[5] == 0");
    }

    // (E3) Structural verification of D for textbook E, nu.
    {
        const double deps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const ElasticIncrementResult r = calc_elasti(deps, 100.0, 0.48);
        const double E      = 100.0;
        const double nu     = 0.48;
        const double lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
        const double mu     = E / (2.0 * (1.0 + nu));

        // Normal diagonal blocks: lambda + 2 mu.
        for (int i = 0; i < 3; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "E3 D[%d][%d] == lambda + 2 mu", i, i);
            expect_near(r.D[i][i], lambda + 2.0 * mu, 1e-10, label);
        }
        // Normal off-diagonal: lambda.
        expect_near(r.D[0][1], lambda, 1e-10, "E3 D[0][1] == lambda");
        expect_near(r.D[1][0], lambda, 1e-10, "E3 D[1][0] == lambda");
        expect_near(r.D[1][2], lambda, 1e-10, "E3 D[1][2] == lambda");
        expect_near(r.D[2][0], lambda, 1e-10, "E3 D[2][0] == lambda");

        // Shear diagonal: G.
        for (int i = 3; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "E3 D[%d][%d] == G", i, i);
            expect_near(r.D[i][i], mu, 1e-10, label);
        }

        // Normal-shear coupling = 0.
        for (int i = 0; i < 3; ++i) {
            for (int j = 3; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "E3 D[%d][%d] == 0 (normal-shear)", i, j);
                expect_near(r.D[i][j], 0.0, 0.0, label);
                std::snprintf(label, sizeof(label),
                              "E3 D[%d][%d] == 0 (shear-normal)", j, i);
                expect_near(r.D[j][i], 0.0, 0.0, label);
            }
        }

        // Shear-shear off-diagonal = 0.
        expect_near(r.D[3][4], 0.0, 0.0, "E3 D[3][4] == 0");
        expect_near(r.D[3][5], 0.0, 0.0, "E3 D[3][5] == 0");
        expect_near(r.D[4][5], 0.0, 0.0, "E3 D[4][5] == 0");

        // Symmetry of D.
        for (int i = 0; i < 6; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "E3 D[%d][%d] == D[%d][%d]", i, j, j, i);
                expect_near(r.D[i][j], r.D[j][i], 0.0, label);
            }
        }
    }

    // (E4) youngel = 0 -> D and d_sig zero (defensive deviation).
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 1.0e-3, 0.0, 0.0};
        const ElasticIncrementResult r = calc_elasti(deps, 0.0, 0.48);
        for (int i = 0; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "E4 d_sig[%d] == 0 (youngel=0)", i);
            expect_near(r.d_sig[i], 0.0, 0.0, label);
        }
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "E4 D[%d][%d] == 0 (youngel=0)", i, j);
                expect_near(r.D[i][j], 0.0, 0.0, label);
            }
        }
    }

    // (E5) youngel = -100 -> D and d_sig zero (the Fortran umat's
    //      test=3 sentinel). Fortran's literal behaviour is to
    //      propagate uninitialized DDtan; Path B refuses.
    {
        const double deps[6] = {-1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0};
        const ElasticIncrementResult r = calc_elasti(deps, -100.0, 0.48);
        for (int i = 0; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "E5 d_sig[%d] == 0 (youngel<0 sentinel)", i);
            expect_near(r.d_sig[i], 0.0, 0.0, label);
        }
        expect_near(r.D[0][0], 0.0, 0.0, "E5 D[0][0] == 0 (youngel<0)");
        expect_near(r.D[3][3], 0.0, 0.0, "E5 D[3][3] == 0 (youngel<0)");
    }

    // ====================================================================
    // norm_res
    // ====================================================================

    auto build_y = [](double sig0, double sig_rest, double q_normal,
                      double q_shear, double void_ratio,
                      double y_out[kYDim]) {
        y_out[0] = sig0;
        y_out[1] = sig_rest;
        y_out[2] = sig_rest;
        y_out[3] = 0.0;
        y_out[4] = 0.0;
        y_out[5] = 0.0;
        y_out[6] = q_normal;
        y_out[7] = q_normal;
        y_out[8] = q_normal;
        y_out[9]  = q_shear;
        y_out[10] = q_shear;
        y_out[11] = q_shear;
        y_out[12] = void_ratio;
    };

    // (N1) Identical y_til == y_hat: norm_R = 0.
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 0.0, 0.0, "N1 identical y -> norm_R == 0");
        expect_eq_bool(r.nan_detected, false,
                       "N1 identical y -> nan_detected == false");
    }

    // (N2) Pure stress residual: y_hat differs from y_til only in slot 0.
    //      sig_hat = (-100, -100, -100, 0, 0, 0); norm_sig (StressLike-dot,
    //      diagonal-only) = sqrt(100^2 + 100^2 + 100^2) = sqrt(30000).
    //      del_sig = (delta, 0, 0, 0, 0, 0). err[0] = delta / sqrt(30000).
    //      All other err = 0 (q identical, void identical).
    //      norm_R = err[0].
    {
        const double delta = 0.1;
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_til[0] = y_hat[0] - delta;  // sig_hat - sig_til = +delta

        const NormResResult r = norm_res(y_til, y_hat);
        const double norm_sig_expected = std::sqrt(3.0 * 100.0 * 100.0);
        const double err0_expected = delta / norm_sig_expected;
        expect_near(r.norm_R, err0_expected, 1e-15,
                    "N2 pure stress residual: norm_R == |delta|/||sig_hat||");
        expect_eq_bool(r.nan_detected, false, "N2 nan_detected == false");
    }

    // (N3) Pure void residual: y_hat differs from y_til only in slot 12.
    //      err[12] = del_void / void_hat. norm_R = err[12].
    {
        const double void_hat = 0.825;
        const double void_til = 0.820;
        const double del_void = void_hat - void_til;
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, void_hat, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_til[12] = void_til;

        const NormResResult r = norm_res(y_til, y_hat);
        const double err12_expected = del_void / void_hat;
        expect_near(r.norm_R, err12_expected, 1e-15,
                    "N3 pure void residual: norm_R == del_void/void_hat");
        expect_eq_bool(r.nan_detected, false, "N3 nan_detected == false");
    }

    // (N4) NaN sentinel: NaN in y_hat[0] -> norm_sig becomes NaN ->
    //      umatisnan_h fires -> norm_R = 1e20, nan_detected = true.
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_hat[0] = std::numeric_limits<double>::quiet_NaN();

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 1.0e20, 0.0,
                    "N4 NaN in y_hat[0] -> norm_R == 1e20 sentinel");
        expect_eq_bool(r.nan_detected, true, "N4 nan_detected == true");
    }

    // (N5) Inf sentinel: +Inf in y_hat[12] -> umatisnan_h fires on void_hat.
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_hat[12] = std::numeric_limits<double>::infinity();

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 1.0e20, 0.0,
                    "N5 +Inf in void_hat -> norm_R == 1e20 sentinel");
        expect_eq_bool(r.nan_detected, true, "N5 nan_detected == true");
    }

    // (N6) Overflow sentinel: finite 1e31 in y_hat[12].
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_hat[12] = 1.0e31;

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 1.0e20, 0.0,
                    "N6 1e31 in void_hat -> norm_R == 1e20 sentinel");
        expect_eq_bool(r.nan_detected, true, "N6 nan_detected == true");
    }

    // (N7) Zero-stress edge: y_hat sig identically zero -> norm_sig = 0.
    //      The Fortran skips the divide; err[0..5] = 0 regardless of
    //      del_sig. Verify no NaN poisoning.
    //
    //      y_hat: sig = 0, q = 0, void = 0.825.
    //      y_til: sig = (0.01, 0.01, 0.01, 0, 0, 0), q identical, void identical.
    //      del_sig != 0 but norm_sig = 0 -> err[0..5] should stay 0.
    //      All other err slots zero too -> norm_R = 0.
    {
        double y_hat[kYDim] = {0.0};
        y_hat[12] = 0.825;
        double y_til[kYDim] = {0.0};
        y_til[0] = 0.01; y_til[1] = 0.01; y_til[2] = 0.01;
        y_til[12] = 0.825;

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 0.0, 0.0,
                    "N7 zero-stress edge: err[0..5] skipped, norm_R == 0");
        expect_eq_bool(r.nan_detected, false,
                       "N7 zero-stress edge: nan_detected == false (no FP poisoning)");
    }

    // (N8) Defensive post-check: finite zero void_hat with non-zero
    //      del_void -> err[12] = del_void/0 = +Inf -> norm_R = Inf.
    //      The Fortran's intermediate-norm scan does not catch this
    //      because void_hat = 0 is finite (passes umatisnan_h). Path B
    //      adds a post-check on norm_R itself; sentinel must fire.
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.0, y_hat);    // void_hat = 0
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_til[12] = -0.005;  // del_void = 0.005 (positive, /0 -> +Inf)

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 1.0e20, 0.0,
                    "N8 void_hat=0 with del_void>0 -> norm_R == 1e20 sentinel");
        expect_eq_bool(r.nan_detected, true,
                       "N8 void_hat=0 with del_void>0 -> nan_detected == true");
    }

    // (N9) Defensive post-check: NaN seeded in y_til (not y_hat).
    //      norm_sig and norm_q (computed from y_hat) stay clean and
    //      pass the intermediate guard, but del_sig[0] = |NaN-finite| =
    //      NaN poisons err[0] and norm_R. Without the post-check the
    //      sentinel would not fire.
    {
        double y_hat[kYDim];
        build_y(-100.0, -100.0, 0.0, 0.0, 0.825, y_hat);
        double y_til[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) y_til[i] = y_hat[i];
        y_til[0] = std::numeric_limits<double>::quiet_NaN();

        const NormResResult r = norm_res(y_til, y_hat);
        expect_near(r.norm_R, 1.0e20, 0.0,
                    "N9 NaN-in-y_til -> norm_R == 1e20 sentinel (post-check)");
        expect_eq_bool(r.nan_detected, true,
                       "N9 NaN-in-y_til -> nan_detected == true");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

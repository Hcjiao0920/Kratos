// Standalone unit-test exerciser for the Path B kernel `get_F_sig_q`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// The function under test is a *composer*: it picks istrain via the
// Fortran-pinned `m_R <= 0.5` rule, calls get_tan, and assembles
// F_sig/F_q. The most useful tests verify the COMPOSITION is correct
// (right operands, right operations, right error transparency), not
// the underlying physics (which test_get_tan already covers).
//
// Coverage:
//   (1) No-IS branch (m_R = 0.4): F_sig = L*deps + N*|deps|_strain.
//       Verified by re-deriving the right-hand side from an
//       independent get_tan call with istrain = 0, then comparing.
//       F_q[0..5] must be zero (H[0..5]==0 in no-IS); F_q[6] must
//       equal (1+v)*(deps_11 + deps_22 + deps_33).
//   (2) IS branch (m_R = 1, textbook): F_sig = M*deps. Same
//       cross-check pattern via independent get_tan(istrain=1) call.
//       F_q[0..5] = H_del * deps. F_q[6] still = (1+v)*tr(deps).
//   (3) Void-evolution scalar identity. The H[6] row is the same
//       formula in both branches and is just trace(deps) scaled by
//       (1+void). Pin this with a non-axisymmetric deps that exercises
//       all three normal slots.
//   (4) Error transparency: r_uc = 0 + m_R > 0.5 must surface
//       error=10 from get_tan with F_sig and F_q left at zero. No
//       partial rate vector escapes.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::RhsResult;
using Kratos::SandHypoCpp::TangentResult;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::get_F_sig_q;
using Kratos::SandHypoCpp::get_tan;
using Kratos::SandHypoCpp::DotProductKind;
using Kratos::SandHypoCpp::dot_vect;
using Kratos::SandHypoCpp::matmul;
using Kratos::SandHypoCpp::kErrorOk;
using Kratos::SandHypoCpp::kStateDim;

namespace {

int g_pass = 0;
int g_fail = 0;

void expect_eq_int(int actual, int expected, const char* label) {
    if (actual == expected) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr, "[FAIL] %s: expected %d, got %d\n",
                     label, expected, actual);
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

Props16 textbook_props() {
    return Props16{ {34.0, 1.0, 3.6e6, 0.43, 0.72, 0.934, 1.2, 0.24, 1.2,
                     1.0,  1.0, 3.3e-4, 0.5, 6.0, 0.0, 0.825} };
}

Parms16 textbook_parms() {
    const auto r = check_parms(textbook_props());
    if (r.error != kErrorOk) {
        std::fprintf(stderr, "[FATAL] textbook_parms: check_parms failed (%s)\n",
                     r.failed ? r.failed : "(unnamed)");
        std::exit(2);
    }
    return r.parms;
}

// Helper: independently compute F_sig the way get_F_sig_q is supposed
// to (call get_tan, then matmul + N*norm_D for no-IS, or matmul for IS).
// This is the cross-check used by cases (1) and (2).
void recompute_F_sig_independent(const double* deps,
                                 const double* sig,
                                 const double* q,
                                 const Parms16& parms,
                                 int istrain,
                                 double F_sig_out[6]) {
    const TangentResult tan_r = get_tan(deps, sig, q, parms, istrain);
    if (istrain == 1) {
        matmul(&tan_r.M[0][0], deps, F_sig_out, 6, 6, 1);
    } else {
        matmul(&tan_r.L[0][0], deps, F_sig_out, 6, 6, 1);
        const double norm_D = std::sqrt(
            dot_vect(DotProductKind::StrainLike, deps, deps, 6));
        for (int i = 0; i < 6; ++i) {
            F_sig_out[i] += tan_r.N[i] * norm_D;
        }
    }
}

// Tolerance for "should match the independent recomputation": both
// sides go through the same matmul/dot/sqrt code paths but with
// independent get_tan calls, so the rebuilt M/L/N/H buffers differ
// in roundoff. A few ulps of relative error is expected.
double match_tolerance(double magnitude) {
    return 1.0e-12 * (std::abs(magnitude) + 1.0);
}

}  // namespace

int main() {
    const Parms16 parms = textbook_parms();
    const double void_ratio = 0.825;
    const double sig_compression[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
    const double deps_small[6]      = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};

    // ====================================================================
    // (1) No-IS branch: m_R = 0.4 < 0.5 forces istrain = 0.
    //     F_sig = L*deps + N*|deps|_strain. F_q[0..5] = 0; F_q[6] = (1+v)*tr(deps).
    // ====================================================================
    {
        Props16 props_low_mR = textbook_props();
        props_low_mR[9] = 0.4;  // m_R = 0.4 -> istrain = 0
        const auto check = check_parms(props_low_mR);
        expect_eq_int(check.error, kErrorOk, "case1 m_R=0.4 accepted");

        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const RhsResult r = get_F_sig_q(deps_small, sig_compression, q,
                                        check.parms);
        expect_eq_int(r.error, kErrorOk, "case1 error == 0");

        // Cross-check F_sig against an independent get_tan + matmul + N*norm_D.
        double F_sig_expected[6];
        recompute_F_sig_independent(deps_small, sig_compression, q,
                                    check.parms, 0, F_sig_expected);
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case1 F_sig[%d] == L*deps + N*|deps| (no-IS)", i);
            expect_near(r.F_sig[i], F_sig_expected[i],
                        match_tolerance(F_sig_expected[i]), label);
        }

        // F_q[0..5] = 0 in no-IS branch (H[0..5] all zero).
        for (int i = 0; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "case1 F_q[%d] == 0 (no-IS)", i);
            expect_near(r.F_q[i], 0.0, 0.0, label);
        }

        // F_q[6] = (1+v) * tr(deps).
        const double tr_deps = deps_small[0] + deps_small[1] + deps_small[2];
        const double F_q6_expected = (1.0 + void_ratio) * tr_deps;
        expect_near(r.F_q[6], F_q6_expected, 1e-15,
                    "case1 F_q[6] == (1+v)*tr(deps)");
    }

    // ====================================================================
    // (2) IS branch: textbook m_R = 1 > 0.5 forces istrain = 1.
    //     F_sig = M*deps. F_q[0..5] = H_del * deps. F_q[6] = (1+v)*tr(deps).
    // ====================================================================
    {
        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const RhsResult r = get_F_sig_q(deps_small, sig_compression, q, parms);
        expect_eq_int(r.error, kErrorOk, "case2 error == 0");

        // Cross-check F_sig against independent M*deps via get_tan.
        double F_sig_expected[6];
        recompute_F_sig_independent(deps_small, sig_compression, q,
                                    parms, 1, F_sig_expected);
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case2 F_sig[%d] == M*deps (IS)", i);
            expect_near(r.F_sig[i], F_sig_expected[i],
                        match_tolerance(F_sig_expected[i]), label);
        }

        // F_q[0..5]: from H_del (zero IS state -> H_del = identity from
        // load<=0 else branch). So F_q[i] = deps[i] for i in 0..5.
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case2 F_q[%d] == deps[%d] (zero-IS H_del = identity)",
                          i, i);
            expect_near(r.F_q[i], deps_small[i], 1e-15, label);
        }

        // F_q[6] = (1+v) * tr(deps), same as case 1.
        const double tr_deps = deps_small[0] + deps_small[1] + deps_small[2];
        const double F_q6_expected = (1.0 + void_ratio) * tr_deps;
        expect_near(r.F_q[6], F_q6_expected, 1e-15,
                    "case2 F_q[6] == (1+v)*tr(deps)");
    }

    // ====================================================================
    // (3) Void-evolution scalar identity. Use a non-axisymmetric deps
    //     that mixes all three normal slots and includes a shear slot.
    //     The void rate is the volumetric strain rate scaled by (1+v),
    //     regardless of branch and regardless of shear-slot content.
    // ====================================================================
    {
        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {-2.0e-3, -1.0e-3, -3.0e-3, 5.0e-4, 0.0, 0.0};
        const double tr_deps = deps[0] + deps[1] + deps[2];
        const double F_q6_expected = (1.0 + void_ratio) * tr_deps;

        // IS branch (textbook m_R = 1).
        {
            const RhsResult r = get_F_sig_q(deps, sig_compression, q, parms);
            expect_eq_int(r.error, kErrorOk, "case3 IS branch error == 0");
            expect_near(r.F_q[6], F_q6_expected, 1e-15,
                        "case3 IS F_q[6] == (1+v)*tr(deps) (deps mixes normals + shear)");
        }

        // No-IS branch (m_R = 0.4).
        {
            Props16 props_low_mR = textbook_props();
            props_low_mR[9] = 0.4;
            const auto check = check_parms(props_low_mR);
            const RhsResult r = get_F_sig_q(deps, sig_compression, q,
                                            check.parms);
            expect_eq_int(r.error, kErrorOk, "case3 no-IS branch error == 0");
            expect_near(r.F_q[6], F_q6_expected, 1e-15,
                        "case3 no-IS F_q[6] == (1+v)*tr(deps)");
        }
    }

    // ====================================================================
    // (4) Error transparency: r_uc = 0 + textbook m_R = 1 must trip
    //     get_tan's IS-branch guard. get_F_sig_q must propagate
    //     error = 10 with F_sig and F_q at zero (no partial rates).
    //
    //     Note: after Codex 2026-05-02 (medium) hoisted "r_uc > 0 if
    //     m_R > 0.5" into check_parms, the combo is rejected at
    //     validation. We bypass check_parms by mutating a validated
    //     Parms16 directly to exercise get_tan's in-kernel
    //     belt-and-suspenders guard. The 4b sanity check still uses
    //     check_parms because m_R=0.4 + r_uc=0 is in the legal regime.
    // ====================================================================
    {
        Parms16 parms_runc_zero = textbook_parms();
        parms_runc_zero[11] = 0.0;          // bypass check_parms cross-check

        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const RhsResult r = get_F_sig_q(deps_small, sig_compression, q,
                                        parms_runc_zero);

        expect_eq_int(r.error, 10,
                      "case4 r_uc=0 + m_R>0.5 -> error=10 propagated from get_tan");

        // F_sig and F_q must remain zero. No partial-rate leak.
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case4 F_sig[%d] == 0 on fatal (no partial rate)", i);
            expect_near(r.F_sig[i], 0.0, 0.0, label);
        }
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case4 F_q[%zu] == 0 on fatal (no partial rate)", i);
            expect_near(r.F_q[i], 0.0, 0.0, label);
        }

        // Sanity cross-check: with m_R = 0.4 (no-IS branch), the
        // r_uc = 0 input is in the legal Fortran regime (the IS
        // branch is gated by m_R > 0.5, so r_uc is unused). check_parms
        // accepts this combo, get_F_sig_q runs to completion.
        Props16 props_low_mR_runc0 = textbook_props();
        props_low_mR_runc0[9]  = 0.4;
        props_low_mR_runc0[11] = 0.0;
        const auto check_lo = check_parms(props_low_mR_runc0);
        expect_eq_int(check_lo.error, kErrorOk,
                      "case4b m_R=0.4 + r_uc=0 still accepted by check_parms");
        const RhsResult r_lo = get_F_sig_q(deps_small, sig_compression, q,
                                           check_lo.parms);
        expect_eq_int(r_lo.error, kErrorOk,
                      "case4b m_R=0.4 + r_uc=0 OK (no-IS path doesn't divide by r_uc)");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

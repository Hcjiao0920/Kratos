// Standalone unit-test exerciser for the Path B kernel `rkf23_update`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// The driver composes every Tier-1 helper (check_RKF, rhs, get_F_sig_q,
// get_tan, norm_res). Tests target the COMPOSER's contract -- nfev
// counting, error-code routing, dtsub writeback ordering, the maxnint
// guard, and the multi-substep accept/reject cycle -- not the physics
// inside the helpers (those are owned by their respective test files).
//
// Coverage (15 cases):
//
// (1) Zero deps + textbook stable state. norm_R = 0 hits the
//     S_hull = 1 branch; strict-< accept fires (0 < err_tol). One
//     substep covers [0, 1]; nfev = 3; y unchanged; r.dtsub = dtime.
//     Also locks the norm_R = 0 branch.
//
// (2) maxnint = 0. The ksubst > maxnint guard fires on the first
//     iteration before any rhs call; nfev = 0; error = 3.
//
// (3) Initial y is tensile (sig_star fails check_RKF(y_k)). Fails
//     at the very first check before any rhs call; nfev = 0; error = 3.
//
// (4) r_uc = 0 with textbook m_R = 1 (IS branch). The first rhs(y_k)
//     hits get_tan's IS-branch r_uc <= 0 guard and returns error = 10.
//     nfev = 1 (counted before the call, per Fortran rhs_h's first
//     line `nfev = nfev + 1`). error = 10 propagates verbatim.
//
// (5) Mid-substep check_RKF(y_2) failure -- locks the gate-2 ordering
//     specifically. Textbook compressive sig + huge tensile deps
//     produces a deterministic kRK_1 push that flips y_2's sig
//     tensile, so check_RKF(y_2) fires and the substep aborts with
//     error = 3 and nfev = 1.
//
// (6) err_tol = 0 with non-trivial deps -> every substep rejects
//     (norm_R < 0 is never true). DT_k shrinks by /4 each iter until
//     below DTmin -> error = 3. Verifies the reject -> DTmin path with
//     a deterministic counter (nfev divisible by 3).
//
// (7) dtsub writeback ordering. Zero deps + dtsub_init = dtime/4 ->
//     two substeps both hit the S_hull = 1 branch. Accept-branch
//     ordering "T_k += DT_k; DT_k = min(4*DT_k, S); dtsub = DT_k*dtime;
//     DT_k = min(1-T_k, DT_k)" means r.dtsub at end of run equals
//     dtime (full step, post-expand pre-clamp), NOT 0 (post-clamp).
//     If the order were reversed, the last accept would write
//     dtsub = 0 because T_k = 1 -> 1-T_k = 0 clamps DT_k to 0.
//
// (8) Oversized dtsub_init silent clamp. dtsub_init = 4*dtime is
//     folded to dtime by the input clamp; the result must equal a
//     fresh dtsub_init = dtime run bit-exactly (no over-integration).
//
// (9) Non-positive dtsub_init -> error = 3, nfev = 0. Both 9a (zero)
//     and 9b (negative) tested.
//
// (10) Reuse chain with an OVERSIZED run1.dtsub. Engineered with small
//      compressive deps + loose err_tol so run1 returns dtsub = 2*dtime
//      (the true overshoot scenario, not a tautology). Run2 with that
//      reused value must match a fresh dtime reference run, proving the
//      clamp prevents over-integration on naive caller reuse.
//
// (11) dtime input contract: dtime <= 0, dtime = NaN, dtime = +Inf, and
//      dtsub_init = NaN must each return error = 3, nfev = 0. Without
//      these guards, dtime < 0 silently produces a successful run with
//      negative r.dtsub via the sign-cancelling DT_k = 1 path.
//
// (12) err_tol input contract: with zero deps (norm_R == 0) plus
//      invalid err_tol (0/NaN/negative/Inf), the strict-< accept fails
//      and the norm_R==0 S_hull=1 branch makes the reject path keep
//      DT_k = 1 forever -- the loop spins until maxnint, doing
//      30000 wasted rhs calls under the default maxnint=10000.
//      The input guard rejects all four invalid err_tol shapes up
//      front with error = 3, nfev = 0.
//
// (13) DTmin input contract: same family of failure as (12) -- non-
//      positive or non-finite DTmin disables the reject-shrunk-DT_k
//      early exit, leading to spin-until-maxnint. Guarded preemptively
//      so the four scalar tuning parameters (dtime, dtsub_init,
//      err_tol, DTmin) form a symmetric input-contract surface.
//
// (14) Void-ratio admissibility through rkf23 (indirect path). y_init
//      with void = 0 trips the new void rule in check_RKF(y_k) on the
//      very first substep; error = 3, nfev = 0, y unchanged.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

using Kratos::SandHypoCpp::FailureClass;
using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::Rkf23Result;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::kErrorOk;
using Kratos::SandHypoCpp::kYDim;
using Kratos::SandHypoCpp::rkf23_update;

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

const char* failure_class_name(FailureClass fc) {
    switch (fc) {
        case FailureClass::Ok:                 return "Ok";
        case FailureClass::RkfReject:          return "RkfReject";
        case FailureClass::ContractViolation:  return "ContractViolation";
        case FailureClass::Fatal:              return "Fatal";
    }
    return "(unknown)";
}

void expect_failure_class(FailureClass actual, FailureClass expected,
                          const char* label) {
    if (actual == expected) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr,
                     "[FAIL] %s: expected failure_class=%s, got %s\n",
                     label, failure_class_name(expected),
                     failure_class_name(actual));
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

void build_y(const double sig[6], const double q[7], double y_out[kYDim]) {
    for (int i = 0; i < 6; ++i) y_out[i]     = sig[i];
    for (int i = 0; i < 7; ++i) y_out[6 + i] = q[i];
}

}  // namespace

int main() {
    const Parms16 parms = textbook_parms();
    const double void_ratio = 0.825;
    const double dtime      = 1.0;     // unit dtime so DT_k == dtsub numerically.
    const double err_tol_default = 1.0e-3;
    const double DTmin_default   = 1.0e-17;
    const int    maxnint_default = 10000;

    const double sig_compression[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
    const double q_zero_is[7]       = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};

    // ====================================================================
    // (1) Zero deps + textbook compressive state -> norm_R = 0, S_hull = 1.
    //     One substep accepts the entire interval; y unchanged; nfev = 3;
    //     r.dtsub = dtime.
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps_zero, parms, dtime);
        expect_eq_int(r.error, kErrorOk, "case1 zero-deps error == 0");
        expect_eq_int(r.nfev,  3,        "case1 zero-deps nfev == 3 (one substep)");
        expect_near(r.dtsub, dtime, 1e-15,
                    "case1 zero-deps r.dtsub == dtime (norm_R=0 -> S_hull=1)");

        // y unchanged: zero ydot means y_hat == y_k; no drift.
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "case1 zero-deps y[%zu] unchanged", i);
            expect_near(r.y[i], y0[i], 0.0, label);
        }
    }

    // ====================================================================
    // (2) maxnint = 0 -> ksubst = 1 > 0 fires on first iter before any
    //     rhs call. nfev = 0; error = 3.
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           /*maxnint=*/0, DTmin_default,
                                           deps, parms, dtime);
        expect_eq_int(r.error, 3, "case2 maxnint=0 error == 3");
        expect_failure_class(r.failure_class, FailureClass::RkfReject,
                             "case2 maxnint=0 failure_class == RkfReject "
                             "(defensive; integrate_step preflights maxnint > 0)");
        expect_eq_int(r.nfev,  0, "case2 maxnint=0 nfev == 0 (no rhs calls)");
        // y on error path stays at y_init.
        for (std::size_t i = 0; i < kYDim; ++i) {
            expect_near(r.y[i], y0[i], 0.0,
                        "case2 maxnint=0 y == y_init on error path");
        }
    }

    // ====================================================================
    // (3) Initial sig is tensile -> check_RKF(y_k) fails on first iter.
    //     nfev = 0; error = 3.
    // ====================================================================
    {
        const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double deps[6]        = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_tensile, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps, parms, dtime);
        expect_eq_int(r.error, 3, "case3 tensile-sig error == 3");
        expect_failure_class(r.failure_class, FailureClass::RkfReject,
                             "case3 tensile-sig failure_class == RkfReject "
                             "(check_RKF(y_k) returns 'pmean' -> stress-domain)");
        expect_eq_int(r.nfev,  0, "case3 tensile-sig nfev == 0 (failed before rhs)");
    }

    // ====================================================================
    // (4) r_uc = 0 + textbook m_R = 1 -> get_tan IS-branch guard fires
    //     on the first rhs(y_k); error = 10 propagates; nfev = 1
    //     (incremented BEFORE the rhs call).
    //
    //     Note: after Codex 2026-05-02 (medium) hoisted "r_uc > 0 if
    //     m_R > 0.5" into check_parms, this combo is now rejected at
    //     validation. We deliberately BYPASS check_parms here by
    //     mutating a validated Parms16 directly -- this exercises
    //     get_tan's in-kernel belt-and-suspenders guard, which is
    //     kept as defense-in-depth for callers that bypass check_parms
    //     (e.g. memcpy, future test seam). Without this test the
    //     belt-and-suspenders layer would have no direct coverage.
    // ====================================================================
    {
        Parms16 parms_runc_zero = textbook_parms();
        parms_runc_zero[11] = 0.0;          // bypass check_parms's cross-check

        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps, parms_runc_zero, dtime);
        expect_eq_int(r.error, 10, "case4 r_uc=0 error == 10 (propagated from get_tan)");
        expect_failure_class(r.failure_class, FailureClass::Fatal,
                             "case4 r_uc=0 failure_class == Fatal");
        expect_eq_int(r.nfev,  1,  "case4 r_uc=0 nfev == 1 (counted before the call)");
    }

    // ====================================================================
    // (5) check_RKF(y_2) failure -- locks the gate-2 path specifically.
    //
    //     Engineering: textbook compressive hydrostatic sig = -100, with
    //     a HUGE tensile strain increment (deps = 1.0 in each normal
    //     slot). For zero-IS textbook params, the tangent stiffness is
    //     stiff enough (L[0][0] ~ 1.3e5, L[0][1] ~ 6e4) that
    //     F_sig[0] = (L[0][0] + 2*L[0][1])*1.0 ~ 2.5e5. With DT_k = 1
    //     and the half-step factor in y_2 = y_k + (DT_k/2)*kRK_1, y_2
    //     gets sig component +1.25e5 -- vastly tensile, so
    //     check_RKF(y_2) MUST fire on the principal-stress branch.
    //
    //     This is a deterministic gate-2 trigger: y_k passes (sig is
    //     compressive at -100), y_2 fails (sig flipped to far-tensile
    //     after k1 push). nfev == 1 LOCKS the gate-2 ordering -- a
    //     regression that skipped or delayed the y_2 check (so y_3 or
    //     y_hat caught it instead) would surface as nfev == 2 or 3.
    //
    //     Note: gates 3 (y_3) and 4 (y_hat) are not individually
    //     exercised -- engineering a state where they fire alone (with
    //     y_2 passing) would require predicting kRK_1 vs kRK_2
    //     interactions that depend sensitively on the constitutive
    //     law. The structural existence of those gates is documented
    //     in the kernel and code-reviewed; this case pins the gate-2
    //     ordering as the most consequential of the three.
    // ====================================================================
    {
        const double deps_huge_tensile[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps_huge_tensile, parms, dtime);
        expect_eq_int(r.error, 3,
                      "case5 gate-2 (y_2 check_RKF) failure error == 3");
        expect_failure_class(r.failure_class, FailureClass::RkfReject,
                             "case5 gate-2 failure_class == RkfReject "
                             "(y_2 trips 'tension' stress-domain label after huge tensile deps)");
        expect_eq_int(r.nfev, 1,
                      "case5 gate-2 ordering locked: nfev == 1 (k1 ran, y_2 failed before k2)");
    }

    // ====================================================================
    // (6) Tiny-but-valid err_tol + non-trivial deps -> every substep
    //     rejects, DT_k shrinks by /4 each iter, eventually below DTmin.
    //
    //     The earlier draft used err_tol = 0 to force every reject.
    //     Once the input guard added in case 12 below rejects err_tol = 0
    //     up front, that approach no longer reaches the substep loop.
    //     Use err_tol = 1e-300 instead: positive and finite (passes the
    //     guard), but small enough vs realistic norm_R (~ 1e-7) that
    //     the strict-< accept never fires and S_hull is essentially
    //     zero, so the reject path computes DT_k_new = DT_k/4 exactly
    //     as before.
    //
    //     Counter check: each iter performs 3 rhs calls. With
    //     dtsub_init = dtime (DT_k = 1) and DTmin = 1e-2 (normalized),
    //     iters needed: DT_k after k iters is 4^-k.
    //       4^-3 = 0.015625 > 0.01,
    //       4^-4 = 0.00390625 < 0.01 -> reject at iter 4.
    //     So 4 iters, each 3 rhs -> nfev = 12.
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);
        const double DTmin_test = 1.0e-2;  // normalized
        const int    maxnint_test = 100;   // generous
        const double err_tol_tiny = 1.0e-300;  // positive, finite, ~ 0 numerically

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_tiny,
                                           maxnint_test, DTmin_test,
                                           deps, parms, dtime);
        expect_eq_int(r.error, 3, "case6 tiny err_tol error == 3 (reject -> DTmin)");
        expect_failure_class(r.failure_class, FailureClass::RkfReject,
                             "case6 DTmin exhaustion failure_class == RkfReject "
                             "(genuine cutback case)");
        expect_eq_int(r.nfev, 12, "case6 tiny err_tol nfev == 12 (4 substeps x 3 rhs)");
    }

    // ====================================================================
    // (7) dtsub writeback ordering. Zero deps + dtsub_init = dtime/4
    //     -> two substeps, each accepts via norm_R = 0 / S_hull = 1.
    //     Final r.dtsub == dtime, NOT 0.
    //
    //     Iter 1: DT_k = 0.25, accept. T_k = 0.25.
    //             DT_k = min(4*0.25, 1) = 1. dtsub = 1*dtime = dtime.
    //             DT_k = min(0.75, 1) = 0.75.
    //     Iter 2: DT_k = 0.75, accept. T_k = 1.
    //             DT_k = min(4*0.75, 1) = 1. dtsub = 1*dtime = dtime.
    //             DT_k = min(0, 1) = 0. Loop exits.
    //
    //     If the writeback order were inverted (clamp before write),
    //     iter 2 would produce dtsub = 0 (because 1-T_k = 0).
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);
        const double dtsub_init = dtime / 4.0;

        const Rkf23Result r = rkf23_update(y0, dtsub_init, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps_zero, parms, dtime);
        expect_eq_int(r.error, kErrorOk, "case7 writeback-order error == 0");
        expect_eq_int(r.nfev,  6,        "case7 writeback-order nfev == 6 (two substeps)");
        expect_near(r.dtsub, dtime, 1e-15,
                    "case7 r.dtsub == dtime (post-expand pre-clamp)");
        // If the write-back were after the clamp, r.dtsub would be 0
        // because the last accept has T_k = 1 -> min(0, 1) = 0.
    }

    // ====================================================================
    // (8) Oversized dtsub_init must NOT over-integrate.
    //     Without an input clamp, dtsub_init = 4*dtime gives DT_k_init = 4,
    //     and the first accepted substep would advance T_k to 4, committing
    //     y at t = 4*dtime instead of t = dtime. Path B clamps the input
    //     silently to dtime; the result must be IDENTICAL to a call with
    //     dtsub_init = dtime.
    //     Codex review 2026-05-01 (high finding).
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r_clean   = rkf23_update(y0, dtime,        err_tol_default,
                                                   maxnint_default, DTmin_default,
                                                   deps_zero, parms, dtime);
        const Rkf23Result r_oversize = rkf23_update(y0, 4.0 * dtime, err_tol_default,
                                                    maxnint_default, DTmin_default,
                                                    deps_zero, parms, dtime);
        expect_eq_int(r_oversize.error, kErrorOk,
                      "case8 oversized dtsub_init error == 0 (clamped, not failed)");
        expect_eq_int(r_oversize.nfev, r_clean.nfev,
                      "case8 nfev matches clean run (clamp prevents extra substeps)");
        expect_near(r_oversize.dtsub, r_clean.dtsub, 1e-15,
                    "case8 r.dtsub matches clean run");
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case8 r.y[%zu] matches clean run (no over-integration)", i);
            expect_near(r_oversize.y[i], r_clean.y[i], 0.0, label);
        }
    }

    // ====================================================================
    // (9) Non-positive dtsub_init -> error = 3, no rhs calls, y unchanged.
    //     Both zero and negative are rejected. The Fortran would enter a
    //     zero-progress substep loop until maxnint fires; Path B refuses up
    //     front and reports the contract violation cleanly.
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        // Sub-case 9a: dtsub_init = 0.
        const Rkf23Result r_zero = rkf23_update(y0, 0.0, err_tol_default,
                                                maxnint_default, DTmin_default,
                                                deps, parms, dtime);
        expect_eq_int(r_zero.error, 3, "case9a dtsub_init=0 -> error=3");
        expect_eq_int(r_zero.nfev,  0, "case9a dtsub_init=0 -> nfev=0");

        // Sub-case 9b: dtsub_init < 0.
        const Rkf23Result r_neg = rkf23_update(y0, -dtime, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, dtime);
        expect_eq_int(r_neg.error, 3, "case9b dtsub_init<0 -> error=3");
        expect_eq_int(r_neg.nfev,  0, "case9b dtsub_init<0 -> nfev=0");

        // y must stay at y_init on the rejection path.
        for (std::size_t i = 0; i < kYDim; ++i) {
            expect_near(r_zero.y[i], y0[i], 0.0,
                        "case9a y unchanged on dtsub<=0 rejection");
            expect_near(r_neg.y[i], y0[i], 0.0,
                        "case9b y unchanged on dtsub<0 rejection");
        }
    }

    // ====================================================================
    // (10) Reuse chain with an OVERSIZED run1.dtsub.
    //
    //      The earlier draft of this case used zero deps + dtsub_init = dtime/4,
    //      which ran into a tautology: with norm_R == 0, S_hull is the
    //      constant 1.0, so 4x expansion gets capped at 1.0 and run1.dtsub
    //      can never exceed dtime. Without the input clamp, run2 with
    //      dtsub_init == dtime would still NOT over-integrate, so the test
    //      could not distinguish clamp-present from clamp-absent.
    //
    //      To genuinely exercise the reuse path, we need run1.dtsub > dtime.
    //      That requires norm_R > 0 (so S_hull is computed via the Hull
    //      formula) AND norm_R << err_tol (so S_hull is large enough that
    //      4*DT_k caps the expansion above 1.0). Setup:
    //        - small compressive deps, smooth constitutive response
    //        - dtsub_init = dtime/2 -> DT_k_init = 0.5, expansion 4x = 2.0
    //        - err_tol = 1.0 (very loose) so norm_R/err_tol is tiny
    //      With these, run1 accepts each substep with maximal expansion and
    //      writes back dtsub = 2 * dtime (the 4x expansion is bound by
    //      4*DT_k = 2, not by S_hull which is much larger).
    //
    //      Then run2 with dtsub_init = run1.dtsub (~= 2*dtime) is the
    //      genuine "reuse" scenario. With the input clamp it folds to
    //      dtime and integrates correctly to t = dtime; without the clamp
    //      it would integrate to t = 2*dtime and produce a different y.
    //
    //      Codex review 2026-05-01 (high finding follow-up).
    // ====================================================================
    {
        const double deps_small[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        const double err_tol_loose  = 1.0;   // generous so norm_R << err_tol
        const double dtsub_init_run1 = dtime / 2.0;
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result run1 = rkf23_update(y0, dtsub_init_run1, err_tol_loose,
                                              maxnint_default, DTmin_default,
                                              deps_small, parms, dtime);
        expect_eq_int(run1.error, kErrorOk, "case10 run#1 error == 0");

        // Setup sanity: run1.dtsub MUST exceed dtime; otherwise the
        // reuse-with-overshoot scenario this case targets is not
        // actually being exercised (the test would degenerate to a
        // tautology like the earlier draft).
        if (run1.dtsub > dtime) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case10 setup: run1.dtsub = %.6g, expected > dtime = %.6g "
                         "(test cannot distinguish clamp-present from clamp-absent)\n",
                         run1.dtsub, dtime);
        }

        // Genuine reuse: feed run1's oversized dtsub back as dtsub_init.
        // With the input clamp this folds to dtime; without it, DT_k_init
        // would be > 1 and the first substep would integrate past
        // T_k = 1, committing a different y than the reference.
        const Rkf23Result run2 = rkf23_update(y0, run1.dtsub, err_tol_loose,
                                              maxnint_default, DTmin_default,
                                              deps_small, parms, dtime);
        const Rkf23Result run_ref = rkf23_update(y0, dtime, err_tol_loose,
                                                 maxnint_default, DTmin_default,
                                                 deps_small, parms, dtime);
        expect_eq_int(run2.error, kErrorOk, "case10 run#2 error == 0");
        expect_eq_int(run2.nfev, run_ref.nfev,
                      "case10 run#2 nfev matches reference (clamp folded oversized to dtime)");
        // y comparison is the primary discriminator: if the clamp were
        // missing, run2 would integrate to t = 2*dtime instead of t = dtime,
        // producing a noticeably different y.
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case10 reuse-chain r.y[%zu] == reference (clamp prevents overshoot)",
                          i);
            expect_near(run2.y[i], run_ref.y[i], 0.0, label);
        }
        expect_near(run2.dtsub, run_ref.dtsub, 1e-15,
                    "case10 run#2 r.dtsub matches reference");
    }

    // ====================================================================
    // (11) dtime input contract: must be finite and strictly positive.
    //
    //      Without a dtime guard:
    //        * dtime < 0 with dtsub_init > 0 silently completes: the
    //          dtsub > dtime clamp folds dtsub_effective to dtime, then
    //          DT_k = dtime/dtime = 1 (sign cancels), the loop runs,
    //          and r.dtsub returns negative.
    //        * dtime == 0 produces DT_k = 0/0 = NaN; the loop wastes
    //          one rhs call before check_RKF(y_2) catches NaN.
    //        * dtime = NaN propagates NaN through the entire driver.
    //
    //      Path B's input guard rejects all three with error = 3, nfev = 0,
    //      y unchanged. The dtsub_init guard separately covers NaN
    //      dtsub_init -- exercised here too via the finiteness check.
    //
    //      Codex review 2026-05-01 (medium follow-up to the dtsub guard).
    // ====================================================================
    {
        const double deps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        // (11a) dtime = 0.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, /*dtime=*/0.0);
            expect_eq_int(r.error, 3, "case11a dtime=0 -> error=3");
            expect_failure_class(r.failure_class, FailureClass::ContractViolation,
                                 "case11a dtime=0 failure_class == ContractViolation "
                                 "(scalar contract guard, NOT cutback)");
            expect_eq_int(r.nfev,  0, "case11a dtime=0 -> nfev=0 (guard before any work)");
        }

        // (11b) dtime < 0.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, /*dtime=*/-1.0);
            expect_eq_int(r.error, 3, "case11b dtime<0 -> error=3");
            expect_eq_int(r.nfev,  0, "case11b dtime<0 -> nfev=0");
        }

        // (11c) dtime = NaN.
        {
            const double dtime_nan = std::numeric_limits<double>::quiet_NaN();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, dtime_nan);
            expect_eq_int(r.error, 3, "case11c dtime=NaN -> error=3");
            expect_eq_int(r.nfev,  0, "case11c dtime=NaN -> nfev=0");
        }

        // (11d) dtime = +Inf -- finiteness check catches.
        {
            const double dtime_inf = std::numeric_limits<double>::infinity();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, dtime_inf);
            expect_eq_int(r.error, 3, "case11d dtime=+Inf -> error=3");
            expect_eq_int(r.nfev,  0, "case11d dtime=+Inf -> nfev=0");
        }

        // (11e) Bonus: dtsub_init = NaN (additional finiteness coverage on
        //       the dtsub_init guard, paired with the dtime checks here).
        {
            const double dtsub_nan = std::numeric_limits<double>::quiet_NaN();
            const Rkf23Result r = rkf23_update(y0, dtsub_nan, err_tol_default,
                                               maxnint_default, DTmin_default,
                                               deps, parms, dtime);
            expect_eq_int(r.error, 3, "case11e dtsub_init=NaN -> error=3");
            expect_eq_int(r.nfev,  0, "case11e dtsub_init=NaN -> nfev=0");
        }
    }

    // ====================================================================
    // (12) err_tol input contract: must be finite and strictly positive.
    //
    //      Without an err_tol guard, the failure mode flagged by Codex
    //      is uniquely nasty: with zero deps (norm_R == 0) and
    //      err_tol = 0 / NaN / negative, the strict-< accept check
    //      fails (no value < {0,NaN,negative}), the S_hull = 1 (constant)
    //      branch fires because norm_R == 0, and the reject path's
    //      `max(DT_k/4, 1)` keeps DT_k = 1 forever. T_k never advances
    //      and ksubst increments each iteration, doing 3 rhs calls per
    //      iter, until ksubst > maxnint -- with the default
    //      maxnint = 10000, that is 30000 wasted rhs calls before
    //      error = 3 finally fires.
    //
    //      Path B's input guard rejects all four invalid err_tol
    //      shapes up front: error = 3, nfev = 0. The zero-deps state
    //      is used here specifically to exercise the norm_R == 0 spin
    //      path that the existing case 6 (non-zero deps) does NOT
    //      reach.
    //
    //      Codex review 2026-05-01 (medium follow-up).
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        // (12a) err_tol = 0 with zero deps -- the specific
        //       norm_R==0-spin failure path Codex flagged.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, /*err_tol=*/0.0,
                                               maxnint_default, DTmin_default,
                                               deps_zero, parms, dtime);
            expect_eq_int(r.error, 3, "case12a err_tol=0 + zero deps -> error=3");
            expect_eq_int(r.nfev,  0, "case12a err_tol=0 + zero deps -> nfev=0 (no spin)");
        }

        // (12b) err_tol < 0.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, /*err_tol=*/-1.0,
                                               maxnint_default, DTmin_default,
                                               deps_zero, parms, dtime);
            expect_eq_int(r.error, 3, "case12b err_tol<0 -> error=3");
            expect_eq_int(r.nfev,  0, "case12b err_tol<0 -> nfev=0");
        }

        // (12c) err_tol = NaN.
        {
            const double err_tol_nan = std::numeric_limits<double>::quiet_NaN();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_nan,
                                               maxnint_default, DTmin_default,
                                               deps_zero, parms, dtime);
            expect_eq_int(r.error, 3, "case12c err_tol=NaN -> error=3");
            expect_eq_int(r.nfev,  0, "case12c err_tol=NaN -> nfev=0");
        }

        // (12d) err_tol = +Inf.
        {
            const double err_tol_inf = std::numeric_limits<double>::infinity();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_inf,
                                               maxnint_default, DTmin_default,
                                               deps_zero, parms, dtime);
            expect_eq_int(r.error, 3, "case12d err_tol=+Inf -> error=3");
            expect_eq_int(r.nfev,  0, "case12d err_tol=+Inf -> nfev=0");
        }
    }

    // ====================================================================
    // (13) DTmin input contract: must be finite and strictly positive.
    //
    //      Same family of failure as case 12 (loop spins until maxnint
    //      because the reject path's `if (DT_k < DTmin)` early-exit
    //      can never fire with non-positive or non-finite DTmin).
    //      Caught preemptively here so the input-contract surface for
    //      the four scalar tuning parameters (dtime, dtsub_init,
    //      err_tol, DTmin) is fully symmetric.
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        // (13a) DTmin = 0.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, /*DTmin=*/0.0,
                                               deps, parms, dtime);
            expect_eq_int(r.error, 3, "case13a DTmin=0 -> error=3");
            expect_eq_int(r.nfev,  0, "case13a DTmin=0 -> nfev=0");
        }

        // (13b) DTmin < 0.
        {
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, /*DTmin=*/-1.0,
                                               deps, parms, dtime);
            expect_eq_int(r.error, 3, "case13b DTmin<0 -> error=3");
            expect_eq_int(r.nfev,  0, "case13b DTmin<0 -> nfev=0");
        }

        // (13c) DTmin = NaN.
        {
            const double DTmin_nan = std::numeric_limits<double>::quiet_NaN();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_nan,
                                               deps, parms, dtime);
            expect_eq_int(r.error, 3, "case13c DTmin=NaN -> error=3");
            expect_eq_int(r.nfev,  0, "case13c DTmin=NaN -> nfev=0");
        }

        // (13d) DTmin = +Inf -- finiteness check catches.
        {
            const double DTmin_inf = std::numeric_limits<double>::infinity();
            const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                               maxnint_default, DTmin_inf,
                                               deps, parms, dtime);
            expect_eq_int(r.error, 3, "case13d DTmin=+Inf -> error=3");
            expect_eq_int(r.nfev,  0, "case13d DTmin=+Inf -> nfev=0");
        }
    }

    // ====================================================================
    // (14) Void-ratio admissibility on the substep gate. y_init with
    //      void = 0 must trip check_RKF(y_k) at the very first iteration;
    //      no rhs work is performed, error = 3, y is left at y_init.
    //
    //      Without check_RKF's new void rule, this state would flow into
    //      get_tan(y_k), where `pow(ec / 0, beta)` = +Inf, propagate
    //      through the tangent matrix, and either fail later in the
    //      substep loop via the NaN scan OR (less likely) escape with a
    //      finite-looking but corrupted y_hat.
    //
    //      Codex review 2026-05-01 (medium).
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        // Build y0 manually with void = 0 (build_y helper sets void=0.825).
        double y0[kYDim] = {0.0};
        for (int i = 0; i < 6; ++i) y0[i] = sig_compression[i];
        y0[12] = 0.0;  // non-positive void

        const Rkf23Result r = rkf23_update(y0, dtime, err_tol_default,
                                           maxnint_default, DTmin_default,
                                           deps, parms, dtime);
        expect_eq_int(r.error, 3, "case14 void=0 in y_init -> error=3");
        expect_failure_class(r.failure_class, FailureClass::ContractViolation,
                             "case14 void=0 failure_class == ContractViolation "
                             "(check_RKF returns 'void', corruption-class label)");
        expect_eq_int(r.nfev,  0, "case14 void=0 -> nfev=0 (caught at check_RKF(y_k))");

        // y unchanged on error path.
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case14 r.y[%zu] == y_init[%zu] (void-fail unchanged)", i, i);
            expect_near(r.y[i], y0[i], 0.0, label);
        }
    }

    // ====================================================================
    // (15) Poisoned-residual fail-closed gate (Codex 2026-05-02 P2
    //      finding -- structural lock).
    //
    //      `norm_res` returns `nan_detected = true` and the sentinel
    //      `norm_R = 1.0e20` whenever any intermediate norm
    //      (norm_sig, norm_q, void_hat) is non-finite or exceeds 1e30,
    //      OR the post-check on the returned norm_R itself is
    //      non-finite. Without an explicit gate at the rkf23_update
    //      site, a caller passing `err_tol > 1.0e20` would let the
    //      sentinel value satisfy the strict-< accept comparison and
    //      commit a step whose error estimate was explicitly poisoned.
    //      kernel.cpp:1124-1140 inserts an unconditional reject when
    //      `nres.nan_detected` is true, BEFORE the tolerance compare.
    //
    //      Behavioral note: producing `nan_detected == true` while
    //      `check_RKF(y_hat)` simultaneously passes is rare under the
    //      normal substep evolution (rhs guards / check_RKF's
    //      corruption-class precedence catch most paths first). The
    //      lock at this layer is therefore primarily structural; the
    //      direct mechanics of `nan_detected` are exhaustively
    //      covered by test_rkf_aux N4-N9. Here we only verify the
    //      POSITIVE control: a large `err_tol = 1e25` (the threshold
    //      Codex flagged) does NOT break normal operation in the
    //      absence of poisoning.
    //
    //      If a future change exposes a path producing `nan_detected`
    //      naturally during substepping, add a behavioral negative
    //      test alongside this positive control.
    // ====================================================================
    {
        const double huge_tol = 1.0e25;
        const double deps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double y0[kYDim]; build_y(sig_compression, q_zero_is, y0);

        const Rkf23Result r = rkf23_update(y0, dtime, huge_tol,
                                           maxnint_default, DTmin_default,
                                           deps, parms, dtime);
        expect_eq_int(r.error, kErrorOk,
                      "case15 huge err_tol=1e25 + zero deps -> success "
                      "(nan_detected gate must NOT fire spuriously)");
        // y unchanged.
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[96];
            std::snprintf(label, sizeof(label),
                          "case15 r.y[%zu] == y_init[%zu] (huge err_tol normal-success)",
                          i, i);
            expect_near(r.y[i], y0[i], 0.0, label);
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

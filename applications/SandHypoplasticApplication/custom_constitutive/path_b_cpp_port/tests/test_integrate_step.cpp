// Standalone unit-test exerciser for the Path B kernel `integrate_step`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// `integrate_step` is the top-level step composer: check_RKF gate +
// (RKF main path | calc_elasti fallback) + final get_tan tangent +
// dtsub write-back clamp. Tests target the COMPOSER's control flow --
// branch routing, error-code propagation, dtsub clamping -- not the
// physics already covered by the lower-level test_*.cpp files.
//
// Coverage (15 cases):
//   (1) Valid compressive + zero deps -> zero-deps short circuit
//       (Fortran .for:254 testing=2; rkf23 NOT called); state
//       unchanged; nfev = 0; D = m_R * L (textbook m_R = 1, IS branch);
//       used_elastic = false. The final tangent is still computed via
//       the shared get_tan call after the short circuit.
//   (2) Initial tensile sig -> elastic fallback fires; nfev = 0;
//       q unchanged; sig = sig_n + d_sig from calc_elasti;
//       used_elastic = true; error = 0.
//   (3) RKF error = 3 (huge tensile deps trips a substep gate after
//       k1 succeeds) -> sig and q reverted to inputs; D zero.
//   (4) RKF error = 10 (r_uc = 0 + textbook m_R = 1) -> propagates
//       through; sig and q reverted; D zero.
//   (5) dtsub write-back clamp: setup that makes rkf23 return
//       dtsub > dtime -> integrate_step's output dtsub_next == dtime.
//   (6) Final-tangent branch selection follows parms[9] <= 0.5:
//       parms[9] = 0.4 -> istrain = 0 -> D = L
//       parms[9] = 0.6 -> istrain = 1 -> D = m_R * L
//       Compare integrate_step.D against an independent get_tan call
//       to lock the m_R scaling.
//   (7) Initial-gate fail-closed: only "pmean"/"tension" failure modes
//       trigger elastic fallback. Numerical / void / overflow failures
//       must surface as error = 3 with state reverted to inputs and
//       D = 0; never as a successful elastic step (which would let
//       corrupted state persist invisibly across calls).
//       Codex review 2026-05-01 high finding follow-up.
//   (8) Scalar tuning parameter contract: dtime, dtsub_in, err_tol,
//       and DTmin must each be finite and strictly positive. The
//       elastic-fallback branch bypasses rkf23_update entirely, so
//       without this guard an invalid scalar (e.g. dtsub_in = NaN,
//       dtime < 0) would leak through the elastic path and corrupt
//       r.dtsub_next via the bottom clamp (NaN survives the
//       comparisons; negative dtime pins r.dtsub_next negative).
//       Reject at entry with error = 3 and r.dtsub_next = 0.
//       Codex review 2026-05-01 (medium follow-up to case 7's high).
//   (9) Strain increment contract: every entry of deps_np1 must be
//       finite with |v| <= 1e30. The elastic-fallback branch passes
//       deps_np1 directly to calc_elasti -- which has no error channel
//       -- so a NaN strain increment with tensile initial sig would
//       slip past the inittension routing as a "successful" elastic
//       step with NaN-corrupted r.sig. The RKF main path catches NaN
//       deps via the chain reaction through rhs / get_tan / check_RKF,
//       but the short elastic path bypasses those gates. Reject any
//       non-finite or overflow deps at entry with error = 3.
//       Codex review 2026-05-02 (high finding follow-up).
//  (10) Final-tangent state regression: the tangent published in r.D
//       after a successful RKF step must be computed at the
//       BEGIN-OF-STEP state (sig_n, q_n) -- matching Fortran
//       perturbate_h .for:1466-1469 -- NOT at the post-RKF state
//       (r.sig, r.q). Case (6) used zero deps so y_np1 == y_n and
//       could not distinguish the two states. This case uses non-
//       trivial compressive deps so the state advances meaningfully,
//       then asserts (i) D == m_R * get_tan(deps, sig_n, q_n).L and
//       (ii) get_tan(deps, sig_n, q_n).L differs from
//       get_tan(deps, r.sig, r.q).L by a meaningful margin. Without
//       both assertions, a regression that flips the state argument
//       back to (r.sig, r.q) could pass test (6) silently.
//       Codex review 2026-05-02 (high finding -- direct regression lock).
//  (15) FailureClass disambiguation: comprehensive 4-way demonstration
//       that the new `r.failure_class` enum field carries the routing
//       signal that `r.error` alone cannot. The bridge needs to
//       distinguish retryable RKF rejects (cutback outer dtime) from
//       contract violations (escalate to KRATOS_ERROR) without doing
//       its own intra-kernel preflight (which it cannot, per the
//       Codex 2026-05-02 P1 high finding -- check_RKF derives
//       min_principal internally and the bridge cannot easily
//       replicate that work). This case fires four scenarios in one
//       block, each with the same `error` field value where
//       applicable, and asserts `failure_class` distinguishes them.
//  (14) Derived-stress-overflow regression (Codex 2026-05-02 P1
//       high finding). With sig = [0,0,0,1e30,1e30,1e30] every
//       component lies on the strict `> 1e30` per-component
//       boundary (passes integrate_step's input pre-flight) but
//       principal_stresses_3 derives eigenvalues ~2e30, which
//       trips both the "tension" threshold and the umatisnan_h
//       overflow sentinel inside check_RKF. Pre-fix, "tension" was
//       set first and preserved by first-failure-wins; integrate_
//       step routed to elastic fallback and emitted error=0 with
//       overflow-finite r.sig. Post-fix, the new corruption-class
//       precedence ("nan" overrides "tension") makes integrate_step
//       fail closed with error=3 and preserve r.sig.
//
//       This is the integrate_step-level lock; the check_RKF-level
//       lock for the same precedence is in test_check_rkf R13.
//  (13) error=10 rollback contract through validated-but-bypass parms.
//       integrate_step has a state-revert + D=0 contract for the
//       rkf23-side error=10 path (rkf23_update returns error=10,
//       integrate_step must propagate with sig=sig_n, q=q_n, D=0).
//       Codex 2026-05-02 P2 finding: this contract had no DIRECT
//       coverage after case 4 was repurposed to test check_parms
//       rejection. fb_temp1<0 is a parms-only get_tan failure mode
//       that check_parms does NOT validate (purely the combination
//       phi/ei0/ed0/ec0/alpha that yields a negative coefficient),
//       so the rollback path IS reachable through validated parms;
//       leaving it untested is a real coverage gap.
//
//       This case takes the simpler approach used by test_rkf23 /
//       test_F_sig_q / test_rhs: bypass check_parms's cross-check
//       by directly mutating parms[11]=0 on a validated Parms16,
//       which reaches get_tan's r_uc belt-and-suspenders guard,
//       which emits error=10 from rkf23_update, which integrate_step
//       must propagate per its rollback contract.
//  (12) First-call dtsub normalization (Codex 2026-05-02 high finding).
//       Fortran umat .for:243-244 normalizes `dtsub <= 0` and
//       `dtsub > dtime` to `dtime` before integration; Kratos / Abaqus
//       zero-initialize persisted statev, so the very first call to
//       `integrate_step` arrives with dtsub_in == 0. Earlier code
//       rejected this as a contract violation -- which would deadlock
//       the bridge (caller cuts back, reads zero from statev, retries,
//       reads zero again, no progress). The kernel now mirrors Fortran:
//         - finite dtsub_in <= 0 OR > dtime -> normalize to dtime
//         - non-finite (NaN/Inf) dtsub_in -> still hard reject (case 8a/f)
//       Coverage:
//         (12a) compressive sig + dtsub_in = 0 -> RKF main path success;
//               equivalent to dtsub_in = dtime.
//         (12b) tensile sig + dtsub_in = 0 -> elastic fallback success;
//               r.dtsub_next == dtime.
//         (12c) dtsub_in = -1.0 -> normalized; equivalent to dtsub_in = dtime.
//         (12d) dtsub_in = 2 * dtime -> normalized; equivalent to dtsub_in = dtime.
//         (12e) NaN dtsub_in still hard-rejects (delta vs case 8a: this
//               case explicitly asserts the boundary between "normalize"
//               and "reject" semantics in one suite, so a regression
//               that flips either side fails locally).
//  (11) Strain norm overflow guard: even when every deps_np1[i] passes
//       the per-component |v| <= 1e30 check, pairs of large finite
//       components can push ||deps||_strain above 1e30. Fortran
//       umat (.for:201-209) computes norm_D = sqrt(dot_vect_h(2,
//       deps, deps, 6)) FIRST and runs umatisnan_h on the scalar
//       norm before any branch decision. Without a matching norm
//       guard, deps = (8e29, 8e29, 0, 0, 0, 0) + tensile sig would
//       route to the elastic fallback, calc_elasti would produce
//       d_sig ~ 1e32 magnitude (finite but overflow-scale), and
//       integrate_step would return error == 0 with corrupted r.sig.
//       Codex review 2026-05-02 (high finding -- norm-level follow-up).

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

using Kratos::SandHypoCpp::FailureClass;
using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::StepResult;
using Kratos::SandHypoCpp::TangentResult;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::get_tan;
using Kratos::SandHypoCpp::integrate_step;
using Kratos::SandHypoCpp::kErrorOk;
using Kratos::SandHypoCpp::kStateDim;

namespace {

int g_pass = 0;
int g_fail = 0;

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

void expect_eq_int(int actual, int expected, const char* label) {
    if (actual == expected) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr, "[FAIL] %s: expected %d, got %d\n",
                     label, expected, actual);
    }
}

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

}  // namespace

int main() {
    const Parms16 parms = textbook_parms();
    const double void_ratio = 0.825;
    const double dtime      = 1.0;

    const double sig_compression[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
    const double q_zero_is[7]       = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};

    // ====================================================================
    // (1) Valid compressive + zero deps -> zero-deps short circuit.
    //     Mirrors Fortran .for:254 (testing=2): rkf23_update is NOT
    //     called when ||deps||_strain == 0. Final state == initial
    //     state, nfev == 0, used_elastic == false. The shared final-
    //     tangent get_tan call still fires (Fortran's perturbate_h at
    //     .for:328 runs regardless of testing). With textbook m_R = 1,
    //     D = 1 * L = L (post-fs scale).
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const StepResult r = integrate_step(sig_compression, q_zero_is,
                                            deps_zero, parms, dtime, dtime);
        expect_eq_int(r.error, kErrorOk, "case1 valid compressive error == 0");
        expect_failure_class(r.failure_class, FailureClass::Ok,
                             "case1 failure_class == Ok (success)");
        expect_eq_int(r.nfev, 0, "case1 nfev == 0 (zero-deps short circuit)");
        expect_eq_bool(r.used_elastic, false, "case1 used_elastic == false");

        for (int i = 0; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "case1 r.sig[%d] unchanged", i);
            expect_near(r.sig[i], sig_compression[i], 0.0, label);
        }
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "case1 r.q[%zu] unchanged", i);
            expect_near(r.q[i], q_zero_is[i], 0.0, label);
        }

        // dtsub clamp: rkf23 returns dtime, clamp keeps it.
        expect_near(r.dtsub_next, dtime, 1e-15, "case1 r.dtsub_next == dtime");

        // D should be non-trivial (non-zero diagonal, symmetric).
        bool diag_nonzero = (r.D[0][0] > 0.0 && r.D[3][3] > 0.0);
        if (diag_nonzero) ++g_pass;
        else { ++g_fail; std::fprintf(stderr,
            "[FAIL] case1 D diagonal: D[0][0]=%g D[3][3]=%g\n", r.D[0][0], r.D[3][3]); }
    }

    // ====================================================================
    // (2) Initial tensile sig -> elastic fallback fires.
    //     check_RKF on (sig=10) state: pmean = -(3*9)/3 = -9 < 0.25 -> fail.
    //     calc_elasti(deps, 100, 0.48) computes elastic D and d_sig.
    //     With small compressive deps, d_sig has small compressive components.
    // ====================================================================
    {
        const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double deps[6]        = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                            deps, parms, dtime, dtime);
        expect_eq_int(r.error, kErrorOk, "case2 elastic fallback error == 0");
        expect_failure_class(r.failure_class, FailureClass::Ok,
                             "case2 failure_class == Ok (elastic fallback IS success)");
        expect_eq_int(r.nfev, 0, "case2 nfev == 0 (no rhs calls in elastic path)");
        expect_eq_bool(r.used_elastic, true, "case2 used_elastic == true");

        // q unchanged.
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case2 r.q[%zu] unchanged on elastic path", i);
            expect_near(r.q[i], q_zero_is[i], 0.0, label);
        }

        // sig = sig_n + d_sig. d_sig[i] for i in 0..2 is compressive (negative)
        // because deps was compressive and elastic D is positive-definite.
        for (int i = 0; i < 3; ++i) {
            if (r.sig[i] < sig_tensile[i]) {
                ++g_pass;
            } else {
                ++g_fail;
                std::fprintf(stderr,
                             "[FAIL] case2 r.sig[%d] not more compressive: sig_n=%g, sig_new=%g\n",
                             i, sig_tensile[i], r.sig[i]);
            }
        }
        // Shear slots untouched.
        for (int i = 3; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "case2 r.sig[%d] == 0 (no shear deps)", i);
            expect_near(r.sig[i], 0.0, 0.0, label);
        }

        // D should be the elastic stiffness (D[0][0] > 0 etc.).
        if (r.D[0][0] > 0.0 && r.D[3][3] > 0.0) ++g_pass;
        else { ++g_fail; std::fprintf(stderr,
            "[FAIL] case2 elastic D positive-diagonal: D[0][0]=%g\n", r.D[0][0]); }
    }

    // ====================================================================
    // (3) RKF error = 3: huge tensile deps trips a substep gate.
    //     This is case 5 from test_rkf23 ported up: deterministically
    //     fails check_RKF(y_2) after k1 was computed. integrate_step
    //     must propagate error = 3 with sig and q reverted to inputs
    //     and D = 0.
    // ====================================================================
    {
        const double deps_huge_tensile[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
        const StepResult r = integrate_step(sig_compression, q_zero_is,
                                            deps_huge_tensile, parms,
                                            dtime, dtime);
        expect_eq_int(r.error, 3, "case3 RKF reject error == 3");
        expect_failure_class(r.failure_class, FailureClass::RkfReject,
                             "case3 failure_class == RkfReject (retryable cutback)");
        expect_eq_bool(r.used_elastic, false,
                       "case3 used_elastic == false (we DID enter RKF before reject)");

        // sig and q reverted to inputs.
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case3 r.sig[%d] reverted to sig_n", i);
            expect_near(r.sig[i], sig_compression[i], 0.0, label);
        }
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case3 r.q[%zu] reverted to q_n", i);
            expect_near(r.q[i], q_zero_is[i], 0.0, label);
        }

        // D zero on error path.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case3 r.D[%d][%d] == 0 on error", i, j);
                expect_near(r.D[i][j], 0.0, 0.0, label);
            }
        }
    }

    // ====================================================================
    // (4) Cross-parameter constraint hoisted into `check_parms`:
    //     r_uc = 0 + m_R > 0.5 is now rejected at validation, BEFORE
    //     `integrate_step` can run.
    //
    //     History: this case used to exercise the rkf23-side error=10
    //     rollback in integrate_step by feeding r_uc=0 + textbook m_R=1
    //     through check_parms (Fortran-style accept) and watching
    //     get_tan emit error=10 mid-substep. Codex 2026-05-02 (medium)
    //     flagged the resulting path-dependence: tensile sig + same
    //     parms returned `error = 0` via the elastic fallback (no r_uc
    //     dependency), then a later compressive step would emit fatal
    //     "out of nowhere". The fix hoisted "r_uc > 0 if m_R > 0.5"
    //     into check_parms, eliminating path-dependence and rejecting
    //     the bad material at validation time.
    //
    //     Consequences for integrate_step coverage:
    //       - The rkf23-side error=10 rollback is now structurally
    //         unreachable through validated parms. get_tan's r_uc
    //         guard remains as belt-and-suspenders.
    //       - The final-get_tan rollback after a successful RKF was
    //         already documented as unreachable (case 4's old
    //         comment); that remains true.
    //       - The integrate_step error=10 revert + D=0 contract is
    //         symmetric to the error=3 contract tested by case 3;
    //         test_rkf23 still covers error=10 emission directly.
    //
    //     This case now LOCKS the new rejection at check_parms (so a
    //     future regression that re-accepts the combo is caught) AND
    //     documents that integrate_step is no longer the place where
    //     this verdict is rendered.
    // ====================================================================
    {
        Props16 props_runc_zero = textbook_props();
        props_runc_zero[11] = 0.0;          // r_uc = 0
        // textbook m_R = 1.0 (parms[9]), well above the 0.5 IS threshold
        const auto check = check_parms(props_runc_zero);
        expect_eq_int(check.error, 10,
                      "case4 r_uc=0 + m_R>0.5 now REJECTED by check_parms");

        const bool name_ok = (check.failed != nullptr)
                          && (std::strcmp(check.failed, "r_uc") == 0);
        if (name_ok) ++g_pass;
        else { ++g_fail; std::fprintf(stderr,
            "[FAIL] case4 failed=%s, expected r_uc\n",
            check.failed ? check.failed : "(null)"); }

        expect_near(check.bad_value, 0.0, 0.0,
                    "case4 bad_value == 0.0 (the rejected r_uc)");

        // Negative control: same r_uc=0 with m_R = 0.5 (boundary, IS off
        // per strict `>` rule) must still be accepted -- proves the
        // cross-check is gated by the IS-active condition, not by r_uc
        // alone (which would over-reject Fortran-legal materials).
        Props16 props_is_off = textbook_props();
        props_is_off[11] = 0.0;             // r_uc = 0
        props_is_off[9]  = 0.5;             // m_R at IS threshold
        const auto check_off = check_parms(props_is_off);
        expect_eq_int(check_off.error, kErrorOk,
                      "case4 r_uc=0 + m_R=0.5 (IS off) still accepted");
    }

    // ====================================================================
    // (5) dtsub write-back clamp.
    //     Setup borrowed from test_rkf23 case 10: small compressive
    //     deps + dtsub_init = dtime/2 + loose err_tol -> rkf23 returns
    //     r.dtsub = 2 * dtime. integrate_step must clamp to dtime.
    //
    //     Without the clamp, integrate_step would publish r.dtsub_next
    //     = 2 * dtime, which (when persisted by the caller and reused
    //     as the next dtsub_in) would trigger rkf23_update's input
    //     clamp on the next call -- functionally OK but semantically
    //     ugly, and a wrapper-side clamp matches the Fortran
    //     write-back style.
    // ====================================================================
    {
        const double deps_small[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        const double err_tol_loose = 1.0;
        const double dtsub_in_test = dtime / 2.0;
        const StepResult r = integrate_step(sig_compression, q_zero_is,
                                            deps_small, parms,
                                            dtime, dtsub_in_test,
                                            err_tol_loose);
        expect_eq_int(r.error, kErrorOk, "case5 dtsub clamp error == 0");
        expect_near(r.dtsub_next, dtime, 1e-15,
                    "case5 r.dtsub_next clamped to dtime");
    }

    // ====================================================================
    // (6) Final-tangent branch selection: D = L for no-IS, D = m_R * L
    //     for IS. Use zero deps so final state == initial state, which
    //     lets us compute a reference L from a direct get_tan call on
    //     the same inputs.
    //
    //     parms[9] = 0.4 (m_R <= 0.5) -> integrate_step's istrain = 0,
    //         D[i][j] should equal L[i][j] (no scaling).
    //     parms[9] = 0.6 (m_R > 0.5)  -> istrain = 1,
    //         D[i][j] should equal m_R * L[i][j].
    // ====================================================================
    {
        const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        // (6a) no-IS branch: m_R = 0.4.
        {
            Props16 props_lo_mR = textbook_props();
            props_lo_mR[9] = 0.4;
            const auto check = check_parms(props_lo_mR);
            expect_eq_int(check.error, kErrorOk, "case6a m_R=0.4 accepted");

            const StepResult r = integrate_step(sig_compression, q_zero_is,
                                                deps_zero, check.parms,
                                                dtime, dtime);
            expect_eq_int(r.error, kErrorOk, "case6a no-IS RKF success");

            // Independent get_tan with istrain = 0 to get reference L.
            const TangentResult ref = get_tan(deps_zero, sig_compression,
                                              q_zero_is, check.parms,
                                              /*istrain=*/0);
            expect_eq_int(ref.error, kErrorOk, "case6a reference get_tan ok");

            // D should equal L (post-fs scaled).
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double expected = ref.L[i][j];
                    const double tol = 1e-9 * (std::abs(expected) + 1.0);
                    char label[80];
                    std::snprintf(label, sizeof(label),
                                  "case6a D[%d][%d] == L[%d][%d] (no-IS, scale=1)",
                                  i, j, i, j);
                    expect_near(r.D[i][j], expected, tol, label);
                }
            }
        }

        // (6b) IS branch: m_R = 0.6.
        {
            Props16 props_hi_mR = textbook_props();
            props_hi_mR[9] = 0.6;
            const auto check = check_parms(props_hi_mR);
            expect_eq_int(check.error, kErrorOk, "case6b m_R=0.6 accepted");

            const StepResult r = integrate_step(sig_compression, q_zero_is,
                                                deps_zero, check.parms,
                                                dtime, dtime);
            expect_eq_int(r.error, kErrorOk, "case6b IS RKF success");

            const TangentResult ref = get_tan(deps_zero, sig_compression,
                                              q_zero_is, check.parms,
                                              /*istrain=*/1);
            expect_eq_int(ref.error, kErrorOk, "case6b reference get_tan ok");

            const double m_R = check.parms[9];
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double expected = m_R * ref.L[i][j];
                    const double tol = 1e-9 * (std::abs(expected) + 1.0);
                    char label[80];
                    std::snprintf(label, sizeof(label),
                                  "case6b D[%d][%d] == m_R * L[%d][%d] (IS, scale=0.6)",
                                  i, j, i, j);
                    expect_near(r.D[i][j], expected, tol, label);
                }
            }
        }

        // (6c) Killer non-degenerate test (Codex 2026-05-02 P2 follow-up):
        //      cases 6a/6b use deps_zero + q_zero_is, so norm_del = 0,
        //      rho = 0, AA = 0, and M happens to equal m_R * L (no IS
        //      contribution). A regression that erroneously returned
        //      `tan_final.M` instead of `m_R * tan_final.L` would slip
        //      through 6b silently. To distinguish, drive the IS
        //      branch into a regime where M and m_R * L diverge:
        //        - nonzero q (norm_del > 0 -> rho > 0)
        //        - m_R != m_T (otherwise mm_temp1 = m_R*fs and AA
        //          collapses again)
        //      Setup mirrors test_get_tan case 6: m_R = 5, m_T = 2,
        //      del = (1e-4, 0, 0, 0, 0, 0). Use deps aligned with del
        //      so load > 0 and AA[*][0] = mm_temp2*Leta + mm_temp3*N
        //      contributes (verified non-zero in test_get_tan case 9).
        //
        //      Killer assertions:
        //        (i)  D == m_R * tan_at_yn.L (the perturbate_h convention)
        //        (ii) D != tan_at_yn.M       (the would-be-bug return)
        //        (iii) tan_at_yn.M differs MEANINGFULLY from m_R * L
        //              somewhere on the j == 0 column -- proves the
        //              two-state distinction is real, not FP noise.
        {
            Props16 props_split_mR = textbook_props();
            props_split_mR[9]  = 5.0;   // m_R
            props_split_mR[10] = 2.0;   // m_T
            const auto check6c = check_parms(props_split_mR);
            expect_eq_int(check6c.error, kErrorOk,
                          "case6c m_R=5 m_T=2 accepted");

            const double del0   = 1.0e-4;
            const double q_is[7] = {del0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                    /*void=*/0.825};
            const double deps_aligned[6] = {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};

            const StepResult r = integrate_step(sig_compression, q_is,
                                                deps_aligned, check6c.parms,
                                                dtime, dtime);
            expect_eq_int(r.error, kErrorOk, "case6c RKF success");
            expect_eq_bool(r.used_elastic, false,
                           "case6c used_elastic == false");

            // Begin-of-step tangent reference at y_n (= sig_compression, q_is).
            const TangentResult tan_at_yn =
                get_tan(deps_aligned, sig_compression, q_is,
                        check6c.parms, /*istrain=*/1);
            expect_eq_int(tan_at_yn.error, kErrorOk,
                          "case6c reference tan_at_yn ok");

            const double m_R6c = check6c.parms[9];

            // (i) D == m_R * tan.L (perturbate_h convention).
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double expected = m_R6c * tan_at_yn.L[i][j];
                    const double tol = 1e-9 * (std::abs(expected) + 1.0);
                    char lbl[120];
                    std::snprintf(lbl, sizeof(lbl),
                                  "case6c (i) D[%d][%d] == m_R * L[%d][%d] "
                                  "(perturbate_h, NOT M)", i, j, i, j);
                    expect_near(r.D[i][j], expected, tol, lbl);
                }
            }

            // (iii) Anti-tautology: tan_at_yn.M MUST differ from
            //       m_R * tan_at_yn.L somewhere meaningfully. If they
            //       happened to coincide here (zero AA / degenerate
            //       parms), assertion (i) would still pass via (M ==
            //       m_R*L) and the regression wouldn't be caught.
            double max_rel_diff = 0.0;
            int max_i = -1, max_j = -1;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double M_val = tan_at_yn.M[i][j];
                    const double mR_L  = m_R6c * tan_at_yn.L[i][j];
                    const double scale = std::abs(M_val) + std::abs(mR_L) + 1.0;
                    const double rel   = std::abs(M_val - mR_L) / scale;
                    if (rel > max_rel_diff) {
                        max_rel_diff = rel;
                        max_i = i;
                        max_j = j;
                    }
                }
            }
            if (max_rel_diff > 1.0e-6) {
                ++g_pass;
            } else {
                ++g_fail;
                std::fprintf(stderr,
                             "[FAIL] case6c (iii) anti-tautology: tan_at_yn.M "
                             "matches m_R*L (max rel diff %.3g at [%d][%d]); "
                             "the M-vs-m_R*L distinction is degenerate so "
                             "assertion (i) does not actually catch a regression "
                             "that returned M instead of m_R*L\n",
                             max_rel_diff, max_i, max_j);
            }

            // (ii) D MUST differ from tan_at_yn.M. The killer: if the
            //      kernel ever started returning M as the final tangent
            //      (regression of the perturbate_h fix), this fires
            //      regardless of whether (i) coincidentally still holds.
            bool D_differs_from_M = false;
            for (int i = 0; i < 6 && !D_differs_from_M; ++i) {
                for (int j = 0; j < 6 && !D_differs_from_M; ++j) {
                    const double D_val = r.D[i][j];
                    const double M_val = tan_at_yn.M[i][j];
                    const double scale = std::abs(D_val) + std::abs(M_val) + 1.0;
                    if (std::abs(D_val - M_val) / scale > 1.0e-6) {
                        D_differs_from_M = true;
                    }
                }
            }
            if (D_differs_from_M) {
                ++g_pass;
            } else {
                ++g_fail;
                std::fprintf(stderr,
                             "[FAIL] case6c (ii) regression: r.D matches "
                             "tan_at_yn.M -- the kernel may have started "
                             "returning M instead of m_R*L (perturbate_h "
                             "convention violated).\n");
            }
        }
    }

    // ====================================================================
    // (7) Initial-gate fail-closed routing.
    //
    //     check_RKF has 4 failure modes ("pmean", "tension", "nan",
    //     "void"). integrate_step should ONLY treat "pmean" and
    //     "tension" as elastic-fallback (Fortran inittension) cases.
    //     Anything else is a caller-side contract violation -- silent
    //     "elastic recovery" would let corrupted state persist
    //     invisibly across step boundaries (Codex 2026-05-01 high).
    //
    //     For each invalid input, assert: error = 3, used_elastic =
    //     false, nfev = 0, sig and q reverted to inputs (which may be
    //     the corrupted values themselves -- the caller knows what
    //     they passed; we just don't compute new state from rotten
    //     input).
    // ====================================================================
    {
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};

        auto run_and_check = [&](const double sig_in[6],
                                 const double q_in[7],
                                 const char* label_prefix) {
            const StepResult r = integrate_step(sig_in, q_in, deps,
                                                parms, dtime, dtime);
            char lbl[160];
            std::snprintf(lbl, sizeof(lbl), "%s error == 3", label_prefix);
            expect_eq_int(r.error, 3, lbl);
            std::snprintf(lbl, sizeof(lbl),
                          "%s failure_class == ContractViolation", label_prefix);
            expect_failure_class(r.failure_class,
                                 FailureClass::ContractViolation, lbl);
            std::snprintf(lbl, sizeof(lbl), "%s used_elastic == false", label_prefix);
            expect_eq_bool(r.used_elastic, false, lbl);
            std::snprintf(lbl, sizeof(lbl), "%s nfev == 0", label_prefix);
            expect_eq_int(r.nfev, 0, lbl);
            // r.sig / r.q stay at inputs (which may BE the corrupted
            // values; we don't sanitize). D zero on error path.
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    char dlbl[120];
                    std::snprintf(dlbl, sizeof(dlbl),
                                  "%s D[%d][%d] == 0", label_prefix, i, j);
                    expect_near(r.D[i][j], 0.0, 0.0, dlbl);
                }
            }
        };

        // (7a) void = 0 (finite, non-positive). check_RKF flags "void".
        {
            const double q_bad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            run_and_check(sig_compression, q_bad, "case7a void=0");
        }

        // (7b) void < 0 (finite, negative).
        {
            const double q_bad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.5};
            run_and_check(sig_compression, q_bad, "case7b void<0");
        }

        // (7c) NaN in a state slot (not void). check_RKF flags "nan"
        //      via the umatisnan_h scan over y[0..12].
        {
            double q_bad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.825};
            q_bad[1] = std::numeric_limits<double>::quiet_NaN();
            run_and_check(sig_compression, q_bad, "case7c NaN-in-q[1]");
        }

        // (7d) +Inf in stress slot. check_RKF flags "nan".
        {
            double sig_bad[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
            sig_bad[0] = std::numeric_limits<double>::infinity();
            run_and_check(sig_bad, q_zero_is, "case7d +Inf-in-sig");
        }

        // (7e) Finite overflow > 1e30 in q (umatisnan_h overflow check).
        {
            double q_bad[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.825};
            q_bad[3] = 1.0e31;
            run_and_check(sig_compression, q_bad, "case7e |q|>1e30");
        }

        // (7f) Sanity: an actual inittension-style failure (tensile
        //      mean stress) MUST still go through the elastic
        //      fallback, NOT the new fail-closed path. This locks the
        //      pmean/tension exemption -- if someone narrows the
        //      exemption later, the existing case 2 also fails.
        {
            const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                                deps, parms, dtime, dtime);
            expect_eq_int(r.error, kErrorOk,
                          "case7f tensile sig (pmean failure) -> error=0 (elastic)");
            expect_eq_bool(r.used_elastic, true,
                           "case7f tensile sig -> used_elastic=true (NOT fail-closed)");
        }
    }

    // ====================================================================
    // (8) Scalar tuning parameter contract.
    //
    //     The wrapper validates dtime, dtsub_in, err_tol, DTmin at
    //     entry, BEFORE the y_init pre-flight and the check_RKF gate.
    //     The elastic-fallback branch bypasses rkf23_update entirely,
    //     so without this entry-side check an invalid scalar would
    //     leak through the elastic path and corrupt r.dtsub_next:
    //       - dtsub_in = NaN: clamp `<= 0` and `>= dtime` are both
    //         false on NaN; r.dtsub_next stays NaN.
    //       - dtime < 0: any positive dtsub_in satisfies `>= dtime`
    //         and r.dtsub_next gets pinned to the negative dtime.
    //
    //     For each invalid scalar, the wrapper must:
    //       - return error = 3
    //       - leave r.sig at sig_n, r.q at q_n (no work done)
    //       - set r.dtsub_next = 0 explicitly (the bottom clamp can
    //         NOT sanitize NaN, so the entry-side guard handles this
    //         directly)
    //       - set used_elastic = false, nfev = 0, D = 0
    //
    //     Codex review 2026-05-01 (medium follow-up).
    // ====================================================================
    {
        // Tensile sig used to confirm the scalar guard FIRES BEFORE
        // the elastic-fallback decision. Without the guard, this
        // tensile sig would have routed to elastic and the corrupt
        // scalar would have leaked through.
        const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double deps[6]        = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        auto check_scalar_rejected =
            [&](double dt, double dts, double et, double DTmin_in,
                const char* lbl) {
                const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                                    deps, parms,
                                                    dt, dts, et,
                                                    /*maxnint=*/100,
                                                    DTmin_in);
                char buf[160];
                std::snprintf(buf, sizeof(buf), "%s error == 3", lbl);
                expect_eq_int(r.error, 3, buf);
                std::snprintf(buf, sizeof(buf),
                              "%s failure_class == ContractViolation", lbl);
                expect_failure_class(r.failure_class,
                                     FailureClass::ContractViolation, buf);
                std::snprintf(buf, sizeof(buf), "%s nfev == 0", lbl);
                expect_eq_int(r.nfev, 0, buf);
                std::snprintf(buf, sizeof(buf), "%s used_elastic == false", lbl);
                expect_eq_bool(r.used_elastic, false, buf);
                // r.dtsub_next must be a SAFE finite value (0), not
                // NaN. NaN-poisoned dtsub_next would survive into the
                // caller's persisted state and corrupt the next step.
                std::snprintf(buf, sizeof(buf),
                              "%s dtsub_next == 0 (NOT NaN/negative)", lbl);
                expect_near(r.dtsub_next, 0.0, 0.0, buf);
                // sig and q stay at inputs.
                for (int i = 0; i < 6; ++i) {
                    char slbl[80];
                    std::snprintf(slbl, sizeof(slbl),
                                  "%s r.sig[%d] preserved", lbl, i);
                    expect_near(r.sig[i], sig_tensile[i], 0.0, slbl);
                }
            };

        const double NaN_d = std::numeric_limits<double>::quiet_NaN();
        const double Inf_d = std::numeric_limits<double>::infinity();

        // (8a) dtsub_in = NaN. Tensile sig would have gone elastic;
        //      scalar guard MUST fire first.
        check_scalar_rejected(/*dt=*/1.0, /*dts=*/NaN_d,
                              /*et=*/1.0e-3, /*DTmin=*/1.0e-17,
                              "case8a dtsub_in=NaN+tensile");

        // (8b) dtime < 0. Same setup.
        check_scalar_rejected(/*dt=*/-1.0, /*dts=*/0.5,
                              /*et=*/1.0e-3, /*DTmin=*/1.0e-17,
                              "case8b dtime<0+tensile");

        // (8c) dtime = 0.
        check_scalar_rejected(/*dt=*/0.0, /*dts=*/0.5,
                              /*et=*/1.0e-3, /*DTmin=*/1.0e-17,
                              "case8c dtime=0");

        // (8d) err_tol = NaN.
        check_scalar_rejected(/*dt=*/1.0, /*dts=*/0.5,
                              /*et=*/NaN_d, /*DTmin=*/1.0e-17,
                              "case8d err_tol=NaN");

        // (8e) DTmin <= 0.
        check_scalar_rejected(/*dt=*/1.0, /*dts=*/0.5,
                              /*et=*/1.0e-3, /*DTmin=*/0.0,
                              "case8e DTmin=0");

        // (8f) dtsub_in = +Inf.
        check_scalar_rejected(/*dt=*/1.0, /*dts=*/Inf_d,
                              /*et=*/1.0e-3, /*DTmin=*/1.0e-17,
                              "case8f dtsub_in=+Inf");

        // (8g) maxnint <= 0 -> ContractViolation, NOT RkfReject.
        //      Codex 2026-05-02 P2 finding: integrate_step previously
        //      did NOT validate maxnint, so rkf23_update would return
        //      error=3 on the first ksubst > maxnint check, and the
        //      wrapper would classify it as FailureClass::RkfReject
        //      (telling the bridge to cutback). A smaller outer dtime
        //      cannot fix an invalid substep limit -- bridge would
        //      loop forever. The fix adds maxnint > 0 to the entry
        //      scalar contract (cpp:1217-1233), routing this to
        //      ContractViolation.
        //
        //      Use compressive sig + small deps so RKF would actually
        //      run if maxnint were valid. Failure must come from the
        //      entry guard, not from rkf23 internals.
        {
            const double sig_compression3[6] = {-100.0, -100.0, -100.0,
                                                0.0, 0.0, 0.0};
            const double deps_small[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                          0.0, 0.0, 0.0};

            auto check_maxnint_rejected = [&](int bad_maxnint, const char* lbl) {
                const StepResult r = integrate_step(sig_compression3, q_zero_is,
                                                    deps_small, parms,
                                                    /*dtime=*/1.0,
                                                    /*dtsub_in=*/0.5,
                                                    /*err_tol=*/1.0e-3,
                                                    bad_maxnint,
                                                    /*DTmin=*/1.0e-17);
                char buf[160];
                std::snprintf(buf, sizeof(buf), "%s error == 3", lbl);
                expect_eq_int(r.error, 3, buf);
                std::snprintf(buf, sizeof(buf),
                              "%s failure_class == ContractViolation (NOT RkfReject)",
                              lbl);
                expect_failure_class(r.failure_class,
                                     FailureClass::ContractViolation, buf);
                std::snprintf(buf, sizeof(buf), "%s nfev == 0 (rejected at entry)", lbl);
                expect_eq_int(r.nfev, 0, buf);
                std::snprintf(buf, sizeof(buf), "%s used_elastic == false", lbl);
                expect_eq_bool(r.used_elastic, false, buf);
                std::snprintf(buf, sizeof(buf), "%s dtsub_next == 0", lbl);
                expect_near(r.dtsub_next, 0.0, 0.0, buf);
            };

            check_maxnint_rejected(/*maxnint=*/0,  "case8g maxnint=0");
            check_maxnint_rejected(/*maxnint=*/-1, "case8g maxnint=-1");
        }
    }

    // ====================================================================
    // (9) Strain increment contract: deps_np1 must be finite + |v| <= 1e30.
    //
    //     The most consequential failure path Codex flagged: tensile
    //     initial sig + NaN deps. Without the deps guard:
    //       - sig is tensile -> check_RKF flags "pmean"
    //       - is_inittension = true -> elastic fallback
    //       - calc_elasti(NaN deps) -> matmul propagates NaN into d_sig
    //       - r.sig = sig_n + NaN = NaN
    //       - return error=0, used_elastic=true, r.sig contains NaN
    //     With the guard, the deps NaN is caught at entry, returning
    //     error=3 BEFORE any branching decision.
    //
    //     Cover both branches:
    //       (a-c) tensile sig + bad deps: would have hit elastic path.
    //       (d) healthy compressive sig + bad deps: would have hit RKF.
    //           Even though the RKF chain reaction would catch deps NaN
    //           via check_RKF mid-substep, we want UNIFORM contract
    //           enforcement (caller can't tell which path will fire).
    // ====================================================================
    {
        const double NaN_d = std::numeric_limits<double>::quiet_NaN();
        const double Inf_d = std::numeric_limits<double>::infinity();
        const double sig_tensile[6]     = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double sig_compression2[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};

        auto check_deps_rejected =
            [&](const double sig_in[6], const double deps_in[6],
                const char* lbl) {
                const StepResult r = integrate_step(sig_in, q_zero_is,
                                                    deps_in, parms,
                                                    dtime, dtime);
                char buf[160];
                std::snprintf(buf, sizeof(buf), "%s error == 3", lbl);
                expect_eq_int(r.error, 3, buf);
                std::snprintf(buf, sizeof(buf),
                              "%s failure_class == ContractViolation", lbl);
                expect_failure_class(r.failure_class,
                                     FailureClass::ContractViolation, buf);
                std::snprintf(buf, sizeof(buf), "%s nfev == 0", lbl);
                expect_eq_int(r.nfev, 0, buf);
                std::snprintf(buf, sizeof(buf), "%s used_elastic == false", lbl);
                expect_eq_bool(r.used_elastic, false, buf);
                // The killer assertion: r.sig must NOT contain NaN.
                // Without the deps guard, the elastic path would have
                // matmul'd the NaN deps into d_sig -> r.sig.
                for (int i = 0; i < 6; ++i) {
                    char slbl[80];
                    std::snprintf(slbl, sizeof(slbl),
                                  "%s r.sig[%d] finite (preserved)", lbl, i);
                    if (std::isfinite(r.sig[i]) && r.sig[i] == sig_in[i]) {
                        ++g_pass;
                    } else {
                        ++g_fail;
                        std::fprintf(stderr,
                                     "[FAIL] %s: got %g, expected %g (finite)\n",
                                     slbl, r.sig[i], sig_in[i]);
                    }
                }
            };

        // (9a) tensile sig + deps[0] = NaN -- the exact path Codex described.
        {
            double deps[6] = {NaN_d, 0.0, 0.0, 0.0, 0.0, 0.0};
            check_deps_rejected(sig_tensile, deps,
                                "case9a tensile sig + deps[0]=NaN (elastic path bypass)");
        }

        // (9b) tensile sig + deps[3] = +Inf (off-diagonal).
        {
            double deps[6] = {0.0, 0.0, 0.0, Inf_d, 0.0, 0.0};
            check_deps_rejected(sig_tensile, deps,
                                "case9b tensile sig + deps[3]=+Inf");
        }

        // (9c) tensile sig + deps[2] = 1e31 (overflow boundary).
        {
            double deps[6] = {0.0, 0.0, 1.0e31, 0.0, 0.0, 0.0};
            check_deps_rejected(sig_tensile, deps,
                                "case9c tensile sig + deps[2]=1e31");
        }

        // (9d) Healthy compressive sig + NaN deps. RKF main path would
        //      have caught this via check_RKF mid-substep, but the
        //      contract should fail at entry uniformly.
        {
            double deps[6] = {NaN_d, NaN_d, NaN_d, 0.0, 0.0, 0.0};
            check_deps_rejected(sig_compression2, deps,
                                "case9d compressive sig + all-NaN deps (RKF path uniform)");
        }

        // (9e) Boundary sanity: deps = exactly 1e30 should be accepted
        //      (umatisnan_h-style strict `>`). Use compressive sig and
        //      small magnitude in other slots so the run can succeed.
        //      Pure check that the boundary is `>` not `>=`.
        {
            double deps[6] = {1.0e30, 0.0, 0.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression2, q_zero_is,
                                                deps, parms, dtime, dtime);
            // We only assert the deps guard did NOT fire at entry.
            // The actual physics may reject (RKF cutback, tangent
            // failure, etc.) for such an extreme deps -- that's fine,
            // just NOT the entry guard.
            // Equivalently: error must NOT be 3 with nfev = 0 caused
            // by the entry guard. nfev > 0 OR error != 3 OR used_elastic
            // both indicate the entry guard let it through.
            const bool guard_did_not_fire =
                (r.error != 3) || (r.nfev > 0) || r.used_elastic;
            if (guard_did_not_fire) {
                ++g_pass;
            } else {
                ++g_fail;
                std::fprintf(stderr,
                             "[FAIL] case9e deps=1e30 boundary: entry guard fired "
                             "(error=%d, nfev=%d, used_elastic=%d) but should not "
                             "(strict > 1e30)\n",
                             r.error, r.nfev, static_cast<int>(r.used_elastic));
            }
        }
    }

    // ====================================================================
    // (10) Final-tangent state regression: nonzero deps so y_np1 != y_n.
    //
    //      Bug being locked (Codex 2026-05-02 high): integrate_step used
    //      to call get_tan(deps_np1, r.sig, r.q, ...) -- the post-RKF
    //      state -- to compute the final tangent. Fortran perturbate_h
    //      .for:1466-1469 explicitly extracts sig/q from y_n (begin of
    //      step), not y_np1.
    //
    //      Case (6) used zero deps, so y_n == y_np1 and could not
    //      distinguish the two states. This case picks deps small enough
    //      that RKF accepts in one substep (so the test stays simple)
    //      but large enough that |L(y_n) - L(y_np1)| is well above
    //      floating noise.
    //
    //      Two assertions:
    //        (i)  r.D matches m_R * get_tan(deps, sig_n, q_n, IS).L
    //             (textbook m_R = 1, so D == L_at_yn directly).
    //        (ii) get_tan(deps, sig_n, q_n).L differs from
    //             get_tan(deps, r.sig, r.q).L by > 1e-6 relative on
    //             at least one component. This proves the test isn't
    //             tautological -- the two states ARE materially
    //             different. Without (ii), a regression that swaps
    //             the state back could still pass (i) if the two
    //             tangents happened to coincide.
    // ====================================================================
    {
        // Compressive deps that the textbook RKF accepts cleanly. We
        // checked offline: this magnitude advances sig by ~5 kPa, which
        // is well above the FP precision of get_tan inputs.
        const double deps[6] = {-1.0e-4, -1.0e-4, -1.0e-4, 0.0, 0.0, 0.0};

        const StepResult r = integrate_step(sig_compression, q_zero_is,
                                            deps, parms, dtime, dtime);
        expect_eq_int(r.error, kErrorOk, "case10 RKF main path success");
        expect_eq_bool(r.used_elastic, false, "case10 used_elastic == false");

        // Sanity: state actually advanced.
        bool state_advanced = false;
        for (int i = 0; i < 3; ++i) {
            if (std::abs(r.sig[i] - sig_compression[i]) > 1.0) {
                state_advanced = true;
                break;
            }
        }
        if (state_advanced) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case10 sanity: state did not advance meaningfully "
                         "(sig_n[0]=%g, sig_new[0]=%g)\n",
                         sig_compression[0], r.sig[0]);
        }

        // Reference tangents at both states. textbook m_R = 1 -> istrain = 1.
        const TangentResult tan_at_yn = get_tan(deps, sig_compression,
                                                q_zero_is, parms,
                                                /*istrain=*/1);
        expect_eq_int(tan_at_yn.error, kErrorOk,
                      "case10 reference get_tan(y_n) ok");

        const TangentResult tan_at_ynp1 = get_tan(deps, r.sig, r.q,
                                                  parms, /*istrain=*/1);
        expect_eq_int(tan_at_ynp1.error, kErrorOk,
                      "case10 reference get_tan(y_np1) ok");

        // (i) r.D must match L at y_n (m_R = 1).
        const double m_R = parms[9];
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double expected = m_R * tan_at_yn.L[i][j];
                const double tol = 1e-9 * (std::abs(expected) + 1.0);
                char label[120];
                std::snprintf(label, sizeof(label),
                              "case10(i) D[%d][%d] == m_R * L_at_yn[%d][%d]",
                              i, j, i, j);
                expect_near(r.D[i][j], expected, tol, label);
            }
        }

        // (ii) Anti-tautology: L(y_n) and L(y_np1) must differ
        //      meaningfully. If they coincide, the test above is empty.
        double max_rel_diff = 0.0;
        int max_i = -1, max_j = -1;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double a = tan_at_yn.L[i][j];
                const double b = tan_at_ynp1.L[i][j];
                const double scale = std::abs(a) + std::abs(b) + 1.0;
                const double rel = std::abs(a - b) / scale;
                if (rel > max_rel_diff) {
                    max_rel_diff = rel;
                    max_i = i;
                    max_j = j;
                }
            }
        }
        if (max_rel_diff > 1.0e-6) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case10(ii) anti-tautology: L_at_yn ~== L_at_ynp1 "
                         "(max rel diff %.3g at [%d][%d]); test (i) is empty.\n",
                         max_rel_diff, max_i, max_j);
        }

        // (iii) The killer: r.D must NOT match L at y_np1. If a regression
        //       swaps the state argument, this assertion will trip even
        //       if (ii) somehow degenerates.
        bool D_differs_from_ynp1 = false;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double diff = std::abs(r.D[i][j] - m_R * tan_at_ynp1.L[i][j]);
                const double scale = std::abs(r.D[i][j]) +
                                     std::abs(m_R * tan_at_ynp1.L[i][j]) + 1.0;
                if (diff / scale > 1.0e-6) {
                    D_differs_from_ynp1 = true;
                    break;
                }
            }
            if (D_differs_from_ynp1) break;
        }
        if (D_differs_from_ynp1) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case10(iii) bug regression: r.D matches L_at_ynp1, "
                         "which is the WRONG state. Final tangent must use sig_n/q_n.\n");
        }
    }

    // ====================================================================
    // (11) Strain norm overflow guard: deps with each component below the
    //      per-component sentinel but ||deps||_strain above 1e30.
    //
    //      Setup: deps = (8e29, 8e29, 0, 0, 0, 0).
    //        Per-component: max |v| = 8e29 < 1e30 -> per-component check passes.
    //        StrainLike norm: sqrt((8e29)^2 + (8e29)^2 + 0 + 0 + 0 + 0)
    //                       = sqrt(1.28e60) ~= 1.131e30 > 1e30 -> norm guard fires.
    //
    //      Bug being locked: without the norm guard, this exact deps with
    //      tensile initial sig (which routes to elastic fallback) would
    //      slip past the entry contract, calc_elasti would matmul those
    //      8e29 components into d_sig (~ E * 8e29 ~ 1e32 magnitude),
    //      and r.sig would be returned overflow-scale with error = 0,
    //      used_elastic = true. Mirrors Fortran umat .for:201-209.
    //
    //      Three sub-checks:
    //        (a) Tensile sig + overflow-norm deps -> error = 3 (NOT
    //            elastic success). The killer assertion.
    //        (b) Compressive sig + overflow-norm deps -> error = 3.
    //            Uniform contract regardless of which branch would
    //            have run.
    //        (c) Boundary sanity: deps = (1e30, 0, 0, 0, 0, 0) has
    //            norm = 1e30 exactly. fails_umatisnan_h uses strict
    //            `>`, so the entry guard must NOT fire here. (The
    //            run may still error out downstream on physics; we
    //            only check the entry guard didn't trip.)
    // ====================================================================
    {
        const double sig_tensile[6]      = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double sig_compression3[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};

        // (11a) tensile + overflow-norm deps. The path Codex described.
        {
            const double deps[6] = {8.0e29, 8.0e29, 0.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                                deps, parms, dtime, dtime);
            expect_eq_int(r.error, 3,
                          "case11a tensile+norm-overflow deps -> error=3");
            expect_failure_class(r.failure_class,
                                 FailureClass::ContractViolation,
                                 "case11a failure_class == ContractViolation "
                                 "(NOT RkfReject -- bridge MUST escalate, not cutback)");
            expect_eq_int(r.nfev, 0, "case11a nfev == 0");
            expect_eq_bool(r.used_elastic, false,
                           "case11a used_elastic == false (NOT corrupted elastic)");
            // The killer assertion: r.sig must NOT be overflow-corrupted.
            // Without the norm guard, calc_elasti would have produced
            // d_sig ~ 1e32 and r.sig = sig_n + d_sig ~ 1e32.
            for (int i = 0; i < 6; ++i) {
                char lbl[80];
                std::snprintf(lbl, sizeof(lbl),
                              "case11a r.sig[%d] preserved (no overflow corruption)", i);
                if (std::isfinite(r.sig[i]) && std::abs(r.sig[i]) < 1.0e15
                    && r.sig[i] == sig_tensile[i]) {
                    ++g_pass;
                } else {
                    ++g_fail;
                    std::fprintf(stderr,
                                 "[FAIL] %s: got %g, expected %g\n",
                                 lbl, r.sig[i], sig_tensile[i]);
                }
            }
        }

        // (11b) compressive + overflow-norm deps. Uniform contract.
        {
            const double deps[6] = {8.0e29, 8.0e29, 0.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression3, q_zero_is,
                                                deps, parms, dtime, dtime);
            expect_eq_int(r.error, 3,
                          "case11b compressive+norm-overflow deps -> error=3");
            expect_failure_class(r.failure_class,
                                 FailureClass::ContractViolation,
                                 "case11b failure_class == ContractViolation "
                                 "(uniform contract regardless of stress sign)");
            expect_eq_int(r.nfev, 0, "case11b nfev == 0");
            expect_eq_bool(r.used_elastic, false,
                           "case11b used_elastic == false");
            for (int i = 0; i < 6; ++i) {
                char lbl[80];
                std::snprintf(lbl, sizeof(lbl),
                              "case11b r.sig[%d] preserved", i);
                expect_near(r.sig[i], sig_compression3[i], 0.0, lbl);
            }
        }

        // (11c) Boundary sanity: norm = exactly 1e30 must NOT trip the
        //       entry guard (umatisnan_h uses strict `>`). Use a single
        //       1e30 component with all others zero, so norm == 1e30.
        {
            const double deps[6] = {1.0e30, 0.0, 0.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression3, q_zero_is,
                                                deps, parms, dtime, dtime);
            // We only assert: the entry guard did NOT fire on the norm
            // (strict `>` semantics). Physics may reject downstream;
            // that's fine. "Entry guard fired" looks like:
            //   error = 3 AND nfev = 0 AND used_elastic = false.
            // Equivalently, anything else means the guard let it through.
            const bool entry_guard_did_not_fire =
                (r.error != 3) || (r.nfev > 0) || r.used_elastic;
            if (entry_guard_did_not_fire) {
                ++g_pass;
            } else {
                ++g_fail;
                std::fprintf(stderr,
                             "[FAIL] case11c norm=1e30 boundary: entry guard "
                             "fired (error=%d, nfev=%d, used_elastic=%d) but "
                             "should not (strict > 1e30)\n",
                             r.error, r.nfev, static_cast<int>(r.used_elastic));
            }
        }
    }

    // ====================================================================
    // (12) First-call dtsub normalization. Mirrors Fortran umat
    //      .for:243-244: finite `dtsub <= 0` or `dtsub > dtime` is
    //      silently set to `dtime` before integration. The previous
    //      kernel rejected dtsub_in <= 0 as a contract violation, which
    //      would deadlock any bridge that persists dtsub in zero-init
    //      statev (Codex 2026-05-02 high finding -- "first call with
    //      zero dtsub is rejected instead of initialized").
    // ====================================================================
    {
        const double sig_tensile[6]      = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        const double sig_compression4[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        const double deps_zero[6]        = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const double deps_small[6]       = {-1.0e-3, -1.0e-3, -1.0e-3,
                                            0.0, 0.0, 0.0};

        // (12a) compressive sig + dtsub_in = 0 -> RKF main path success.
        //       Reference: same call with dtsub_in = dtime must produce
        //       identical sig, q, D and dtsub_next (proves the
        //       normalization is exactly `dtsub_in = dtime`, not some
        //       other clamp like 0.5*dtime).
        {
            const StepResult r0 = integrate_step(sig_compression4, q_zero_is,
                                                 deps_zero, parms,
                                                 dtime, /*dtsub_in=*/0.0);
            expect_eq_int(r0.error, kErrorOk,
                          "case12a dtsub_in=0 + compressive -> error=0");
            expect_eq_bool(r0.used_elastic, false,
                           "case12a used_elastic == false (RKF main path)");
            // Reference: explicit dtsub_in = dtime.
            const StepResult ref = integrate_step(sig_compression4, q_zero_is,
                                                  deps_zero, parms,
                                                  dtime, dtime);
            expect_eq_int(ref.error, kErrorOk, "case12a reference run ok");
            expect_near(r0.dtsub_next, ref.dtsub_next, 1e-15,
                        "case12a dtsub_in=0 dtsub_next == reference");
            for (int i = 0; i < 6; ++i) {
                char lbl[80];
                std::snprintf(lbl, sizeof(lbl),
                              "case12a r.sig[%d] == reference", i);
                expect_near(r0.sig[i], ref.sig[i], 1e-12, lbl);
            }
            expect_eq_int(r0.nfev, ref.nfev, "case12a nfev == reference");
        }

        // (12b) tensile sig + dtsub_in = 0 -> elastic fallback success.
        //       The exact path Codex described: bridge first-call deadlock
        //       under tensile initial state. Must NOT reject.
        {
            const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                                deps_small, parms,
                                                dtime, /*dtsub_in=*/0.0);
            expect_eq_int(r.error, kErrorOk,
                          "case12b dtsub_in=0 + tensile -> error=0 (elastic)");
            expect_eq_bool(r.used_elastic, true,
                           "case12b used_elastic == true");
            expect_eq_int(r.nfev, 0, "case12b nfev == 0 (elastic path)");
            // After normalization, dtsub_in == dtime, and the elastic
            // path leaves r.dtsub_next at the input value (then the
            // bottom clamp keeps it == dtime).
            expect_near(r.dtsub_next, dtime, 1e-15,
                        "case12b r.dtsub_next == dtime (normalized)");
        }

        // (12c) dtsub_in < 0 -> normalize to dtime.
        {
            const StepResult r = integrate_step(sig_compression4, q_zero_is,
                                                deps_zero, parms,
                                                dtime, /*dtsub_in=*/-1.0);
            expect_eq_int(r.error, kErrorOk,
                          "case12c dtsub_in=-1 -> error=0 (normalized)");
            expect_eq_bool(r.used_elastic, false,
                           "case12c used_elastic == false");
        }

        // (12d) dtsub_in > dtime -> normalize to dtime.
        //       Note: rkf23_update has its own internal clamp for this,
        //       but integrate_step's outer normalization should fire
        //       first (independent of the elastic path).
        {
            const StepResult r = integrate_step(sig_compression4, q_zero_is,
                                                deps_zero, parms,
                                                dtime, /*dtsub_in=*/2.0 * dtime);
            expect_eq_int(r.error, kErrorOk,
                          "case12d dtsub_in=2*dtime -> error=0 (normalized)");
            expect_eq_bool(r.used_elastic, false,
                           "case12d used_elastic == false");
        }

        // (12e) Boundary check: NaN dtsub_in MUST still reject (delta
        //       vs the normalize cases above). Locks "non-finite stays
        //       hard-rejected after the first-call normalization
        //       relaxation".
        {
            const double NaN_d = std::numeric_limits<double>::quiet_NaN();
            const StepResult r = integrate_step(sig_compression4, q_zero_is,
                                                deps_zero, parms,
                                                dtime, /*dtsub_in=*/NaN_d);
            expect_eq_int(r.error, 3,
                          "case12e dtsub_in=NaN STILL rejected (no normalize)");
            expect_failure_class(r.failure_class, FailureClass::ContractViolation,
                                 "case12e failure_class == ContractViolation (NaN dtsub)");
            expect_near(r.dtsub_next, 0.0, 0.0,
                        "case12e dtsub_next == 0 on NaN reject");
        }

        // (12f) RKF-path coverage: dtsub_in normalization with NONZERO
        //       deps MUST actually flow into rkf23_update (not the
        //       zero-deps short circuit added by case 1). Cases 12a/c/d
        //       all use deps_zero, so post the zero-deps short-circuit
        //       fix the normalized dtsub never reaches rkf23 there. This
        //       sub-case feeds compressive deps so rkf23 does run, then
        //       compares the {dtsub_in=0, dtsub_in=-1, dtsub_in=2*dtime}
        //       runs against a {dtsub_in=dtime} reference -- proving
        //       the normalization actually feeds dtime into the RKF
        //       integrator, not just into the elastic / short-circuit
        //       branches. Codex review 2026-05-02 (P2 finding:
        //       dtsub-normalization-on-RKF-path coverage gap).
        {
            const double deps_compress[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                             0.0, 0.0, 0.0};

            // Reference: explicit dtsub_in = dtime, RKF success.
            const StepResult ref = integrate_step(sig_compression4, q_zero_is,
                                                  deps_compress, parms,
                                                  dtime, dtime);
            expect_eq_int(ref.error, kErrorOk,
                          "case12f reference run (dtsub_in=dtime, nonzero deps) ok");
            expect_eq_bool(ref.used_elastic, false,
                           "case12f reference used_elastic == false");
            // Sanity: nfev > 0 -> rkf23 actually entered. If this is
            // ever 0 the test setup degenerated to the zero-deps short
            // circuit again and the rest of the case is empty.
            if (ref.nfev > 0) ++g_pass;
            else { ++g_fail; std::fprintf(stderr,
                "[FAIL] case12f reference nfev > 0 (rkf23 entered): got %d\n",
                ref.nfev); }

            auto check_normalized_run = [&](double dtsub_attack,
                                            const char* label_prefix) {
                const StepResult r = integrate_step(sig_compression4, q_zero_is,
                                                    deps_compress, parms,
                                                    dtime, dtsub_attack);
                char lbl[160];
                std::snprintf(lbl, sizeof(lbl), "%s error == 0 (normalized + RKF)",
                              label_prefix);
                expect_eq_int(r.error, kErrorOk, lbl);
                std::snprintf(lbl, sizeof(lbl), "%s used_elastic == false", label_prefix);
                expect_eq_bool(r.used_elastic, false, lbl);

                // The killer assertion: nfev MUST match the reference.
                // If integrate_step had handed an unnormalized dtsub_in
                // (e.g. left it at 0 or -1) into rkf23, rkf23's own
                // internal clamp would still rescue the run, but the
                // resulting nfev pattern would differ. Matching nfev
                // proves the normalization happened at the wrapper
                // before rkf23 saw it.
                std::snprintf(lbl, sizeof(lbl),
                              "%s nfev == reference (proves dtsub normalized at wrapper)",
                              label_prefix);
                expect_eq_int(r.nfev, ref.nfev, lbl);

                // sig and dtsub_next match reference bit-tight.
                for (int i = 0; i < 6; ++i) {
                    char slbl[120];
                    std::snprintf(slbl, sizeof(slbl), "%s r.sig[%d] == reference",
                                  label_prefix, i);
                    expect_near(r.sig[i], ref.sig[i], 1e-12, slbl);
                }
                std::snprintf(lbl, sizeof(lbl), "%s r.dtsub_next == reference",
                              label_prefix);
                expect_near(r.dtsub_next, ref.dtsub_next, 1e-15, lbl);
            };

            // (12f-i)  dtsub_in = 0 (the bridge first-call default)
            check_normalized_run(/*dtsub_in=*/0.0,
                                 "case12f-i dtsub_in=0+RKF");
            // (12f-ii) dtsub_in = -1.0
            check_normalized_run(/*dtsub_in=*/-1.0,
                                 "case12f-ii dtsub_in=-1+RKF");
            // (12f-iii) dtsub_in = 2*dtime (above-range)
            check_normalized_run(/*dtsub_in=*/2.0 * dtime,
                                 "case12f-iii dtsub_in=2*dtime+RKF");
        }
    }

    // ====================================================================
    // (13) error=10 rollback contract: rkf23_update emits error=10,
    //      integrate_step must propagate with state revert + D=0.
    //
    //      Setup: bypass check_parms by mutating parms[11]=0 on a
    //      validated Parms16. This is the same pattern test_rkf23
    //      case 4, test_F_sig_q case 4, and test_rhs case 2 use to
    //      cover get_tan's r_uc belt-and-suspenders guard at the
    //      lower-level entry points. integrate_step needs its own
    //      direct coverage because the composer-level rollback
    //      (sig=sig_n, q=q_n, D=0) is independent of how rkf23_update
    //      itself behaves.
    //
    //      Run with healthy compressive sig + small compressive deps
    //      so rkf23 actually invokes rhs (and through it get_tan,
    //      which trips the r_uc guard on the IS branch). Tensile sig
    //      would route to elastic fallback and bypass rkf23 entirely.
    //
    //      Lock the full output contract: error=10, sig==sig_n,
    //      q==q_n, all D entries zero, used_elastic=false, nfev > 0
    //      (rkf23 entered before failing -- distinguishes this from
    //      the input-corrupt path where nfev==0).
    //
    //      Codex review 2026-05-02 (P2 finding: error=10 rollback
    //      coverage gap).
    // ====================================================================
    {
        Parms16 parms_runc_zero = textbook_parms();
        parms_runc_zero[11] = 0.0;          // bypass check_parms cross-check

        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        const StepResult r = integrate_step(sig_compression, q_zero_is,
                                            deps, parms_runc_zero,
                                            dtime, dtime);
        expect_eq_int(r.error, 10,
                      "case13 rkf23-side error=10 propagated");
        expect_failure_class(r.failure_class, FailureClass::Fatal,
                             "case13 failure_class == Fatal (KRATOS_ERROR escalation)");
        expect_eq_bool(r.used_elastic, false,
                       "case13 used_elastic == false (RKF entered, NOT elastic)");

        // nfev must be > 0: distinguishes this rollback from the
        // input-corrupt path (which fails with nfev==0 before any
        // rhs call). Locks "rkf23 actually entered before failing".
        if (r.nfev > 0) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case13 nfev > 0 (rkf23 entered): got nfev=%d\n",
                         r.nfev);
        }

        // sig and q reverted to inputs.
        for (int i = 0; i < 6; ++i) {
            char lbl[80];
            std::snprintf(lbl, sizeof(lbl),
                          "case13 r.sig[%d] reverted to sig_n on error=10", i);
            expect_near(r.sig[i], sig_compression[i], 0.0, lbl);
        }
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char lbl[80];
            std::snprintf(lbl, sizeof(lbl),
                          "case13 r.q[%zu] reverted to q_n on error=10", i);
            expect_near(r.q[i], q_zero_is[i], 0.0, lbl);
        }

        // D must be all zero on the error=10 path.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char lbl[80];
                std::snprintf(lbl, sizeof(lbl),
                              "case13 r.D[%d][%d] == 0 on error=10", i, j);
                expect_near(r.D[i][j], 0.0, 0.0, lbl);
            }
        }
    }

    // ====================================================================
    // (14) Derived-stress-overflow regression. The exact attack input
    //      Codex described in the 2026-05-02 P1 high finding.
    //
    //      Pre-fix path:
    //        sig=[0,0,0,1e30,1e30,1e30] -> per-component boundary check
    //        passes -> check_RKF: pmean ok, principal_stresses_3 emits
    //        s[2]~2e30 -> "tension" set first -> umatisnan_h(min_principal)
    //        triggers but failed-string preserved as "tension" ->
    //        integrate_step routes to elastic -> r.sig overflow-corrupted,
    //        error=0, used_elastic=true.
    //
    //      Post-fix path (corruption-class precedence in check_RKF):
    //        same up through derivation; then "nan" UNCONDITIONALLY
    //        overrides "tension" -> integrate_step's
    //        is_inittension classifier sees "nan" (not pmean/tension)
    //        -> fail-closed branch -> error=3, sig preserved.
    //
    //      Killer assertions:
    //        (a) error == 3 (NOT 0).
    //        (b) used_elastic == false (NOT silent elastic).
    //        (c) r.sig == sig_attack (preserved, even though it was
    //            corrupt -- caller's responsibility, we don't sanitize
    //            inputs we just refuse to compute on them).
    //        (d) D entries all finite, NOT propagating overflow into
    //            the tangent. (D should be all zero on the error=3 path.)
    //
    //      Boundary sanity: a "near-attack" input with shears at 5e29
    //      (each below the per-component sentinel) should produce
    //      |min_principal| around 1e30 -- the test in test_check_rkf
    //      R13 already locks the check_RKF-side derivation; here we
    //      only assert the integrate_step-side composition.
    // ====================================================================
    {
        const double sig_attack[6] = {0.0, 0.0, 0.0,
                                      1.0e30, 1.0e30, 1.0e30};
        const double deps_small_sane[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                           0.0, 0.0, 0.0};

        const StepResult r = integrate_step(sig_attack, q_zero_is,
                                            deps_small_sane, parms,
                                            dtime, dtime);

        expect_eq_int(r.error, 3,
                      "case14 derived-overflow attack -> error=3 (fail-closed)");
        expect_failure_class(r.failure_class, FailureClass::ContractViolation,
                             "case14 failure_class == ContractViolation (NOT retryable, "
                             "DISTINCT from case3 RKF reject despite shared error=3)");
        expect_eq_bool(r.used_elastic, false,
                       "case14 used_elastic == false (NOT silently elastic)");

        // r.sig preserved at the corrupt input (we don't sanitize).
        for (int i = 0; i < 6; ++i) {
            char lbl[80];
            std::snprintf(lbl, sizeof(lbl),
                          "case14 r.sig[%d] preserved at input value", i);
            expect_near(r.sig[i], sig_attack[i], 0.0, lbl);
        }

        // D must be all zero (NOT overflow-propagated).
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char lbl[80];
                std::snprintf(lbl, sizeof(lbl),
                              "case14 r.D[%d][%d] == 0 (no overflow leak)", i, j);
                expect_near(r.D[i][j], 0.0, 0.0, lbl);
            }
        }

        // Sanity contrast: a healthy compressive sig with the same
        // small deps must NOT trip the fail-closed path -- proves the
        // attack-rejection is specific to the overflow case, not a
        // blanket reject of small deps.
        const double sig_healthy[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        const StepResult r_ok = integrate_step(sig_healthy, q_zero_is,
                                               deps_small_sane, parms,
                                               dtime, dtime);
        expect_eq_int(r_ok.error, kErrorOk,
                      "case14 sanity: healthy sig + same deps -> error=0");
        expect_eq_bool(r_ok.used_elastic, false,
                       "case14 sanity: healthy sig used_elastic == false");
    }

    // ====================================================================
    // (15) FailureClass disambiguation. Exercise all four enum values
    //      from one block to demonstrate the bridge-routing contract:
    //
    //        Ok                 -> commit
    //        RkfReject          -> caller cuts back outer dtime
    //        ContractViolation  -> caller escalates (KRATOS_ERROR)
    //        Fatal              -> caller escalates (KRATOS_ERROR)
    //
    //      The killer assertion is the SAME-error-code-DIFFERENT-class
    //      contrast: case 3 (huge tensile deps) and case 14 (derived
    //      overflow) both return `error == 3`, but the former is
    //      retryable (RkfReject) while the latter is not
    //      (ContractViolation). With only the `error` field a bridge
    //      cannot distinguish them.
    //
    //      Codex review 2026-05-02 (P1 high follow-up: error=3
    //      ambiguity).
    // ====================================================================
    {
        // (15a) Ok: healthy compressive sig + zero deps -> short circuit.
        {
            const double deps_zero[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression, q_zero_is,
                                                deps_zero, parms,
                                                dtime, dtime);
            expect_eq_int(r.error, kErrorOk, "case15a Ok: error == 0");
            expect_failure_class(r.failure_class, FailureClass::Ok,
                                 "case15a failure_class == Ok");
        }

        // (15b) Ok via elastic fallback: tensile sig + small deps.
        {
            const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
            const double deps_small[6]  = {-1.0e-3, -1.0e-3, -1.0e-3,
                                           0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_tensile, q_zero_is,
                                                deps_small, parms,
                                                dtime, dtime);
            expect_eq_int(r.error, kErrorOk,
                          "case15b Ok via elastic: error == 0");
            expect_failure_class(r.failure_class, FailureClass::Ok,
                                 "case15b failure_class == Ok (elastic IS success)");
            expect_eq_bool(r.used_elastic, true,
                           "case15b used_elastic == true (sanity)");
        }

        // (15c) RkfReject: huge tensile deps trips check_RKF mid-substep.
        //       error == 3, retryable.
        {
            const double deps_huge_tensile[6] = {1.0, 1.0, 1.0, 0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression, q_zero_is,
                                                deps_huge_tensile, parms,
                                                dtime, dtime);
            expect_eq_int(r.error, 3, "case15c RkfReject: error == 3");
            expect_failure_class(r.failure_class, FailureClass::RkfReject,
                                 "case15c failure_class == RkfReject (cutback OK)");
        }

        // (15d) ContractViolation: derived-overflow attack from case 14.
        //       Same error == 3 as (15c) but DIFFERENT failure_class.
        //       The killer disambiguation.
        {
            const double sig_attack[6] = {0.0, 0.0, 0.0,
                                          1.0e30, 1.0e30, 1.0e30};
            const double deps_small_sane[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                               0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_attack, q_zero_is,
                                                deps_small_sane, parms,
                                                dtime, dtime);
            expect_eq_int(r.error, 3,
                          "case15d ContractViolation: error == 3 (same as 15c!)");
            expect_failure_class(r.failure_class,
                                 FailureClass::ContractViolation,
                                 "case15d failure_class == ContractViolation "
                                 "(NOT cutback -- distinguishes from 15c via enum)");
        }

        // (15e) ContractViolation via entry pre-flight: NaN sig.
        {
            double sig_nan[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
            sig_nan[0] = std::numeric_limits<double>::quiet_NaN();
            const double deps_small[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                          0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_nan, q_zero_is,
                                                deps_small, parms,
                                                dtime, dtime);
            expect_eq_int(r.error, 3, "case15e ContractViolation: error == 3");
            expect_failure_class(r.failure_class,
                                 FailureClass::ContractViolation,
                                 "case15e failure_class == ContractViolation (NaN sig)");
        }

        // (15f) Fatal: bypass check_parms with parms[11]=0 to reach
        //       get_tan's r_uc belt-and-suspenders, which emits
        //       error == 10 from rkf23_update.
        {
            Parms16 parms_runc_zero = textbook_parms();
            parms_runc_zero[11] = 0.0;          // bypass cross-check
            const double deps_small[6] = {-1.0e-3, -1.0e-3, -1.0e-3,
                                          0.0, 0.0, 0.0};
            const StepResult r = integrate_step(sig_compression, q_zero_is,
                                                deps_small, parms_runc_zero,
                                                dtime, dtime);
            expect_eq_int(r.error, 10, "case15f Fatal: error == 10");
            expect_failure_class(r.failure_class, FailureClass::Fatal,
                                 "case15f failure_class == Fatal (escalate)");
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

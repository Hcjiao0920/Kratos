// Standalone unit-test exerciser for the Path B kernel `principal_stresses_3`
// and `check_RKF`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// Coverage:
//   principal_stresses_3:
//     (P1) Hydrostatic (diagonal): all eigenvalues equal.
//     (P2) Diagonal triaxial (different normals, no shear).
//     (P3) Mixed-sign diagonal (axisymmetric extension: +2, -1, -1).
//     (P4) Pure shear (off-diagonal only): eigenvalues +-S, 0.
//     (P5) Trace-and-sum sanity on a non-degenerate symmetric input.
//
//   check_RKF:
//     (R1) Healthy compressive state -> error=0.
//     (R2) Tensile mean stress -> error=1, failed="pmean".
//     (R3) Compressive mean but one tensile principal -> error=1,
//          failed="tension". Hand-construct a sig such that
//          tr(sig - p_t I) < 0 (so pmean check passes) but the
//          maximum eigenvalue of sig_star is positive.
//     (R4) NaN in y -> error=1, failed="nan".
//     (R5) Class precedence: pmean fails AND y has NaN. The
//          corruption class (nan/void) overrides the stress-domain
//          class (pmean/tension), so failed must be "nan".
//          Codex 2026-05-02 P1 follow-up.
//     (R6) Fortran umatisnan_h overflow guards: +Inf in q, -Inf in
//          stress, and finite 1e31 sentinel all fail as "nan".
//     (R7) Fortran umatisnan_h threshold is strict: exactly 1e30
//          remains admissible if the stress state is otherwise valid.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>

using Kratos::SandHypoCpp::CheckRkfResult;
using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::PrincipalStresses;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::check_RKF;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::kErrorOk;
using Kratos::SandHypoCpp::kYDim;
using Kratos::SandHypoCpp::principal_stresses_3;

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

void expect_failed(const CheckRkfResult& r, const char* expected,
                   const char* label) {
    if (r.failed != nullptr && std::strcmp(r.failed, expected) == 0) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr,
                     "[FAIL] %s: failed=%s, expected %s\n",
                     label,
                     r.failed ? r.failed : "(null)",
                     expected);
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
    // ====================================================================
    // principal_stresses_3
    // ====================================================================

    // (P1) Hydrostatic: eigenvalues all equal to the diagonal value.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        const PrincipalStresses ps = principal_stresses_3(sig);
        expect_near(ps.s[0], -100.0, 0.0, "P1 hydro s[0]==-100");
        expect_near(ps.s[1], -100.0, 0.0, "P1 hydro s[1]==-100");
        expect_near(ps.s[2], -100.0, 0.0, "P1 hydro s[2]==-100");
    }

    // (P2) Diagonal triaxial: eigenvalues are diagonal entries sorted.
    //      sig = (-100, -50, -25, 0, 0, 0). Sorted: (-100, -50, -25).
    {
        const double sig[6] = {-100.0, -50.0, -25.0, 0.0, 0.0, 0.0};
        const PrincipalStresses ps = principal_stresses_3(sig);
        expect_near(ps.s[0], -100.0, 0.0, "P2 triax s[0]==-100");
        expect_near(ps.s[1], -50.0,  0.0, "P2 triax s[1]==-50");
        expect_near(ps.s[2], -25.0,  0.0, "P2 triax s[2]==-25");
    }

    // (P3) Mixed-sign diagonal (axisymmetric extension): (+2, -1, -1).
    //      Sorted: (-1, -1, +2).
    {
        const double sig[6] = {2.0, -1.0, -1.0, 0.0, 0.0, 0.0};
        const PrincipalStresses ps = principal_stresses_3(sig);
        expect_near(ps.s[0], -1.0, 0.0, "P3 axi-ext s[0]==-1");
        expect_near(ps.s[1], -1.0, 0.0, "P3 axi-ext s[1]==-1");
        expect_near(ps.s[2],  2.0, 0.0, "P3 axi-ext s[2]==+2");
    }

    // (P4) Pure shear: sig = (0,0,0, 1, 0, 0). Eigenvalues: -1, 0, +1.
    {
        const double sig[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        const PrincipalStresses ps = principal_stresses_3(sig);
        expect_near(ps.s[0], -1.0, 1e-14, "P4 pure-shear s[0]==-1");
        expect_near(ps.s[1],  0.0, 1e-14, "P4 pure-shear s[1]==0");
        expect_near(ps.s[2],  1.0, 1e-14, "P4 pure-shear s[2]==+1");
    }

    // (P5) Non-degenerate trace-and-sum sanity: trace(sig) = sum(s).
    //      sig = (10, 5, 2, 1, 0.5, 0.25). Doesn't matter what the
    //      eigenvalues are individually, but sum must equal trace.
    {
        const double sig[6] = {10.0, 5.0, 2.0, 1.0, 0.5, 0.25};
        const PrincipalStresses ps = principal_stresses_3(sig);
        const double sum_s = ps.s[0] + ps.s[1] + ps.s[2];
        const double trace = sig[0] + sig[1] + sig[2];
        expect_near(sum_s, trace, 1e-13, "P5 sum(eigs) == trace(sig)");

        // Ascending order sanity.
        if (ps.s[0] <= ps.s[1] && ps.s[1] <= ps.s[2]) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] P5 ascending order: s = (%.6g, %.6g, %.6g)\n",
                         ps.s[0], ps.s[1], ps.s[2]);
        }
    }

    // ====================================================================
    // check_RKF
    // ====================================================================
    const Parms16 parms = textbook_parms();
    const double p_t       = parms[1];   // = 1.0 for textbook
    const double minstress = p_t / 4.0;  // = 0.25

    auto build_y = [](const double sig[6], double y_out[kYDim]) {
        for (int i = 0; i < 6; ++i) y_out[i] = sig[i];
        for (std::size_t i = 6; i < kYDim; ++i) y_out[i] = 0.0;
        y_out[12] = 0.825;  // void ratio
    };

    // (R1) Healthy compressive state. sig = (-100, -100, -100, 0, 0, 0).
    //      sig_star = (-101, -101, -101, 0, 0, 0). pmean = +101 > 0.25.
    //      Principal stresses of sig_star = (-101, -101, -101); max = -101;
    //      min_principal = -(-101) = 101 > 0.25.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 0, "R1 healthy compressive error==0");
        if (r.failed == nullptr) ++g_pass;
        else { ++g_fail; std::fprintf(stderr, "[FAIL] R1 failed=%s, expected null\n", r.failed); }
        expect_near(r.pmean,         101.0, 1e-12, "R1 pmean==101");
        expect_near(r.min_principal, 101.0, 1e-12, "R1 min_principal==101");
    }

    // (R2) Tensile mean stress: sig = (+10, +10, +10, 0, 0, 0).
    //      sig_star = (+9, +9, +9, 0, 0, 0). pmean = -9 < 0.25 -> "pmean".
    {
        const double sig[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R2 tensile-mean error==1");
        expect_failed(r, "pmean",  "R2 failed==pmean");
        expect_near(r.pmean, -9.0, 1e-12, "R2 pmean==-9 (negative)");
    }

    // (R3) Compressive mean but one tensile principal.
    //      Pick sig = (-100, -100, +50, 0, 0, 0).
    //      sig_star = (-101, -101, +49, 0, 0, 0).
    //      pmean = -((-101) + (-101) + 49) / 3 = -(-153)/3 = 51 > 0.25 OK.
    //      Principal stresses of sig_star (diagonal): (-101, -101, +49).
    //      max(S) = 49; min_principal = -49 < 0.25 -> "tension".
    {
        const double sig[6] = {-100.0, -100.0, 50.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R3 mixed-principal error==1");
        expect_failed(r, "tension", "R3 failed==tension (pmean ok, principal bad)");
        expect_near(r.pmean,         51.0, 1e-12, "R3 pmean==51 (compressive)");
        expect_near(r.min_principal, -49.0, 1e-12, "R3 min_principal==-49");
    }

    // (R4) NaN in y. Healthy sig but corrupt y[6] (an IS slot).
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[6] = std::numeric_limits<double>::quiet_NaN();
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R4 NaN-in-y error==1");
        expect_failed(r, "nan", "R4 failed==nan");
    }

    // (R5) Class precedence: tensile mean stress AND NaN in y.
    //      Corruption class (nan) overrides stress-domain class
    //      (pmean) so integrate_step does not silently absorb a
    //      numerically corrupt state as elastic. Codex 2026-05-02
    //      P1 follow-up: previous semantics ("first-failure-wins
    //      across classes") had a leak where a derived overflow
    //      labelled as "tension" first could mask "nan" later.
    {
        const double sig[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[7] = std::numeric_limits<double>::quiet_NaN();
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R5 multi-failure error==1");
        expect_failed(r, "nan",
                      "R5 corruption-class overrides stress-domain: failed==nan (not pmean)");
    }

    // (R6a) Fortran umatisnan_h also rejects +Inf, even in an IS slot
    //       that does not affect pmean/principal stress.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[6] = std::numeric_limits<double>::infinity();
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R6a +Inf-in-q error==1");
        expect_failed(r, "nan", "R6a failed==nan");
    }

    // (R6b) -Inf in a stress slot must also be rejected by the gate.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[0] = -std::numeric_limits<double>::infinity();
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R6b -Inf-in-stress error==1");
        expect_failed(r, "nan", "R6b failed==nan");
    }

    // (R6c) Finite overflow sentinels beyond +/-1.d30 mirror umatisnan_h.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[8] = 1.0e31;
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R6c finite-1e31 error==1");
        expect_failed(r, "nan", "R6c failed==nan");
    }

    // (R7) Fortran uses `.gt. 1.d30`, so exactly +1.d30 is not rejected.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[8] = 1.0e30;
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 0, "R7 finite-1e30 accepted");
    }

    // (R8) Sanity boundary: state right at the threshold.
    //      Use sig such that sig_star has pmean exactly equal to minstress.
    //      pmean = -tr(sig_star)/3 = minstress = p_t/4 = 0.25.
    //      tr(sig_star) = -3*0.25 = -0.75. With sig_star = (s,s,s,0,0,0),
    //      s = -0.25. So sig = (-0.25 + p_t, ...) = (0.75, 0.75, 0.75).
    //      Fortran uses `.le.` so EQUAL also fails (error=1).
    {
        const double sig[6] = {0.75, 0.75, 0.75, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R8 boundary pmean==minstress fails (Fortran .le.)");
        expect_near(r.pmean, minstress, 1e-14, "R8 pmean exactly at minstress");
        expect_failed(r, "pmean", "R8 boundary failed==pmean");
    }

    // (R9) Void-ratio admissibility: y[12] == 0 (finite, not caught by
    //      the umatisnan_h scan). The new void-rule must reject with
    //      failed = "void". Codex review 2026-05-01 (medium).
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[12] = 0.0;
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R9 void=0 -> error=1");
        expect_failed(r, "void", "R9 failed=='void'");
    }

    // (R10) Void-ratio admissibility: y[12] = -0.5 (finite, negative).
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[12] = -0.5;
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R10 void<0 -> error=1");
        expect_failed(r, "void", "R10 failed=='void'");
    }

    // (R11) First-failure-wins precedence: NaN void should still be
    //       reported as "nan" (caught by the umatisnan_h scan first),
    //       NOT "void". The void rule only fires for finite non-positive.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig, y);
        y[12] = std::numeric_limits<double>::quiet_NaN();
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R11 NaN-void error == 1");
        expect_failed(r, "nan",
                      "R11 NaN-void failed=='nan' (umatisnan_h takes precedence over void rule)");
    }

    // (R12) Class precedence: pmean fails AND void = 0. void is in
    //       the corruption class, so it overrides the stress-domain
    //       label "pmean". Codex 2026-05-02 P1 follow-up.
    {
        const double sig_tensile[6] = {10.0, 10.0, 10.0, 0.0, 0.0, 0.0};
        double y[kYDim]; build_y(sig_tensile, y);
        y[12] = 0.0;
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R12 pmean+void multi-fail error == 1");
        expect_failed(r, "void",
                      "R12 corruption-class overrides stress-domain: failed=='void' (not 'pmean')");
    }

    // (R13) Derived-overflow attack: sig has zero normals + 1e30
    //       shears, each component on the strict `> 1e30` boundary.
    //       The principal_stresses_3 derivation produces eigenvalues
    //       at ~2e30 magnitude (matrix has structure ~1e30*(J-I)
    //       where J is the all-ones matrix), so min_principal trips
    //       BOTH the "tension" threshold AND the umatisnan_h
    //       overflow sentinel.
    //
    //       Pre-2026-05-02 (first-failure-wins): "tension" set
    //       first, then "nan" preserved-but-masked -> failed ends
    //       up as "tension" -> integrate_step routes to elastic
    //       fallback -> overflow-finite sig accepted as success.
    //
    //       Post-2026-05-02 (class precedence): "nan" overrides
    //       "tension" -> integrate_step fails closed with error=3.
    //       This test locks the override at the check_RKF layer;
    //       the integrate_step-side regression is in
    //       test_integrate_step case 14.
    //
    //       Codex review 2026-05-02 (P1 high finding: derived
    //       overflow can be mistaken for elastic fallback).
    {
        const double sig_attack[6] = {0.0, 0.0, 0.0, 1.0e30, 1.0e30, 1.0e30};
        double y[kYDim]; build_y(sig_attack, y);
        const CheckRkfResult r = check_RKF(y, parms);
        expect_eq_int(r.error, 1, "R13 derived-overflow attack error == 1");
        expect_failed(r, "nan",
                      "R13 derived-overflow: failed=='nan' (overrides any 'tension')");
        // Sanity: this attack actually does trip min_principal beyond
        // 1e30. If a future change to principal_stresses_3 keeps
        // min_principal in range for this input, this test loses its
        // teeth -- the assertion below catches that drift.
        if (std::abs(r.min_principal) > 1.0e30) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] R13 sanity: |min_principal|=%g should exceed 1e30 "
                         "for this attack input (test no longer exercises the override)\n",
                         std::abs(r.min_principal));
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

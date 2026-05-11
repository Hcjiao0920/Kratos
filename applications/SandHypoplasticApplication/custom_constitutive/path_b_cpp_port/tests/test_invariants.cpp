// Standalone unit-test exerciser for the Path B kernel invariant
// helpers `inv_eps` and `inv_sig`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// Coverage:
//   inv_eps:
//     - Hydrostatic strain (eps_s == 0 -> sin3t = -1 early-return).
//     - Zero strain (degenerate hydrostatic).
//     - Axisymmetric extension and compression (sin3t = +/-1).
//     - Pure engineering shear (sin3t = 0).
//
//   inv_sig:
//     - Hydrostatic stress (qq = 0; norm2eta == 0 -> cos3t = -1).
//     - Triaxial axisymmetric compression (cos3t = -1).
//     - I1 == 0 axisymmetric deviatoric (triggers the tiny-divide
//       branch; cos3t numerically clamps to -1).
//     - I1 == 0 pure shear (tiny-divide branch with off-diagonal-only
//       Lode-angle contribution).
//     - Zero stress (sentinel: I1 == 0 with norm2eta > tiny driven
//       solely by the -1/3 normalization, exercising the cos3t clamp).
//
// Hand-computed expected values where exact; tight tolerances
// (~1e-13) for cases involving sqrt and division.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>

using Kratos::SandHypoCpp::StrainInvariants;
using Kratos::SandHypoCpp::StressInvariants;
using Kratos::SandHypoCpp::inv_eps;
using Kratos::SandHypoCpp::inv_sig;

namespace {

int g_pass = 0;
int g_fail = 0;

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
    // inv_eps: 5 cases
    // ====================================================================

    // (1) Hydrostatic strain: eps = a*[1,1,1,0,0,0].
    //     `a` MUST be binary-exact (a power of two times an integer)
    //     so that eps_v/3 returns to `a` without roundoff and edev is
    //     identically zero. The Fortran early return uses strict FP
    //     equality `eps_s.eq.zero`, so a non-exact `a` (e.g., 0.01)
    //     leaves roundoff in edev, eps_s ~ 2.5e-18, the early return
    //     does NOT fire, and sin3t falls through to the formula and
    //     clamps to a roundoff-determined sign. Our port mirrors this
    //     verbatim; this test pins down the early-return path itself.
    //     Expected: eps_v = 3a, edev = 0, eps_s = 0, sin3t = -1.
    {
        const double a = 0.5;
        const double eps[6] = {a, a, a, 0.0, 0.0, 0.0};
        const StrainInvariants r = inv_eps(eps);
        expect_near(r.eps_v, 3.0 * a, 0.0,         "eps hydrostatic eps_v=3a");
        expect_near(r.eps_s, 0.0,     0.0,         "eps hydrostatic eps_s=0 (exact)");
        expect_near(r.sin3t, -1.0,    0.0,         "eps hydrostatic sin3t=-1 (early return)");
    }

    // (2) Zero strain: degenerate hydrostatic.
    {
        const double eps[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const StrainInvariants r = inv_eps(eps);
        expect_near(r.eps_v, 0.0,  0.0, "eps zero eps_v=0");
        expect_near(r.eps_s, 0.0,  0.0, "eps zero eps_s=0");
        expect_near(r.sin3t, -1.0, 0.0, "eps zero sin3t=-1");
    }

    // (3) Axisymmetric extension: eps = [1, -0.5, -0.5, 0, 0, 0].
    //     Hand-computed: eps_v = 0; edev = eps; norm2 = 1.5;
    //                    eps_s = sqrt(2/3 * 1.5) = 1; sin3t = +1.
    {
        const double eps[6] = {1.0, -0.5, -0.5, 0.0, 0.0, 0.0};
        const StrainInvariants r = inv_eps(eps);
        expect_near(r.eps_v, 0.0, 0.0,    "eps axi-ext eps_v=0");
        expect_near(r.eps_s, 1.0, 1e-15,  "eps axi-ext eps_s=1");
        expect_near(r.sin3t, 1.0, 1e-13,  "eps axi-ext sin3t=+1");
    }

    // (4) Axisymmetric compression: eps = [-1, 0.5, 0.5, 0, 0, 0].
    //     Mirror image of (3); sin3t = -1.
    {
        const double eps[6] = {-1.0, 0.5, 0.5, 0.0, 0.0, 0.0};
        const StrainInvariants r = inv_eps(eps);
        expect_near(r.eps_v, 0.0,  0.0,    "eps axi-comp eps_v=0");
        expect_near(r.eps_s, 1.0,  1e-15,  "eps axi-comp eps_s=1");
        expect_near(r.sin3t, -1.0, 1e-13,  "eps axi-comp sin3t=-1");
    }

    // (5) Pure engineering shear: eps = [0,0,0,gamma,0,0] with gamma=1.
    //     Hand-computed: eps_v = 0; edev = [0,0,0,0.5,0,0];
    //                    norm2 = 2 * 0.25 = 0.5; eps_s = sqrt(2/3 * 0.5) = 1/sqrt(3);
    //                    edev2(i) = 0 for i != 1,2; tredev3 = 0; sin3t = 0.
    {
        const double eps[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        const StrainInvariants r = inv_eps(eps);
        expect_near(r.eps_v, 0.0,                      0.0,    "eps shear eps_v=0");
        expect_near(r.eps_s, 1.0 / std::sqrt(3.0),     1e-15,  "eps shear eps_s=1/sqrt(3)");
        expect_near(r.sin3t, 0.0,                      1e-15,  "eps shear sin3t=0");
    }

    // ====================================================================
    // inv_sig: 5 cases
    // ====================================================================

    // (1) Hydrostatic compression: sig = [-100,-100,-100,0,0,0].
    //     I1 = -300, pp = -100, sdev = 0, qq = 0;
    //     eta = sig/I1 = (1/3)(1,1,1,0,0,0); eta_d = 0; norm2eta = 0
    //     -> cos3t = -1 (early-return);
    //     I2 = 0.5*(sig:sig - I1^2) = 0.5*(30000 - 90000) = -30000  [Fortran sign];
    //     I3 = (-100)^3 = -1e6.
    {
        const double sig[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
        const StressInvariants r = inv_sig(sig);
        expect_near(r.I1,    -300.0,    0.0,    "sig hydro I1");
        expect_near(r.pp,    -100.0,    1e-15,  "sig hydro pp=-100");
        expect_near(r.qq,    0.0,       0.0,    "sig hydro qq=0");
        expect_near(r.cos3t, -1.0,      0.0,    "sig hydro cos3t=-1 (norm2eta<tiny early return)");
        expect_near(r.I2,    -30000.0,  1e-9,   "sig hydro I2=-30000 [Fortran sign 0.5*(sig:sig-I1^2)]");
        expect_near(r.I3,    -1.0e6,    1e-6,   "sig hydro I3=-1e6");
    }

    // (2) Triaxial axisymmetric compression: sig = [-100,-50,-50,0,0,0].
    //     I1 = -200, pp = -200/3;
    //     sdev = [-100/3, 50/3, 50/3, 0,0,0]; norm2 = 15000/9;
    //     qq = sqrt(3/2 * 15000/9) = 50;
    //     I2 = 0.5*(sig:sig - I1^2) = 0.5*(15000 - 40000) = -12500  [Fortran sign];
    //     I3 = (-100)*(-50)*(-50) = -250000;
    //     eta_d = [1/6, -1/12, -1/12, 0,0,0]; norm2eta = 1/24;
    //     tretadev3 = 1/288; cos3t = -1 (computed exactly).
    {
        const double sig[6] = {-100.0, -50.0, -50.0, 0.0, 0.0, 0.0};
        const StressInvariants r = inv_sig(sig);
        expect_near(r.I1,    -200.0,           0.0,    "sig tri-comp I1=-200");
        expect_near(r.pp,    -200.0 / 3.0,     1e-13,  "sig tri-comp pp=-200/3");
        expect_near(r.qq,    50.0,             1e-12,  "sig tri-comp qq=50");
        expect_near(r.cos3t, -1.0,             1e-13,  "sig tri-comp cos3t=-1");
        expect_near(r.I2,    -12500.0,         1e-10,  "sig tri-comp I2=-12500 [Fortran sign 0.5*(sig:sig-I1^2)]");
        expect_near(r.I3,    -250000.0,        1e-9,   "sig tri-comp I3=-250000");

        // Re-derive I2 directly from inputs (sig:sig with the
        // Voigt-stress weight of 2 on off-diagonals, then I1^2). If a
        // future change "fixes" the sign in inv_sig to the textbook
        // convention, this assertion fails alongside the numeric
        // expectation above, forcing a deliberate cross-cutting update
        // that also has to address Fortran and Path A.
        const double sigsig = sig[0]*sig[0] + sig[1]*sig[1] + sig[2]*sig[2]
                            + 2.0 * (sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5]);
        const double i2_from_formula = 0.5 * (sigsig - r.I1 * r.I1);
        expect_near(r.I2, i2_from_formula, 1e-10,
                    "sig tri-comp I2 == 0.5*(sig:sig-I1^2) recomputed from inputs [Fortran sign]");
    }

    // (3) I1 == 0 axisymmetric deviatoric: sig = [2,-1,-1,0,0,0].
    //     I1 = 0 (TRIGGERS tiny-divide branch); pp = 0;
    //     sdev = sig; norm2 = 6; qq = sqrt(3/2 * 6) = 3;
    //     I2 = 0.5*(sig:sig - I1^2) = 0.5*(6 - 0) = 3  [Fortran sign];
    //     I3 = 2*(1-0) = 2;
    //     eta = sig / 1e-18 (huge); eta_d ~ eta;
    //     norm2eta = 6 * 1e36; tretadev3 = 6 * 1e54;
    //     cos3t = -sqrt6 * 6e54 / (6e36)^1.5 = -1 algebraically.
    {
        const double sig[6] = {2.0, -1.0, -1.0, 0.0, 0.0, 0.0};
        const StressInvariants r = inv_sig(sig);
        expect_near(r.I1,    0.0,    0.0,    "sig I1=0 dev I1=0");
        expect_near(r.pp,    0.0,    0.0,    "sig I1=0 dev pp=0");
        expect_near(r.qq,    3.0,    1e-14,  "sig I1=0 dev qq=3");
        expect_near(r.cos3t, -1.0,   1e-13,  "sig I1=0 dev cos3t=-1 (tiny-divide branch)");
        expect_near(r.I2,    3.0,    0.0,    "sig I1=0 dev I2=3 [Fortran sign 0.5*(sig:sig-I1^2)]");
        expect_near(r.I3,    2.0,    0.0,    "sig I1=0 dev I3=2");
    }

    // (4) I1 == 0 pure shear: sig = [0,0,0,1,0,0].
    //     Exercises the tiny-divide branch with the off-diagonal-only
    //     Lode-angle path. Hand-computed result for cos3t in this
    //     branch tends toward zero as tiny -> 0; we assert this is
    //     bounded (|cos3t| <= 1) rather than pinning a roundoff-
    //     dominated value, since the I1==0 branch is a numerical
    //     placeholder per the header docstring.
    //     Other invariants: I1 = 0, pp = 0, norm2 = 2, qq = sqrt(3),
    //     I2 = 0.5*(sig:sig - I1^2) = 0.5*(2 - 0) = 1  [Fortran sign]; I3 = 0.
    {
        const double sig[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        const StressInvariants r = inv_sig(sig);
        expect_near(r.I1,    0.0,            0.0,    "sig I1=0 shear I1=0");
        expect_near(r.pp,    0.0,            0.0,    "sig I1=0 shear pp=0");
        expect_near(r.qq,    std::sqrt(3.0), 1e-15,  "sig I1=0 shear qq=sqrt(3)");
        expect_near(r.I2,    1.0,            0.0,    "sig I1=0 shear I2=1 [Fortran sign 0.5*(sig:sig-I1^2)]");
        expect_near(r.I3,    0.0,            0.0,    "sig I1=0 shear I3=0");
        // cos3t in this branch is dominated by the tiny-divide rescale;
        // assert only that it is finite and within the clamp bounds.
        const bool cos3t_bounded = std::isfinite(r.cos3t)
                                   && r.cos3t >= -1.0 - 1e-15
                                   && r.cos3t <=  1.0 + 1e-15;
        if (cos3t_bounded) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] sig I1=0 shear cos3t bounded: got %.17g\n",
                         r.cos3t);
        }
    }

    // (5) Zero stress: sig = 0. Exercises the cos3t clamp branch.
    //     I1 = 0 -> tiny-divide; eta = 0; eta_d = -(1/3)*(1,1,1,0,0,0);
    //     norm2eta = 3*(1/3)^2 = 1/3 (NOT < tiny);
    //     tretadev3 = 3*(-1/3)*(1/9) + 0 = -1/9;
    //     numer = -sqrt6 * (-1/9) = sqrt6/9;
    //     denom = (1/3)^1.5 = 1/(3*sqrt(3));
    //     cos3t pre-clamp = sqrt(6)/9 * 3*sqrt(3) = sqrt(18)/3 = sqrt(2) ~= 1.414;
    //     clamp triggers -> cos3t = +1 exactly.
    //     Other invariants: pp = 0, qq = 0, I2 = 0, I3 = 0.
    {
        const double sig[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const StressInvariants r = inv_sig(sig);
        expect_near(r.I1,    0.0, 0.0, "sig zero I1=0");
        expect_near(r.pp,    0.0, 0.0, "sig zero pp=0");
        expect_near(r.qq,    0.0, 0.0, "sig zero qq=0");
        expect_near(r.I2,    0.0, 0.0, "sig zero I2=0");
        expect_near(r.I3,    0.0, 0.0, "sig zero I3=0");
        // Clamp branch must produce exactly +1 (the clamp is
        // sign-preserving and computed via x/|x| so post-clamp value
        // is bit-exact +1 regardless of pre-clamp roundoff).
        expect_near(r.cos3t, 1.0, 0.0, "sig zero cos3t clamps to +1");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

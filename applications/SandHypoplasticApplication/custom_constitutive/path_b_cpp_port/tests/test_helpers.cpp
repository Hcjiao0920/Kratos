// Standalone unit-test exerciser for the Path B kernel linear-algebra
// helpers `dot_vect` and `matmul`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// Coverage:
//   dot_vect: n=3 plain; n=6 with shear-only, mixed, diagonal-only,
//             and asymmetric inputs; all three kinds; n=0 and n=1
//             edge cases. Hand-computed exact values, tol=0.
//   matmul:   2x3 * 3x2 (non-square layout sentinel); 3x3 identity
//             both sides; 1x3 * 3x1 scalar; 6x6 identity (kernel
//             usage size); 6x6 permutation * 6x1 to verify row/col
//             indexing without leaning on diagonal symmetry.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>

using Kratos::SandHypoCpp::DotProductKind;
using Kratos::SandHypoCpp::dot_vect;
using Kratos::SandHypoCpp::matmul;

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
    // dot_vect
    // ====================================================================

    // Plain, n=3.  [1,2,3] . [4,5,6] = 4 + 10 + 18 = 32.
    {
        const double a[3] = {1.0, 2.0, 3.0};
        const double b[3] = {4.0, 5.0, 6.0};
        expect_near(dot_vect(DotProductKind::Plain, a, b, 3), 32.0, 0.0,
                    "Plain n=3 [1,2,3].[4,5,6]=32");
    }

    // Plain, n=6 mixed: a=b=[1,2,3,4,5,6]. Self-dot = 1+4+9+16+25+36 = 91.
    {
        const double v[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        expect_near(dot_vect(DotProductKind::Plain, v, v, 6), 91.0, 0.0,
                    "Plain n=6 [1..6] self-dot = 91 (no weighting)");
    }

    // Shear-only Voigt vector exercises only off-diagonal slot (index 3):
    //   v = [0,0,0,1,0,0]
    // StressLike: weight 2 on index 3 -> 2.
    // StrainLike: weight 0.5 -> 0.5.
    // Plain:      weight 1   -> 1.
    {
        const double s[6] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        expect_near(dot_vect(DotProductKind::StressLike, s, s, 6), 2.0, 0.0,
                    "StressLike shear-12 self-dot = 2");
        expect_near(dot_vect(DotProductKind::StrainLike, s, s, 6), 0.5, 0.0,
                    "StrainLike shear-12 self-dot = 0.5");
        expect_near(dot_vect(DotProductKind::Plain, s, s, 6), 1.0, 0.0,
                    "Plain shear-12 self-dot = 1");
    }

    // Mixed Voigt vector exercises diagonal AND off-diagonal at once:
    //   v = [1,1,1,2,3,4]
    //   diagonal contribution = 1+1+1 = 3
    //   off-diagonal raw      = 4+9+16 = 29
    // StressLike: 3 + 2*29   = 61
    // StrainLike: 3 + 0.5*29 = 17.5
    // Plain:      3 + 29     = 32
    {
        const double v[6] = {1.0, 1.0, 1.0, 2.0, 3.0, 4.0};
        expect_near(dot_vect(DotProductKind::StressLike, v, v, 6), 61.0, 0.0,
                    "StressLike mixed [1,1,1,2,3,4] self-dot = 61");
        expect_near(dot_vect(DotProductKind::StrainLike, v, v, 6), 17.5, 0.0,
                    "StrainLike mixed [1,1,1,2,3,4] self-dot = 17.5");
        expect_near(dot_vect(DotProductKind::Plain, v, v, 6), 32.0, 0.0,
                    "Plain mixed [1,1,1,2,3,4] self-dot = 32");
    }

    // Diagonal-only vector: kind has no effect.
    //   v = [1,2,3,0,0,0]; self-dot = 1+4+9 = 14.
    {
        const double v[6] = {1.0, 2.0, 3.0, 0.0, 0.0, 0.0};
        expect_near(dot_vect(DotProductKind::StressLike, v, v, 6), 14.0, 0.0,
                    "StressLike diagonal-only = 14");
        expect_near(dot_vect(DotProductKind::StrainLike, v, v, 6), 14.0, 0.0,
                    "StrainLike diagonal-only = 14");
        expect_near(dot_vect(DotProductKind::Plain, v, v, 6), 14.0, 0.0,
                    "Plain diagonal-only = 14");
    }

    // Asymmetric a vs b, off-diagonal only:
    //   a = [1,0,0,2,0,0], b = [0,0,0,3,0,0]
    //   raw a.b = 2*3 = 6
    // StressLike: 2*6 = 12; StrainLike: 0.5*6 = 3; Plain: 6.
    {
        const double a[6] = {1.0, 0.0, 0.0, 2.0, 0.0, 0.0};
        const double b[6] = {0.0, 0.0, 0.0, 3.0, 0.0, 0.0};
        expect_near(dot_vect(DotProductKind::StressLike, a, b, 6), 12.0, 0.0,
                    "StressLike a.b shear-cross = 12");
        expect_near(dot_vect(DotProductKind::StrainLike, a, b, 6), 3.0, 0.0,
                    "StrainLike a.b shear-cross = 3");
        expect_near(dot_vect(DotProductKind::Plain, a, b, 6), 6.0, 0.0,
                    "Plain a.b shear-cross = 6");
    }

    // n=0 -> 0 in all kinds.
    {
        double dummy = 0.0;  // not read since n=0
        expect_near(dot_vect(DotProductKind::Plain, &dummy, &dummy, 0), 0.0, 0.0,
                    "Plain n=0 = 0");
        expect_near(dot_vect(DotProductKind::StressLike, &dummy, &dummy, 0), 0.0, 0.0,
                    "StressLike n=0 = 0");
    }

    // n=1: only diagonal slot; kind has no effect.
    {
        const double a[1] = {3.0};
        const double b[1] = {5.0};
        expect_near(dot_vect(DotProductKind::StressLike, a, b, 1), 15.0, 0.0,
                    "StressLike n=1 = 15");
        expect_near(dot_vect(DotProductKind::Plain, a, b, 1), 15.0, 0.0,
                    "Plain n=1 = 15");
    }

    // n=4 isolates exactly one off-diagonal slot (index 3) on top of three
    // diagonal slots; kind weighting must apply only to that one.
    //   v = [1,1,1,2]; self-dot:
    //   StressLike: 3 + 2*4 = 11; StrainLike: 3 + 0.5*4 = 5; Plain: 3+4=7.
    {
        const double v[4] = {1.0, 1.0, 1.0, 2.0};
        expect_near(dot_vect(DotProductKind::StressLike, v, v, 4), 11.0, 0.0,
                    "StressLike n=4 mixed = 11");
        expect_near(dot_vect(DotProductKind::StrainLike, v, v, 4), 5.0, 0.0,
                    "StrainLike n=4 mixed = 5");
        expect_near(dot_vect(DotProductKind::Plain, v, v, 4), 7.0, 0.0,
                    "Plain n=4 mixed = 7");
    }

    // ====================================================================
    // matmul
    // ====================================================================

    // 2x3 * 3x2 (non-square; would catch a row/col-major swap).
    //   A = [[1,2,3], [4,5,6]]      flat: {1,2,3,4,5,6}
    //   B = [[1,0], [0,1], [1,1]]   flat: {1,0,0,1,1,1}
    //   AB = [[1+0+3, 0+2+3], [4+0+6, 0+5+6]]
    //      = [[4, 5], [10, 11]]     flat: {4,5,10,11}
    {
        const double a[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        const double b[6] = {1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        double c[4] = {-99.0, -99.0, -99.0, -99.0};
        matmul(a, b, c, 2, 3, 2);
        expect_near(c[0], 4.0,  0.0, "matmul 2x3*3x2 c[0]=4");
        expect_near(c[1], 5.0,  0.0, "matmul 2x3*3x2 c[1]=5");
        expect_near(c[2], 10.0, 0.0, "matmul 2x3*3x2 c[2]=10");
        expect_near(c[3], 11.0, 0.0, "matmul 2x3*3x2 c[3]=11");
    }

    // 3x3 identity * A = A.
    {
        const double I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        const double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        double C[9] = {0};
        matmul(I, A, C, 3, 3, 3);
        for (int i = 0; i < 9; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "I*A: c[%d]=A[%d]=%g",
                          i, i, A[i]);
            expect_near(C[i], A[i], 0.0, label);
        }
    }

    // 3x3 A * identity = A (right-side identity).
    {
        const double A[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
        const double I[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        double C[9] = {0};
        matmul(A, I, C, 3, 3, 3);
        for (int i = 0; i < 9; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label), "A*I: c[%d]=A[%d]=%g",
                          i, i, A[i]);
            expect_near(C[i], A[i], 0.0, label);
        }
    }

    // 1x3 row vector times 3x1 column vector = 1x1 scalar.
    //   A = [1, 2, 3], B = [4; 5; 6]. AB = 1*4 + 2*5 + 3*6 = 32.
    {
        const double A[3] = {1.0, 2.0, 3.0};
        const double B[3] = {4.0, 5.0, 6.0};
        double C[1] = {-99.0};
        matmul(A, B, C, 1, 3, 1);
        expect_near(C[0], 32.0, 0.0, "matmul 1x3*3x1 = 32");
    }

    // 6x6 identity * H = H (kernel usage size; sanity check on the
    // 6-dim diagonal pattern).
    {
        double I[36] = {0};
        for (int k = 0; k < 6; ++k) I[k*6 + k] = 1.0;
        double H[36];
        for (int i = 0; i < 36; ++i) H[i] = static_cast<double>(i + 1);
        double C[36] = {0};
        matmul(I, H, C, 6, 6, 6);
        for (int i = 0; i < 36; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "I_6*H: c[%d]=H[%d]=%g", i, i, H[i]);
            expect_near(C[i], H[i], 0.0, label);
        }
    }

    // 6x6 permutation matrix times 6x1 vector. Catches a transposed
    // row/col formula where I*x would still pass.
    //
    // Permutation P maps input -> output as:
    //   row 0 -> input 1   (P[0,1]=1)
    //   row 1 -> input 2
    //   row 2 -> input 0
    //   row 3 -> input 4
    //   row 4 -> input 5
    //   row 5 -> input 3
    // So P * [1,2,3,4,5,6]^T = [2,3,1,5,6,4]^T.
    {
        double P[36] = {0};
        P[0*6 + 1] = 1.0;
        P[1*6 + 2] = 1.0;
        P[2*6 + 0] = 1.0;
        P[3*6 + 4] = 1.0;
        P[4*6 + 5] = 1.0;
        P[5*6 + 3] = 1.0;
        const double v[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
        double out[6] = {0};
        matmul(P, v, out, 6, 6, 1);
        const double expected[6] = {2.0, 3.0, 1.0, 5.0, 6.0, 4.0};
        for (int i = 0; i < 6; ++i) {
            char label[32];
            std::snprintf(label, sizeof(label), "P*v[%d]=%g", i, expected[i]);
            expect_near(out[i], expected[i], 0.0, label);
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

// Standalone unit-test exerciser for the Path B kernel `get_tan`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// Coverage (4 cases, two per istrain branch):
//
//   istrain = 0 (no-IS branch; selected by callers when m_R <= 0.5 in
//   the Fortran umat at .for:602-606):
//     (1) Hydrostatic compression with zero IS state. Verifies M is
//         left zero, L is symmetric and axially symmetric under
//         hydrostatic stress, N is axially symmetric, and H carries
//         only the void-evolution row.
//
//   istrain = 1 (IS branch; selected when m_R > 0.5):
//     (2) Hydrostatic compression with zero IS state. norm_del = 0
//         drives rho=0 and load=0, so the load<=0 sub-branch fires
//         and AA = 0; with textbook m_T = m_R = 1 this gives
//         M == m_R * L (== L). H_del comes from the load<=0 else
//         branch, i.e. identity; H_e is the void evolution row.
//     (3) Nonzero IS state aligned with deps -> load > 0 sub-branch.
//         Verifies H_del[0][0] = 1 - rho^beta_r against a hand value;
//         with m_T = m_R = 1 the off-axis elements of M still equal
//         the corresponding L elements (AA[i][j] = mm_temp3*N*eta_delta
//         is rank-1 only in column 0).
//     (4) Nonzero IS state ANTI-aligned with deps -> load < 0
//         sub-branch. H_del falls back to identity; with m_T = m_R = 1
//         AA = 0 again so M == L.
//
// In addition every case writes its full M / L / N / H arrays to a CSV
// in tests/path_b_results/ for later A-vs-B cross-validation.
//
// Known coverage gap (P3, deferred -- Codex 2026-05-02 P2 finding):
// per-entry numerical assertions against Fortran / Path A golden
// fixtures are not yet present. Structural assertions (symmetry,
// hand values for individual entries, post-fs ratio sanity, the
// rho>1 clamp / fb_temp1<0 fatal / AA[*][0] non-degenerate cases
// added below) catch most mistranslations but NOT a sign error or
// coefficient typo that preserves structure. See
// `path_b_cpp_port/README.md` "Known coverage gaps" for the
// planned harness design (out-of-process Fortran/Path A driver
// emits canonical M/L/N/H values to a data file; this suite would
// then read and assert per-entry agreement). Path B's isolation
// rule forbids reading Path A code directly, so the fixture file
// itself is the reference, not a runtime cross-call.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <string>

using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::TangentResult;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::get_tan;
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

void dump_tangent(std::ostream& out, const char* label,
                  const TangentResult& r) {
    out << "# " << label << "\n";
    out << "# error = " << r.error << "\n";
    out << "M = [\n";
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            out << r.M[i][j];
            out << (j < 5 ? "\t" : "\n");
        }
    }
    out << "]\nL = [\n";
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            out << r.L[i][j];
            out << (j < 5 ? "\t" : "\n");
        }
    }
    out << "]\nN = [";
    for (int i = 0; i < 6; ++i) {
        out << r.N[i] << (i < 5 ? "\t" : "");
    }
    out << "]\nH = [\n";
    for (std::size_t i = 0; i < kStateDim; ++i) {
        for (int j = 0; j < 6; ++j) {
            out << r.H[i][j];
            out << (j < 5 ? "\t" : "\n");
        }
    }
    out << "]\n\n";
}

}  // namespace

int main(int argc, char** argv) {
    const std::string outdir = (argc > 1) ? argv[1]
                                          : std::string("..\\..\\..\\tests\\path_b_results");
    const std::string outpath = outdir + "\\get_tan_textbook.txt";
    std::ofstream out(outpath);
    if (out) {
        out << std::scientific << std::setprecision(17);
        out << "# Path B kernel get_tan output for textbook params\n";
        out << "# Stress sig in kPa; Voigt order is Abaqus [11,22,33,12,13,23].\n\n";
    }

    const Parms16 parms = textbook_parms();
    const double void_ratio = 0.825;
    const double sig_compression[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};
    const double deps_small[6]      = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};

    // ====================================================================
    // (1) istrain = 0, hydrostatic compression, zero IS state
    // ====================================================================
    {
        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const TangentResult r = get_tan(deps_small, sig_compression, q, parms, 0);
        if (out) dump_tangent(out, "case 1: istrain=0, hydrostatic, zero IS", r);

        expect_eq_int(r.error, kErrorOk, "case1 error == 0");

        // M is identically zero in the no-IS branch.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case1 M[%d][%d] == 0 (no-IS)", i, j);
                expect_near(r.M[i][j], 0.0, 0.0, label);
            }
        }

        // L symmetry.
        for (int i = 0; i < 6; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case1 L symmetric: L[%d][%d]=L[%d][%d]", i, j, j, i);
                expect_near(r.L[i][j], r.L[j][i], 1e-10, label);
            }
        }

        // Axial symmetry of L diagonal under hydrostatic stress.
        expect_near(r.L[0][0], r.L[1][1], 1e-10, "case1 L[0][0] == L[1][1]");
        expect_near(r.L[1][1], r.L[2][2], 1e-10, "case1 L[1][1] == L[2][2]");
        expect_near(r.L[3][3], r.L[4][4], 1e-10, "case1 L[3][3] == L[4][4]");
        expect_near(r.L[4][4], r.L[5][5], 1e-10, "case1 L[4][4] == L[5][5]");

        // Cross-coupling between normal slots is axially symmetric.
        expect_near(r.L[0][1], r.L[0][2], 1e-10, "case1 L[0][1] == L[0][2]");
        expect_near(r.L[0][1], r.L[1][2], 1e-10, "case1 L[0][1] == L[1][2]");

        // Hydrostatic stress -> no normal-shear coupling and no
        // shear-shear off-diagonal coupling.
        for (int i = 0; i < 3; ++i) {
            for (int j = 3; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case1 L[%d][%d] == 0 (hydro)", i, j);
                expect_near(r.L[i][j], 0.0, 1e-10, label);
            }
        }
        expect_near(r.L[3][4], 0.0, 1e-10, "case1 L[3][4] == 0 (hydro)");
        expect_near(r.L[3][5], 0.0, 1e-10, "case1 L[3][5] == 0 (hydro)");
        expect_near(r.L[4][5], 0.0, 1e-10, "case1 L[4][5] == 0 (hydro)");

        // N axial structure.
        expect_near(r.N[0], r.N[1], 1e-10, "case1 N[0] == N[1]");
        expect_near(r.N[1], r.N[2], 1e-10, "case1 N[1] == N[2]");
        expect_near(r.N[3], 0.0,    1e-10, "case1 N[3] == 0 (hydro)");
        expect_near(r.N[4], 0.0,    1e-10, "case1 N[4] == 0 (hydro)");
        expect_near(r.N[5], 0.0,    1e-10, "case1 N[5] == 0 (hydro)");

        // H structure: rows 0..5 zero, row 6 = (1+v, 1+v, 1+v, 0, 0, 0).
        for (std::size_t i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case1 H[%zu][%d] == 0 (no-IS)", i, j);
                expect_near(r.H[i][j], 0.0, 0.0, label);
            }
        }
        const double one_plus_v = 1.0 + void_ratio;
        expect_near(r.H[6][0], one_plus_v, 0.0, "case1 H[6][0] == 1+v");
        expect_near(r.H[6][1], one_plus_v, 0.0, "case1 H[6][1] == 1+v");
        expect_near(r.H[6][2], one_plus_v, 0.0, "case1 H[6][2] == 1+v");
        expect_near(r.H[6][3], 0.0,        0.0, "case1 H[6][3] == 0");
        expect_near(r.H[6][4], 0.0,        0.0, "case1 H[6][4] == 0");
        expect_near(r.H[6][5], 0.0,        0.0, "case1 H[6][5] == 0");
    }

    // ====================================================================
    // (2) istrain = 1, hydrostatic compression, zero IS state
    //     With m_T = m_R = 1 and rho = 0, the IS contribution AA
    //     vanishes and M reduces to m_R * L (and L is post-scaled by
    //     fs at the very end of the routine, so r.M == r.L exactly).
    // ====================================================================
    {
        const double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const TangentResult r = get_tan(deps_small, sig_compression, q, parms, 1);
        if (out) dump_tangent(out, "case 2: istrain=1, hydrostatic, zero IS", r);

        expect_eq_int(r.error, kErrorOk, "case2 error == 0");

        // m_R = 1.0 (textbook) so M[i][j] == L[i][j] within FP roundoff.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[80];
                std::snprintf(label, sizeof(label),
                              "case2 M[%d][%d] == m_R*L[%d][%d] (m_R=1, rho=0)",
                              i, j, i, j);
                expect_near(r.M[i][j], r.L[i][j], 1e-9, label);
            }
        }

        // H rows 0..5: identity from the load<=0 else branch.
        for (std::size_t i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double expected = (static_cast<int>(i) == j) ? 1.0 : 0.0;
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case2 H[%zu][%d] = identity[%zu][%d]", i, j, i, j);
                expect_near(r.H[i][j], expected, 0.0, label);
            }
        }

        // H row 6 same as case 1.
        const double one_plus_v = 1.0 + void_ratio;
        expect_near(r.H[6][0], one_plus_v, 0.0, "case2 H[6][0] == 1+v");
        expect_near(r.H[6][1], one_plus_v, 0.0, "case2 H[6][1] == 1+v");
        expect_near(r.H[6][2], one_plus_v, 0.0, "case2 H[6][2] == 1+v");
        expect_near(r.H[6][3], 0.0,        0.0, "case2 H[6][3] == 0");
        expect_near(r.H[6][4], 0.0,        0.0, "case2 H[6][4] == 0");
        expect_near(r.H[6][5], 0.0,        0.0, "case2 H[6][5] == 0");
    }

    // ====================================================================
    // (3) istrain = 1, nonzero IS state aligned with deps -> load > 0.
    //     del aligned with axis 1, deps aligned with axis 1, both
    //     positive. eta_del = eta_eps = (1,0,0,0,0,0), load = 1.
    //     rho = norm_del / r_uc = 1e-4 / 3.3e-4.
    //     H_del[0][0] = 1 - rho^beta_r is the only off-identity slot.
    // ====================================================================
    {
        const double del0 = 1.0e-4;
        const double q[7] = {del0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q, parms, 1);
        if (out) dump_tangent(out, "case 3: istrain=1, nonzero IS, load>0", r);

        expect_eq_int(r.error, kErrorOk, "case3 error == 0");

        // Hand-computed rho^beta_r.
        const double r_uc = parms[11];
        const double beta_r = parms[12];
        const double rho = del0 / r_uc;            // norm_del = del0 (axis-1 only)
        const double rho_beta_r = std::pow(rho, beta_r);
        const double H00_expected = 1.0 - rho_beta_r;

        expect_near(r.H[0][0], H00_expected, 1e-12,
                    "case3 H[0][0] == 1 - rho^beta_r (load>0 IS update)");

        // Off-diagonals of column 0 (rows 1..5): zero (eta_del[i]=0 for i>0).
        for (std::size_t i = 1; i < 6; ++i) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "case3 H[%zu][0] == 0", i);
            expect_near(r.H[i][0], 0.0, 0.0, label);
        }

        // Row 0, columns 1..5: zero (eta_delta[j]=0 for j>0).
        for (int j = 1; j < 6; ++j) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "case3 H[0][%d] == 0", j);
            expect_near(r.H[0][j], 0.0, 0.0, label);
        }

        // Lower-right 5x5 block: identity (IU diag - 0).
        for (std::size_t i = 1; i < 6; ++i) {
            for (int j = 1; j < 6; ++j) {
                const double expected = (static_cast<int>(i) == j) ? 1.0 : 0.0;
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case3 H[%zu][%d] (5x5 ident)", i, j);
                expect_near(r.H[i][j], expected, 0.0, label);
            }
        }

        // Row 6: H_e structure, same as case 1/2.
        const double one_plus_v = 1.0 + void_ratio;
        expect_near(r.H[6][0], one_plus_v, 0.0, "case3 H[6][0] == 1+v");
        expect_near(r.H[6][3], 0.0,        0.0, "case3 H[6][3] == 0");

        // M off-axis-0 columns: with m_T = m_R = 1, AA[i][j] = mm_temp3
        // * N[i] * eta_delta[j], and eta_delta[j] = 0 for j > 0, so
        // M[i][j] = mm_temp1*L_unscaled[i][j] = m_R * fs * L_unscaled
        // = m_R * post-scale L. With m_R = 1, M[i][j] == L[i][j].
        for (int i = 0; i < 6; ++i) {
            for (int j = 1; j < 6; ++j) {
                char label[80];
                std::snprintf(label, sizeof(label),
                              "case3 M[%d][%d] == L[%d][%d] (m_R=m_T=1, j>0 col)",
                              i, j, i, j);
                expect_near(r.M[i][j], r.L[i][j], 1e-8, label);
            }
        }
    }

    // ====================================================================
    // (4) istrain = 1, nonzero IS state ANTI-aligned with deps -> load < 0.
    //     del = +axis-1, deps = -axis-1. load = -1.
    //     H_del = identity (else branch). With m_T = m_R = 1, AA = 0
    //     and M = m_R * L = L.
    // ====================================================================
    {
        const double del0 = 1.0e-4;
        const double q[7] = {del0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {-1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q, parms, 1);
        if (out) dump_tangent(out, "case 4: istrain=1, nonzero IS, load<0", r);

        expect_eq_int(r.error, kErrorOk, "case4 error == 0");

        // H rows 0..5 = identity (else branch).
        for (std::size_t i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double expected = (static_cast<int>(i) == j) ? 1.0 : 0.0;
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case4 H[%zu][%d] = identity (load<0)", i, j);
                expect_near(r.H[i][j], expected, 0.0, label);
            }
        }

        // M == L (m_R = m_T = 1, rho^chi*(m_R-m_T) = 0).
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[80];
                std::snprintf(label, sizeof(label),
                              "case4 M[%d][%d] == L[%d][%d]", i, j, i, j);
                expect_near(r.M[i][j], r.L[i][j], 1e-9, label);
            }
        }
    }

    // ====================================================================
    // (5) Defensive guard: r_uc = 0 + istrain = 1 must return error = 10.
    //     Codex review 2026-05-01 (high) flagged that the original
    //     unconditional `rho = norm_del / r_uc` produces NaN under the
    //     common zero-IS state when r_uc = 0, and the NaN propagates
    //     into M/H while error stays 0. The kernel guards this in
    //     the IS branch as a deliberate deviation from the Fortran.
    //
    //     Note: after Codex 2026-05-02 (medium) hoisted "r_uc > 0 if
    //     m_R > 0.5" into check_parms, the (r_uc=0, m_R=1) combo is
    //     rejected at validation. We bypass check_parms here by
    //     mutating a validated Parms16 directly to exercise the
    //     in-kernel belt-and-suspenders guard. case 5b still uses
    //     check_parms because m_R=0.4 + r_uc=0 is in the legal regime.
    // ====================================================================
    {
        Parms16 parms_runc_zero = textbook_parms();
        parms_runc_zero[11] = 0.0;          // bypass check_parms cross-check

        const double q_zero_is[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const TangentResult r = get_tan(deps_small, sig_compression, q_zero_is,
                                        parms_runc_zero, 1);
        if (out) dump_tangent(out,
                              "case 5: r_uc=0 + istrain=1 (guard fires)", r);

        expect_eq_int(r.error, 10,
                      "case5 r_uc<=0 + istrain=1 -> error=10 (NaN guard fires)");

        // M/H/L/N must all be left at the zero-init value when the guard
        // fires (no spurious NaN in any output).
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[80];
                std::snprintf(label, sizeof(label),
                              "case5 M[%d][%d] == 0 (guard fires before any IS work)",
                              i, j);
                expect_near(r.M[i][j], 0.0, 0.0, label);
            }
        }
        // L is computed BEFORE the guard, so it may already be populated.
        // We only assert M / H stay zero (which is what matters for
        // downstream contamination); skipping L assertion here.

        // Sanity: with istrain == 0, the guard does NOT fire even with
        // r_uc = 0 (no-IS path doesn't divide by r_uc). Use a parms
        // set that check_parms accepts (m_R = 0.4 < 0.5) so this test
        // also documents the legal Fortran regime.
        Props16 props_low_mR_runc0 = textbook_props();
        props_low_mR_runc0[9]  = 0.4;
        props_low_mR_runc0[11] = 0.0;
        const auto check_lo = check_parms(props_low_mR_runc0);
        expect_eq_int(check_lo.error, kErrorOk,
                      "case5b m_R=0.4 + r_uc=0 still accepted by check_parms");
        const TangentResult r0 = get_tan(deps_small, sig_compression, q_zero_is,
                                         check_lo.parms, 0);
        expect_eq_int(r0.error, kErrorOk,
                      "case5b r_uc=0 + istrain=0 OK (guard scoped to IS branch)");
    }

    // ====================================================================
    // (6) m_R != m_T (5 vs 2), nonzero IS aligned with deps -> load > 0.
    //     Catches mistranslations of the rho^chi / m_T / m_R terms in
    //     mm_temp1 that case 3 cannot detect (because case 3 uses
    //     m_T = m_R = 1, which makes mm_temp2 = mm_temp4 = 0 and
    //     mm_temp1 = m_R * fs trivially).
    //
    //     For j > 0, eta_delta[j] = 0 collapses AA[i][j] to zero, so:
    //         M[i][j] = (rho^chi * m_T + (1 - rho^chi) * m_R) * L[i][j]
    //     The scalar ratio M[i][j] / L[i][j] is hand-computable and
    //     is NOT 1 here; it falls between m_T (= 2) and m_R (= 5).
    //     This pins the exact mm_temp1 coefficient against an
    //     independent recomputation.
    // ====================================================================
    {
        Props16 props_split_mR = textbook_props();
        props_split_mR[9]  = 5.0;  // m_R
        props_split_mR[10] = 2.0;  // m_T
        const auto check6 = check_parms(props_split_mR);
        expect_eq_int(check6.error, kErrorOk,
                      "case6 m_R=5 m_T=2 accepted");

        const double del0 = 1.0e-4;
        const double q[7] = {del0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q,
                                        check6.parms, 1);
        if (out) dump_tangent(out,
                              "case 6: m_R=5 m_T=2, nonzero IS, load>0", r);

        expect_eq_int(r.error, kErrorOk, "case6 error == 0");

        const double r_uc6   = check6.parms[11];
        const double chi6    = check6.parms[13];
        const double m_R6    = check6.parms[9];
        const double m_T6    = check6.parms[10];
        const double rho6     = del0 / r_uc6;
        const double rho_chi6 = std::pow(rho6, chi6);
        const double mm_temp1_over_fs = rho_chi6 * m_T6
                                        + (1.0 - rho_chi6) * m_R6;

        // Sanity: ratio should sit between m_T and m_R (rho_chi6 small,
        // ratio ~= m_R = 5 - small correction). Fail loudly if the
        // computation degenerates to 1 or 0 or m_T (which would mean
        // a future implementation changed the term ordering).
        const bool ratio_in_band = mm_temp1_over_fs > 4.5
                                && mm_temp1_over_fs < 5.0;
        if (ratio_in_band) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case6 mm_temp1/fs sanity: got %.10g, expected in (4.5, 5.0)\n",
                         mm_temp1_over_fs);
        }

        // For j > 0 (eta_delta[j] = 0), M[i][j] should equal
        // mm_temp1_over_fs * L[i][j] within FP arithmetic precision.
        // The two sides differ in operation ordering (LHS computed
        // pre-final-scaling, RHS uses post-scaled L), so use a
        // proportional tolerance.
        for (int i = 0; i < 6; ++i) {
            for (int j = 1; j < 6; ++j) {
                const double expected = mm_temp1_over_fs * r.L[i][j];
                const double tol = 1.0e-9 * (std::abs(expected) + 1.0);
                char label[96];
                std::snprintf(label, sizeof(label),
                              "case6 M[%d][%d] = (rho^chi*m_T + (1-rho^chi)*m_R) * L[%d][%d]",
                              i, j, i, j);
                expect_near(r.M[i][j], expected, tol, label);
            }
        }
    }

    // ====================================================================
    // (7) Triaxial axisymmetric stress (axis-1 unique), no-IS branch.
    //     Cases 1-4 all used HYDROSTATIC stress, where eta_dev = 0,
    //     tanpsi = 0, FF = 1 trivially. With non-hydrostatic stress
    //     these intermediates carry real values, exercising the
    //     tanpsi/cos3t branches of the FF computation that case 1
    //     did not reach.
    //
    //     For sig = (-100, -50, -50, 0, 0, 0) the (2 <-> 3) axis swap
    //     leaves the stress invariant. The tangent must respect this:
    //         L[1][1] == L[2][2]   (slot 22-22 vs 33-33)
    //         L[3][3] == L[4][4]   (slot 12-12 vs 13-13)
    //         L[0][1] == L[0][2]   (slot 11-22 vs 11-33)
    //     Slot 5 (= 23) is invariant under the swap, so L[5][5] is
    //     standalone. Slot 0 is also invariant.
    // ====================================================================
    {
        const double sig_triax[6] = {-100.0, -50.0, -50.0, 0.0, 0.0, 0.0};
        const double q_zero_is[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const TangentResult r = get_tan(deps_small, sig_triax, q_zero_is,
                                        parms, 0);
        if (out) dump_tangent(out, "case 7: triaxial axisym, no-IS", r);

        expect_eq_int(r.error, kErrorOk, "case7 error == 0");

        // Axis-2/3 swap symmetry of L.
        expect_near(r.L[1][1], r.L[2][2], 1e-9,
                    "case7 L[1][1]==L[2][2] (axis 2<->3 swap, normal-normal)");
        expect_near(r.L[3][3], r.L[4][4], 1e-9,
                    "case7 L[3][3]==L[4][4] (axis 2<->3 swap, slot 12<->13)");
        expect_near(r.L[0][1], r.L[0][2], 1e-9,
                    "case7 L[0][1]==L[0][2] (axis 2<->3 swap, normal-cross)");

        // L symmetric.
        for (int i = 0; i < 6; ++i) {
            for (int j = i + 1; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case7 L symmetric L[%d][%d]==L[%d][%d]", i, j, j, i);
                expect_near(r.L[i][j], r.L[j][i], 1e-9, label);
            }
        }

        // For axisymmetric stress, the (1,2) shear slot (slot 3) does
        // NOT decouple from the normal slots in general -- but slot 5
        // (component 23) does, because component (2,3) is the ONLY one
        // unchanged by both the normal-axis structure and the 2<->3
        // swap. Assert that L's slot-5 row/column off-diagonals to the
        // normal block are zero.
        for (int j = 0; j < 3; ++j) {
            char label[64];
            std::snprintf(label, sizeof(label),
                          "case7 L[5][%d]==0 (slot 23 decoupled from normals)", j);
            expect_near(r.L[5][j], 0.0, 1e-9, label);
        }

        // No-IS branch: M still zero.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char label[64];
                std::snprintf(label, sizeof(label),
                              "case7 M[%d][%d] == 0 (no-IS)", i, j);
                expect_near(r.M[i][j], 0.0, 0.0, label);
            }
        }
    }

    // ====================================================================
    // (8) Shear-slot IS evolution: del[3] = 1e-4 (component 12).
    //     Cases 3-4 only put IS on a normal slot (axis-1). This case
    //     exercises eta_delta[3..5], which the get_tan code halves
    //     from eta_del. A bug like writing `eta_delta[3] = eta_del[2]`
    //     instead of `0.5 * eta_del[3]` would slip through case 3 (the
    //     halved slots are zero there) but is caught here.
    //
    //     With del = (0,0,0,1e-4,0,0):
    //       norm_del2 = 0.5 * (1e-4)^2  (StrainLike halves slot 3)
    //       norm_del  = sqrt(0.5)*1e-4 = 7.07e-5
    //       eta_del[3] = 1e-4 / 7.07e-5 = sqrt(2)
    //       eta_delta[3] = sqrt(2)/2 = 1/sqrt(2)
    //
    //     With deps[3] = 1e-4 (aligned), load = 0.5 * sqrt(2) * sqrt(2) = 1 > 0.
    //     H_del[3][3] = 1 - rho^beta_r * sqrt(2) * (sqrt(2)/2) = 1 - rho^beta_r.
    //     All other H_del entries follow the identity-with-rank-1-update pattern.
    //
    //     rho = norm_del / r_uc = 7.07e-5 / 3.3e-4 ~= 0.2143.
    // ====================================================================
    {
        const double del3 = 1.0e-4;
        const double q[7] = {0.0, 0.0, 0.0, del3, 0.0, 0.0, void_ratio};
        const double deps[6] = {0.0, 0.0, 0.0, 1.0e-4, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q, parms, 1);
        if (out) dump_tangent(out, "case 8: shear-slot IS, load>0", r);

        expect_eq_int(r.error, kErrorOk, "case8 error == 0");

        const double r_uc8   = parms[11];
        const double beta_r8 = parms[12];
        // norm_del = sqrt(0.5) * del3 (StrainLike norm halves slot 3).
        const double norm_del8 = std::sqrt(0.5) * del3;
        const double rho8        = norm_del8 / r_uc8;
        const double rho_beta_r8 = std::pow(rho8, beta_r8);
        const double H33_expected = 1.0 - rho_beta_r8;

        expect_near(r.H[3][3], H33_expected, 1e-12,
                    "case8 H[3][3] == 1 - rho^beta_r at slot 3 (load>0 IS update)");

        // Identity-with-hole pattern: H[3][3] is the only non-IU slot.
        for (std::size_t i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                if (i == 3 && j == 3) continue;  // checked above
                const double expected = (static_cast<int>(i) == j) ? 1.0 : 0.0;
                char label[80];
                std::snprintf(label, sizeof(label),
                              "case8 H[%zu][%d] (identity outside slot 3,3)", i, j);
                expect_near(r.H[i][j], expected, 1e-15, label);
            }
        }

        // Row 6 same as before.
        const double one_plus_v = 1.0 + void_ratio;
        expect_near(r.H[6][0], one_plus_v, 0.0, "case8 H[6][0] == 1+v");
        expect_near(r.H[6][3], 0.0,        0.0, "case8 H[6][3] == 0");
    }

    // ====================================================================
    // (9) AA[*][0] non-degenerate (Codex 2026-05-02 P2 follow-up).
    //     Case 6 only checked j > 0 columns of M, where eta_delta[j] = 0
    //     collapses AA[i][j] to zero. The j == 0 column is where the
    //     IS branch's rank-1 update LIVES, so a bug like dropping AA
    //     entirely (M[i][j] = mm_temp1 * L[i][j] for ALL j) would slip
    //     through case 6.
    //
    //     Setup mirrors case 6: m_R=5, m_T=2, del=(1e-4,0,0,0,0,0),
    //     deps aligned with del so load > 0. eta_delta = (1,0,0,0,0,0).
    //
    //     Closed-form derivation for the j == 0 column (load > 0):
    //         M[i][0] = mm_temp1 * L_unscaled[i][0] + AA[i][0]
    //         AA[i][0] = mm_temp2 * Leta[i] + mm_temp3 * N_unscaled[i]
    //         (Leta = L_unscaled @ eta_del; with eta_del = (1,0,...),
    //          Leta[i] = L_unscaled[i][0])
    //     Substituting the post-fs scalings at .for:567-577:
    //         L_unscaled[i][0]   = r.L[i][0] / fs
    //         N_unscaled[i]      = r.N[i] / (fs * fd)
    //         mm_temp1 = (rho^chi*m_T + (1-rho^chi)*m_R) * fs
    //         mm_temp2 = rho^chi * (1 - m_T) * fs
    //         mm_temp3 = rho^chi * fs * fd
    //     gives:
    //         M[i][0] = [rho^chi + (1-rho^chi) * m_R] * r.L[i][0]
    //                 +  rho^chi                      * r.N[i]
    //
    //     Killer assertion: M[i][0] differs MEANINGFULLY from
    //     mm_temp1_over_fs * r.L[i][0] -- proves AA actually contributes
    //     at j == 0 and case 6's "j > 0 only" coverage isn't tautological.
    // ====================================================================
    {
        Props16 props_split_mR = textbook_props();
        props_split_mR[9]  = 5.0;  // m_R
        props_split_mR[10] = 2.0;  // m_T
        const auto check9 = check_parms(props_split_mR);
        expect_eq_int(check9.error, kErrorOk, "case9 m_R=5 m_T=2 accepted");

        const double del0    = 1.0e-4;
        const double q[7]    = {del0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q,
                                        check9.parms, 1);
        if (out) dump_tangent(out,
                              "case 9: AA non-degenerate (load>0, j=0 column)", r);

        expect_eq_int(r.error, kErrorOk, "case9 error == 0");

        const double r_uc9     = check9.parms[11];
        const double chi9      = check9.parms[13];
        const double m_R9      = check9.parms[9];
        const double m_T9      = check9.parms[10];
        const double rho9      = del0 / r_uc9;
        const double rho_chi9  = std::pow(rho9, chi9);
        const double mm_temp1_over_fs = rho_chi9 * m_T9
                                        + (1.0 - rho_chi9) * m_R9;
        const double coef_L_full = rho_chi9 + (1.0 - rho_chi9) * m_R9;
        const double coef_N_full = rho_chi9;

        // Direct verification: M[i][0] == coef_L_full * L[i][0] + coef_N_full * N[i].
        for (int i = 0; i < 6; ++i) {
            const double expected = coef_L_full * r.L[i][0]
                                  + coef_N_full * r.N[i];
            const double tol = 1.0e-9 * (std::abs(expected) + 1.0);
            char label[120];
            std::snprintf(label, sizeof(label),
                          "case9 M[%d][0] = (rho^chi + (1-rho^chi)*m_R)*L[%d][0] "
                          "+ rho^chi*N[%d]", i, i, i);
            expect_near(r.M[i][0], expected, tol, label);
        }

        // Killer (anti-tautology): SOMEWHERE on the j==0 column, the
        // full formula must differ from the mm_temp1_over_fs * L[i][0]
        // value (which is what the j > 0 columns get) by enough that
        // a bug dropping the AA term wouldn't pass FP noise.
        double max_rel_diff = 0.0;
        for (int i = 0; i < 6; ++i) {
            const double full_val   = coef_L_full * r.L[i][0]
                                    + coef_N_full * r.N[i];
            const double naive_val  = mm_temp1_over_fs * r.L[i][0];
            const double scale = std::abs(full_val) + std::abs(naive_val) + 1.0;
            const double rel   = std::abs(full_val - naive_val) / scale;
            if (rel > max_rel_diff) max_rel_diff = rel;
        }
        if (max_rel_diff > 1.0e-6) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case9 anti-tautology: AA contribution at j==0 "
                         "is FP-noise-level (max rel diff %.3g); a bug dropping "
                         "AA would not be caught here\n", max_rel_diff);
        }
    }

    // ====================================================================
    // (10) rho > 1 clamp (Codex 2026-05-02 P2 follow-up).
    //      The IS branch at .for:884 / kernel.cpp:487-488 clamps rho to
    //      1 when norm_del / r_uc > 1 (saturation: intergranular strain
    //      maxed out). With rho clamped to 1, rho^chi = 1 and:
    //          mm_temp1 = (1 * m_T + 0 * m_R) * fs = m_T * fs
    //      So for j > 0 columns (eta_delta[j] = 0, AA = 0):
    //          M[i][j] / L[i][j] = m_T   (NOT m_R)
    //
    //      Setup: m_R=5, m_T=2, huge del=(1.0,0,...). norm_del = 1.0,
    //      r_uc = textbook 3.3e-4, so rho = 1.0/3.3e-4 ~ 3030 >> 1
    //      and the clamp definitely fires.
    //
    //      Killer assertion: ratio M/L at j > 0 is m_T (= 2.0), NOT m_R
    //      (= 5.0). Distinct from case 6 which had rho ~0.3 and ratio
    //      close to m_R (~5).
    // ====================================================================
    {
        Props16 props_split_mR = textbook_props();
        props_split_mR[9]  = 5.0;   // m_R
        props_split_mR[10] = 2.0;   // m_T
        const auto check10 = check_parms(props_split_mR);
        expect_eq_int(check10.error, kErrorOk,
                      "case10 m_R=5 m_T=2 accepted");

        // del HUGE so norm_del / r_uc >> 1 -> rho clamped to 1.
        // physically extreme but legal; the clamp test is the point.
        const double del0    = 1.0;
        const double q[7]    = {del0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
        const TangentResult r = get_tan(deps, sig_compression, q,
                                        check10.parms, 1);
        if (out) dump_tangent(out, "case 10: rho>1 clamp", r);

        expect_eq_int(r.error, kErrorOk, "case10 error == 0");

        // After rho clamp, mm_temp1_over_fs = m_T = 2.0.
        const double m_T10               = check10.parms[10];
        const double mm_temp1_over_fs    = m_T10;  // (1*m_T + 0*m_R)
        const double m_R10               = check10.parms[9];

        // For j > 0 (eta_delta = 0, AA = 0): M[i][j] = m_T * L[i][j].
        // Use loose absolute tolerance because L can be zero on some
        // entries; relative tolerance with explicit floor.
        for (int i = 0; i < 6; ++i) {
            for (int j = 1; j < 6; ++j) {
                const double expected = mm_temp1_over_fs * r.L[i][j];
                const double tol = 1.0e-9 * (std::abs(expected) + 1.0);
                char label[120];
                std::snprintf(label, sizeof(label),
                              "case10 M[%d][%d] = m_T * L[%d][%d] (rho clamped to 1)",
                              i, j, i, j);
                expect_near(r.M[i][j], expected, tol, label);
            }
        }

        // Killer: ratio is m_T (=2), NOT m_R (=5). Pick an entry with
        // non-trivial L magnitude and verify it tracks m_T scaling.
        bool any_strong_check = false;
        for (int i = 0; i < 6; ++i) {
            for (int j = 1; j < 6; ++j) {
                if (std::abs(r.L[i][j]) > 1.0e3) {  // pick a meaningful entry
                    const double ratio = r.M[i][j] / r.L[i][j];
                    if (std::abs(ratio - m_T10) < 0.01
                     && std::abs(ratio - m_R10) > 1.0) {
                        any_strong_check = true;
                    }
                }
            }
        }
        if (any_strong_check) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case10 rho-clamp ratio: no strong entry showed "
                         "M/L ~= m_T (=2) and not m_R (=5). Either rho clamp "
                         "is broken OR test setup degenerated.\n");
        }
    }

    // ====================================================================
    // (11) fb_temp1 < 0 fatal (Codex 2026-05-02 P2 follow-up).
    //      Fortran .for:838-841: if `temp1 = 3 + a^2 - a*sqrt(3) *
    //      ((ei0-ed0)/(ec0-ed0))^alpha` is negative, the original code
    //      called `xit_h` ("factor fb not defined"). Path B translates
    //      this to `error = kErrorFatal`. The guard fires BEFORE L/N
    //      are computed (see kernel.cpp:423-428) so all output matrices
    //      stay zero.
    //
    //      Construct fb_temp1 < 0: with textbook phi=34deg, a~2.674,
    //      a^2~7.15, 3+a^2~10.15, a*sqrt(3)~4.63. Need
    //      pow((ei0-ed0)/(ec0-ed0), alpha) > 10.15/4.63 ~= 2.19.
    //      With alpha=0.24: ratio > 2.19^(1/0.24) = 2.19^4.17 ~= 27.
    //      Pick ed0=0.5, ec0=0.51, ei0=1.0 -> ratio = 0.5/0.01 = 50
    //      -> pow(50, 0.24) = exp(0.24 * ln(50)) ~= 2.557
    //      -> fb_temp1 ~= 10.15 - 4.63*2.557 ~= -1.69 < 0.
    //
    //      Verify the constructed parms still pass check_parms (no
    //      negative values, cross-parameter check OK), then the get_tan
    //      call must emit error = kErrorFatal with M/L/N/H all zero.
    // ====================================================================
    {
        Props16 props_fb_neg = textbook_props();
        props_fb_neg[4] = 0.5;    // ed0
        props_fb_neg[5] = 0.51;   // ec0
        props_fb_neg[6] = 1.0;    // ei0
        // alpha (props[7]) = textbook 0.24 keeps fb_temp1 negative.
        const auto check11 = check_parms(props_fb_neg);
        expect_eq_int(check11.error, kErrorOk,
                      "case11 fb_temp1<0 props pass check_parms (no negative)");

        const double q[7]    = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, void_ratio};
        const double deps[6] = {-1.0e-3, 0.0, 0.0, 0.0, 0.0, 0.0};

        // istrain == 0 first (no-IS path also computes fb).
        const TangentResult r0 = get_tan(deps, sig_compression, q,
                                         check11.parms, 0);
        expect_eq_int(r0.error, 10,
                      "case11 fb_temp1<0 + istrain=0 -> error=10 fatal");
        // Guard fires BEFORE L/N computation; all matrices stay zero.
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char lbl1[80], lbl2[80];
                std::snprintf(lbl1, sizeof(lbl1),
                              "case11 r0.L[%d][%d] == 0 (fb_temp1 fatal)", i, j);
                std::snprintf(lbl2, sizeof(lbl2),
                              "case11 r0.M[%d][%d] == 0 (fb_temp1 fatal)", i, j);
                expect_near(r0.L[i][j], 0.0, 0.0, lbl1);
                expect_near(r0.M[i][j], 0.0, 0.0, lbl2);
            }
            char lbl[80];
            std::snprintf(lbl, sizeof(lbl),
                          "case11 r0.N[%d] == 0 (fb_temp1 fatal)", i);
            expect_near(r0.N[i], 0.0, 0.0, lbl);
        }

        // istrain == 1 same outcome.
        const TangentResult r1 = get_tan(deps, sig_compression, q,
                                         check11.parms, 1);
        expect_eq_int(r1.error, 10,
                      "case11 fb_temp1<0 + istrain=1 -> error=10 fatal");
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                char lbl[96];
                std::snprintf(lbl, sizeof(lbl),
                              "case11 r1.L[%d][%d] == 0 (fb_temp1 before L)", i, j);
                expect_near(r1.L[i][j], 0.0, 0.0, lbl);
            }
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

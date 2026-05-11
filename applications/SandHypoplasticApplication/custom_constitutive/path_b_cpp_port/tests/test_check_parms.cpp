// Standalone unit-test exerciser for the Path B kernel `check_parms`.
//
// Compile + run via `build_and_run.bat` (sibling). The test:
//   1. Validates the textbook medium-dense sand parameter set
//      (project_umat_model.md). Expects accept; parms[0] == 34 deg in rad;
//      remaining 15 slots pass through bit-identical.
//   2. Locks every Fortran sentinel rejection rule in
//      `check_parms_h` (.for:389-508), including the `failed` parameter
//      name and `bad_value` returned to the caller:
//        phi <= 0                  -> reject "phi"
//        m_R, m_T, r_uc, beta_r,
//        chi, bulk_w, p_t  < 0     -> reject named slot
//        m_R == 0 (and == 0 for
//        the other "< 0" slots)    -> accepted (Fortran uses strict `<`)
//   3. Writes the validated `parms` for the textbook case to
//      `tests/path_b_results/check_parms_textbook.txt` for cross-path
//      comparison.
//
// Note on non-finite inputs (NaN, +/-Inf): `check_parms` mirrors the
// Fortran's IEEE-blind comparisons, so non-finite values currently pass
// validation. This is intentional for the literal-port phase; the host
// `Check()` layer will add a finite guard before the kernel ships to
// any caller. No tests added here yet.
//
// Cross-parameter constraint (defensive deviation from Fortran, Codex
// 2026-05-02 medium follow-up): when m_R > 0.5 (IS branch active),
// r_uc must be strictly positive. Fortran accepts `r_uc = 0` and lets
// `get_tan_h` NaN-poison via `rho = norm_del / r_uc`; we hoist the
// rejection here so the bad combo never reaches integration. The
// previous in-kernel guard at get_tan is kept as belt-and-suspenders.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <string>

using Kratos::SandHypoCpp::CheckParmsResult;
using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::kErrorFatal;
using Kratos::SandHypoCpp::kErrorOk;

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

// Asserts the success contract: error==Ok, failed==nullptr, bad_value==NaN.
void expect_accepted(const CheckParmsResult& r, const char* label) {
    const bool good = (r.error == kErrorOk)
                      && (r.failed == nullptr)
                      && std::isnan(r.bad_value);
    if (good) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr,
                     "[FAIL] %s: error=%d, failed=%s, bad_value=%.17g\n",
                     label, r.error,
                     r.failed ? r.failed : "(null)",
                     r.bad_value);
    }
}

// Asserts the rejection contract: error==Fatal, failed names the slot,
// bad_value equals the rejected numeric value (bit-exact: no tolerance).
void expect_rejected(const CheckParmsResult& r,
                     const char* expected_name,
                     double expected_bad_value,
                     const char* label) {
    expect_eq_int(r.error, kErrorFatal, label);

    const bool name_ok = (r.failed != nullptr)
                         && (std::strcmp(r.failed, expected_name) == 0);
    if (name_ok) {
        ++g_pass;
    } else {
        ++g_fail;
        std::fprintf(stderr,
                     "[FAIL] %s: failed=%s, expected %s\n",
                     label,
                     r.failed ? r.failed : "(null)",
                     expected_name);
    }

    expect_near(r.bad_value, expected_bad_value, 0.0, label);
}

Props16 textbook_params() {
    // From memory/project_umat_model.md (commented block at .for:122-138):
    //   phi=34, p_t=1, hs=3.6e6, en=0.43, ed0=0.72, ec0=0.934, ei0=1.2,
    //   alpha=0.24, beta=1.2, m_R=1, m_T=1, r_uc=3.3e-4, beta_r=0.5,
    //   chi=6.0, bulk_w=0, e_0=0.825
    return Props16{ {34.0, 1.0, 3.6e6, 0.43, 0.72, 0.934, 1.2, 0.24, 1.2,
                     1.0,  1.0, 3.3e-4, 0.5, 6.0, 0.0, 0.825} };
}

// Slot indices (0-based C++) matching sand_hypoplastic_kernel.h's
// 1-based-Fortran -> 0-based-C++ index map.
constexpr std::size_t kIdxPhi    = 0;
constexpr std::size_t kIdxPt     = 1;
constexpr std::size_t kIdxMR     = 9;
constexpr std::size_t kIdxMT     = 10;
constexpr std::size_t kIdxRuc    = 11;
constexpr std::size_t kIdxBetaR  = 12;
constexpr std::size_t kIdxChi    = 13;
constexpr std::size_t kIdxBulkW  = 14;

}  // namespace

int main(int argc, char** argv) {
    const std::string outdir = (argc > 1) ? argv[1]
                                          : std::string("..\\..\\..\\tests\\path_b_results");

    const Props16 props = textbook_params();
    const CheckParmsResult ok = check_parms(props);

    // --- Acceptance: textbook params ---
    expect_accepted(ok, "textbook accepted");

    const double pi = 4.0 * std::atan(1.0);
    expect_near(ok.parms[kIdxPhi], 34.0 * pi / 180.0, 1e-15,
                "textbook phi (rad)");

    // All other slots must pass through bit-identical.
    for (std::size_t i = 1; i < 16; ++i) {
        char label[64];
        std::snprintf(label, sizeof(label),
                      "textbook parms[%zu] passthrough", i);
        expect_near(ok.parms[i], props[i], 0.0, label);
    }

    // bulk_w == 0 accepted (textbook value, asserted again with full contract).
    expect_accepted(check_parms(props), "bulk_w=0 (textbook) accepted");

    // --- Rejection contract for every guarded slot ---

    // phi == 0 -> rejected. bad_value is the post-deg2rad value (= 0.0).
    {
        Props16 p = props; p[kIdxPhi] = 0.0;
        expect_rejected(check_parms(p), "phi", 0.0, "phi=0");
    }

    // phi < 0 -> rejected. bad_value is the post-deg2rad value.
    {
        Props16 p = props; p[kIdxPhi] = -1.0;
        const double phi_neg_rad = -1.0 * pi / 180.0;
        expect_rejected(check_parms(p), "phi", phi_neg_rad, "phi<0");
    }

    // m_R < 0 rejected.
    {
        Props16 p = props; p[kIdxMR] = -0.5;
        expect_rejected(check_parms(p), "m_R", -0.5, "m_R<0");
    }
    // m_R == 0 accepted.
    {
        Props16 p = props; p[kIdxMR] = 0.0;
        expect_accepted(check_parms(p), "m_R=0 accepted");
    }

    // m_T < 0 rejected.
    {
        Props16 p = props; p[kIdxMT] = -0.25;
        expect_rejected(check_parms(p), "m_T", -0.25, "m_T<0");
    }
    // m_T == 0 accepted.
    {
        Props16 p = props; p[kIdxMT] = 0.0;
        expect_accepted(check_parms(p), "m_T=0 accepted");
    }

    // r_uc < 0 rejected.
    {
        Props16 p = props; p[kIdxRuc] = -1.0e-4;
        expect_rejected(check_parms(p), "r_uc", -1.0e-4, "r_uc<0");
    }
    // r_uc == 0 with m_R > 0.5 (IS active) -> rejected as "r_uc"
    // by the cross-parameter constraint (textbook m_R = 1.0).
    {
        Props16 p = props; p[kIdxRuc] = 0.0;
        expect_rejected(check_parms(p), "r_uc", 0.0,
                        "r_uc=0 with m_R>0.5 (IS active) -> rejected");
    }
    // r_uc == 0 with m_R == 0.5 (boundary, IS NOT active per the
    // > 0.5 rule -- get_tan/rkf23/calc_statev all use strict `>`)
    // -> accepted: r_uc is unused in this regime, no division-by-zero.
    {
        Props16 p = props; p[kIdxRuc] = 0.0; p[kIdxMR] = 0.5;
        expect_accepted(check_parms(p),
                        "r_uc=0 with m_R=0.5 (boundary, IS off) accepted");
    }
    // r_uc == 0 with m_R = 0.4 (IS off) -> accepted: legal Fortran combo.
    {
        Props16 p = props; p[kIdxRuc] = 0.0; p[kIdxMR] = 0.4;
        expect_accepted(check_parms(p),
                        "r_uc=0 with m_R=0.4 (IS off) accepted");
    }
    // r_uc == 0 with m_R = 0 (IS off) -> accepted.
    {
        Props16 p = props; p[kIdxRuc] = 0.0; p[kIdxMR] = 0.0;
        expect_accepted(check_parms(p),
                        "r_uc=0 with m_R=0 accepted");
    }
    // r_uc > 0 with m_R > 0.5 (textbook) is the normal case -- already
    // covered by the textbook acceptance at the top, but assert the
    // cross-check does not spuriously trigger when r_uc is positive.
    {
        Props16 p = props; p[kIdxMR] = 1.0; p[kIdxRuc] = 1.0e-6;
        expect_accepted(check_parms(p),
                        "r_uc>0 with m_R>0.5 accepted (no spurious cross-check)");
    }

    // beta_r < 0 rejected.
    {
        Props16 p = props; p[kIdxBetaR] = -0.5;
        expect_rejected(check_parms(p), "beta_r", -0.5, "beta_r<0");
    }
    // beta_r == 0 accepted.
    {
        Props16 p = props; p[kIdxBetaR] = 0.0;
        expect_accepted(check_parms(p), "beta_r=0 accepted");
    }

    // chi < 0 rejected.
    {
        Props16 p = props; p[kIdxChi] = -2.0;
        expect_rejected(check_parms(p), "chi", -2.0, "chi<0");
    }
    // chi == 0 accepted.
    {
        Props16 p = props; p[kIdxChi] = 0.0;
        expect_accepted(check_parms(p), "chi=0 accepted");
    }

    // bulk_w < 0 rejected.
    {
        Props16 p = props; p[kIdxBulkW] = -1.0e3;
        expect_rejected(check_parms(p), "bulk_w", -1.0e3, "bulk_w<0");
    }

    // p_t < 0 rejected.
    {
        Props16 p = props; p[kIdxPt] = -1.0;
        expect_rejected(check_parms(p), "p_t", -1.0, "p_t<0");
    }
    // p_t == 0 accepted.
    {
        Props16 p = props; p[kIdxPt] = 0.0;
        expect_accepted(check_parms(p), "p_t=0 accepted");
    }

    // --- Ordering: when multiple slots are bad, the first checked wins ---
    // Fortran order: phi, m_R, m_T, r_uc, beta_r, chi, bulk_w, p_t.
    // If both m_T and chi are negative, m_T must be reported first.
    {
        Props16 p = props;
        p[kIdxMT]  = -1.0;
        p[kIdxChi] = -1.0;
        expect_rejected(check_parms(p), "m_T", -1.0,
                        "ordering: m_T reported before chi");
    }

    // --- Demo output: validated parms for the textbook params ---
    const std::string out_path = outdir + "\\check_parms_textbook.txt";
    std::ofstream out(out_path);
    if (!out) {
        std::fprintf(stderr,
                     "[WARN] could not open %s for writing; skipping demo file\n",
                     out_path.c_str());
    } else {
        out << "# Path B kernel: check_parms output for textbook params\n";
        out << "# Source: 00_hypo_is_original.for textbook block (lines 122-138)\n";
        out << "# Note: phi is reported in radians (kernel-internal).\n";
        out << "name\tvalue\n";
        static const char* names[16] = {
            "phi_rad", "p_t",  "hs",     "en",     "ed0",   "ec0",
            "ei0",     "alpha","beta",   "m_R",    "m_T",   "r_uc",
            "beta_r",  "chi",  "bulk_w", "e_0_raw"
        };
        out << std::scientific << std::setprecision(17);
        for (std::size_t i = 0; i < 16; ++i) {
            out << names[i] << "\t" << ok.parms[i] << "\n";
        }
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

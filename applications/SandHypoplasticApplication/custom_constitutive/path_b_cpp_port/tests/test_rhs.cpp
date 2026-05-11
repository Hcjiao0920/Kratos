// Standalone unit-test exerciser for the Path B kernel `rhs`.
//
// Compile + run via `build_and_run.bat` (sibling).
//
// `rhs` is a tiny composer over `get_F_sig_q`: it unpacks
// y[0..12] into (sig, q), calls get_F_sig_q, and packs F_sig + F_q
// back into y_dot[0..12]. The tests target the COMPOSITION
// (unpack/pack indices, error transparency) rather than the physics
// inside get_F_sig_q (which test_F_sig_q already owns).
//
// Coverage:
//   (1) Packing correctness. Build a y with known sig and q,
//       independently call get_F_sig_q to get F_sig and F_q, then
//       call rhs and assert y_dot[0..5] == F_sig and y_dot[6..12]
//       == F_q bit-exactly. A swap of the slice boundary or the
//       state offset (e.g. y[6+i] vs y[7+i]) would surface here.
//   (2) Error transparency. r_uc = 0 + textbook m_R = 1 surfaces
//       error = 10 from get_tan -> get_F_sig_q -> rhs. y_dot must
//       remain zero (no partial rate leak).
//   (3) No-IS slot pattern. Under m_R = 0.4, IS slots y_dot[6..11]
//       must be exactly zero (H[0..5] zero in get_tan's no-IS
//       branch); the void slot y_dot[12] must equal (1+v)*tr(deps).
//   (4) State-index map sanity. Use a non-axisymmetric deps so each
//       normal slot has a different value, and assert
//       y_dot[12] = (1+void) * (deps[0] + deps[1] + deps[2]).
//       This locks the void slot AT INDEX 12 (not 6 or 11), so an
//       off-by-one in the state map is caught immediately.

#include "../sand_hypoplastic_kernel.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>

using Kratos::SandHypoCpp::Parms16;
using Kratos::SandHypoCpp::Props16;
using Kratos::SandHypoCpp::RhsResult;
using Kratos::SandHypoCpp::YdotResult;
using Kratos::SandHypoCpp::check_parms;
using Kratos::SandHypoCpp::get_F_sig_q;
using Kratos::SandHypoCpp::rhs;
using Kratos::SandHypoCpp::kErrorOk;
using Kratos::SandHypoCpp::kStateDim;
using Kratos::SandHypoCpp::kYDim;

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

}  // namespace

int main() {
    const Parms16 parms = textbook_parms();
    const double void_ratio = 0.825;
    const double sig_compression[6] = {-100.0, -100.0, -100.0, 0.0, 0.0, 0.0};

    // Build the canonical y: stress in 0..5, intergranular strain in
    // 6..11 (zero IS state), void ratio in 12.
    auto build_y = [&](const double sig[6],
                       const double del[6],
                       double e,
                       double y_out[kYDim]) {
        for (int i = 0; i < 6; ++i) y_out[i]     = sig[i];
        for (int i = 0; i < 6; ++i) y_out[6 + i] = del[i];
        y_out[12] = e;
    };

    // ====================================================================
    // (1) Packing correctness: y_dot[0..5] == F_sig, y_dot[6..12] == F_q.
    //     Use IS branch (textbook m_R = 1) with a small compressive deps.
    // ====================================================================
    {
        const double del[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};

        double y[kYDim];
        build_y(sig_compression, del, void_ratio, y);

        const YdotResult yr = rhs(y, deps, parms);
        expect_eq_int(yr.error, kErrorOk, "case1 error == 0");

        // Independent reference from get_F_sig_q.
        const double q[kStateDim] = {del[0], del[1], del[2], del[3], del[4],
                                     del[5], void_ratio};
        const RhsResult ref = get_F_sig_q(deps, sig_compression, q, parms);
        expect_eq_int(ref.error, kErrorOk, "case1 reference error == 0");

        // y_dot[0..5] == F_sig.
        for (int i = 0; i < 6; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case1 y_dot[%d] == F_sig[%d] (stress slot pack)", i, i);
            expect_near(yr.y_dot[i], ref.F_sig[i], 0.0, label);
        }
        // y_dot[6..12] == F_q[0..6].
        for (std::size_t i = 0; i < kStateDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case1 y_dot[%zu] == F_q[%zu] (state slot pack)",
                          6 + i, i);
            expect_near(yr.y_dot[6 + i], ref.F_q[i], 0.0, label);
        }
    }

    // ====================================================================
    // (2) Error transparency: r_uc=0 + m_R=1 (IS-branch r_uc guard fires).
    //     y_dot must stay zero on the fatal path.
    //
    //     Note: after Codex 2026-05-02 (medium) hoisted "r_uc > 0 if
    //     m_R > 0.5" into check_parms, this combo is rejected at
    //     validation. We deliberately bypass check_parms by mutating
    //     a validated Parms16 directly to exercise get_tan's in-kernel
    //     belt-and-suspenders guard.
    // ====================================================================
    {
        Parms16 parms_runc_zero = textbook_parms();
        parms_runc_zero[11] = 0.0;          // bypass check_parms cross-check

        const double del[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y[kYDim];
        build_y(sig_compression, del, void_ratio, y);

        const YdotResult yr = rhs(y, deps, parms_runc_zero);
        expect_eq_int(yr.error, 10,
                      "case2 r_uc=0 + m_R>0.5 -> error=10 (propagated)");
        for (std::size_t i = 0; i < kYDim; ++i) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case2 y_dot[%zu] == 0 on fatal (no partial rate)", i);
            expect_near(yr.y_dot[i], 0.0, 0.0, label);
        }
    }

    // ====================================================================
    // (3) No-IS slot pattern: m_R = 0.4 (no-IS branch). y_dot[6..11] = 0,
    //     y_dot[12] = (1+v)*tr(deps).
    // ====================================================================
    {
        Props16 props_low_mR = textbook_props();
        props_low_mR[9] = 0.4;
        const auto check = check_parms(props_low_mR);

        const double del[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const double deps[6] = {-1.0e-3, -1.0e-3, -1.0e-3, 0.0, 0.0, 0.0};
        double y[kYDim];
        build_y(sig_compression, del, void_ratio, y);

        const YdotResult yr = rhs(y, deps, check.parms);
        expect_eq_int(yr.error, kErrorOk, "case3 error == 0 (no-IS)");

        // IS slots zero.
        for (int idx = 6; idx < 12; ++idx) {
            char label[80];
            std::snprintf(label, sizeof(label),
                          "case3 y_dot[%d] == 0 (IS slot in no-IS branch)", idx);
            expect_near(yr.y_dot[idx], 0.0, 0.0, label);
        }

        // Void slot.
        const double tr_deps = deps[0] + deps[1] + deps[2];
        const double y_dot_12_expected = (1.0 + void_ratio) * tr_deps;
        expect_near(yr.y_dot[12], y_dot_12_expected, 1e-15,
                    "case3 y_dot[12] == (1+v)*tr(deps) (void slot)");

        // Stress slots non-zero (sanity: rhs actually computed something).
        bool any_stress_nonzero = false;
        for (int i = 0; i < 6; ++i) {
            if (yr.y_dot[i] != 0.0) { any_stress_nonzero = true; break; }
        }
        if (any_stress_nonzero) {
            ++g_pass;
        } else {
            ++g_fail;
            std::fprintf(stderr,
                         "[FAIL] case3 sanity: y_dot[0..5] all zero (rhs didn't compute stress rate)\n");
        }
    }

    // ====================================================================
    // (4) State-index map sanity. Non-axisymmetric deps + non-textbook
    //     void ratio. The void slot is AT y_dot[12]; an off-by-one in
    //     packing (e.g. F_q packed into y_dot[7..13] instead of y_dot[6..12])
    //     would put the (1+v)*tr(deps) value at the wrong slot and fail.
    // ====================================================================
    {
        const double void_42 = 0.42;     // distinct from the textbook 0.825
        const double del[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        const double deps[6] = {-2.0e-3, -1.0e-3, -3.0e-3, 5.0e-4, 0.0, 0.0};

        double y[kYDim];
        build_y(sig_compression, del, void_42, y);

        const YdotResult yr = rhs(y, deps, parms);  // textbook m_R = 1
        expect_eq_int(yr.error, kErrorOk, "case4 error == 0");

        // Compute (1+v)*tr(deps) using the SAME void value passed in y[12].
        const double tr_deps = deps[0] + deps[1] + deps[2];
        const double y_dot_12_expected = (1.0 + void_42) * tr_deps;
        expect_near(yr.y_dot[12], y_dot_12_expected, 1e-15,
                    "case4 y_dot[12] == (1+0.42)*(deps[0]+deps[1]+deps[2])");

        // Negative cross-check: y_dot[11] (last IS slot) must NOT carry
        // the void rate. With zero-IS + m_R=1, H_del=identity in get_tan's
        // load<=0 else branch, so F_q[5] = deps[5] = 0 here. y_dot[11]
        // therefore equals 0, NOT (1+v)*tr(deps).
        expect_near(yr.y_dot[11], 0.0, 1e-15,
                    "case4 y_dot[11] == 0 (last IS slot, void rate not bleeding into wrong slot)");

        // y_dot[6] should equal deps[0] (zero-IS H_del=identity, F_q[0]=deps[0]).
        expect_near(yr.y_dot[6], deps[0], 1e-15,
                    "case4 y_dot[6] == deps[0] (first IS slot, zero-IS H_del=identity)");
    }

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return (g_fail == 0) ? 0 : 1;
}

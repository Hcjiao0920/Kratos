// Path B: kernel implementation. See sand_hypoplastic_kernel.h.
//
// Each routine mirrors a Fortran subroutine in `00_hypo_is_original.for`:
//   check_parms   <-> check_parms_h    (.for:389-508)
//   dot_vect      <-> dot_vect_h       (.for:532-579)
//   matmul        <-> matmul_h         (.for:1233-1253)
//   inv_eps       <-> inv_eps_h        (.for:1026-1104)
//   inv_sig       <-> inv_sig_h        (.for:1106-1231)
//   get_tan            <-> get_tan_h        (.for:631-997)
//   get_F_sig_q        <-> get_F_sig_q_h    (.for:581-629)
//   rhs                <-> rhs_h            (.for:1511-1567)
//   principal_stresses_3 <- replaces Eig_3a_h / PrnSig_h with closed-form Cardano
//   check_RKF          <-> check_RKF_h      (.for:1768-1831)
//   calc_elasti        <-> calc_elasti_h    (.for:2372-2453)
//   norm_res           <-> norm_res_h       (.for:1342-1421)
//   rkf23_update       <-> rkf23_update_h   (.for:1569-1764)
//   integrate_step     <-> umat top-level driver structure
//                          (.for:19-386, minus testing/pore/statev plumbing)

#include "sand_hypoplastic_kernel.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>

namespace Kratos {
namespace SandHypoCpp {

namespace {

constexpr double kUmatisnanLimit = 1.0e30;

bool fails_umatisnan_h(double value) {
    return !std::isfinite(value)
        || value > kUmatisnanLimit
        || value < -kUmatisnanLimit;
}

}  // namespace

CheckParmsResult check_parms(const Props16& props) {
    CheckParmsResult result;
    result.parms = props;
    result.error = kErrorOk;
    result.failed = nullptr;
    result.bad_value = std::numeric_limits<double>::quiet_NaN();

    // Fortran: pi = four * datan(one); phi = phi_deg * pi / 180.0
    const double pi = 4.0 * std::atan(1.0);
    const double phi_deg = result.parms[0];
    const double phi = phi_deg * pi / 180.0;
    result.parms[0] = phi;

    auto reject = [&result](const char* name, double value) {
        result.error = kErrorFatal;
        result.failed = name;
        result.bad_value = value;
    };

    if (phi <= 0.0)              { reject("phi",    phi);              return result; }
    if (result.parms[9]  < 0.0)  { reject("m_R",    result.parms[9]);  return result; }
    if (result.parms[10] < 0.0)  { reject("m_T",    result.parms[10]); return result; }
    if (result.parms[11] < 0.0)  { reject("r_uc",   result.parms[11]); return result; }
    if (result.parms[12] < 0.0)  { reject("beta_r", result.parms[12]); return result; }
    if (result.parms[13] < 0.0)  { reject("chi",    result.parms[13]); return result; }
    if (result.parms[14] < 0.0)  { reject("bulk_w", result.parms[14]); return result; }
    if (result.parms[1]  < 0.0)  { reject("p_t",    result.parms[1]);  return result; }

    // Cross-parameter constraint: when intergranular strain is active
    // (m_R > 0.5, the same threshold get_tan / rkf23_update / calc_statev
    // use to select the IS branch), r_uc must be strictly positive.
    // The IS branch in get_tan computes `rho = norm_del / r_uc` and the
    // Fortran reference would silently NaN-poison when r_uc = 0; we
    // already reject it inside get_tan as a defensive deviation, but
    // catching it here at validation time is strictly stronger:
    //   - first-failure-wins ordering: bad material rejected before any
    //     integration runs, not at the first step that happens to use
    //     the IS branch;
    //   - eliminates path-dependence in `integrate_step`: previously a
    //     tensile initial sig would route to the elastic fallback
    //     (no r_uc dependency), return error=0, then the next
    //     compressive step would emit error=10 from the IS guard --
    //     the user would see a fatal "out of nowhere" mid-trajectory.
    //     Hoisting here makes the verdict stress-state-independent.
    //   - makes the get_tan defensive guard structurally unreachable
    //     through validated parms (kept as belt-and-suspenders).
    // Codex review 2026-05-02 (medium follow-up: path-dependence).
    if (result.parms[9] > 0.5 && result.parms[11] <= 0.0) {
        reject("r_uc", result.parms[11]);
        return result;
    }

    return result;
}

double dot_vect(DotProductKind kind,
                const double* a,
                const double* b,
                std::size_t n) {
    // Initialized to 1.0 so an out-of-range enum value (e.g. an integer
    // cast outside the enum's enumerators) degrades to Plain semantics
    // rather than reading an uninitialized double. Every defined
    // enumerator is handled below.
    double coeff = 1.0;
    switch (kind) {
        case DotProductKind::StressLike: coeff = 2.0; break;
        case DotProductKind::StrainLike: coeff = 0.5; break;
        case DotProductKind::Plain:      coeff = 1.0; break;
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (i < 3) {
            sum += a[i] * b[i];
        } else {
            sum += coeff * a[i] * b[i];
        }
    }
    return sum;
}

void matmul(const double* a,
            const double* b,
            double* c,
            std::size_t l,
            std::size_t m,
            std::size_t n) {
    // Row-major: a is (l x m), b is (m x n), c is (l x n).
    // Precondition: c does not alias a or b. Not enforced at runtime.
    for (std::size_t i = 0; i < l; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < m; ++k) {
                s += a[i * m + k] * b[k * n + j];
            }
            c[i * n + j] = s;
        }
    }
}

StrainInvariants inv_eps(const double* eps) {
    StrainInvariants r;

    const double onethird  = 1.0 / 3.0;
    const double twothirds = 2.0 / 3.0;
    const double sqrt6     = std::sqrt(6.0);

    // Volumetric strain: eps_v = trace(eps).
    r.eps_v = eps[0] + eps[1] + eps[2];
    const double ev3 = onethird * r.eps_v;

    // Deviator strain. Off-diagonal slots in the input are engineering
    // shear strain (gamma_ij = 2*eps_ij), so divide by two to obtain
    // tensor strain components for the deviator.
    double edev[6];
    edev[0] = eps[0] - ev3;
    edev[1] = eps[1] - ev3;
    edev[2] = eps[2] - ev3;
    edev[3] = eps[3] / 2.0;
    edev[4] = eps[4] / 2.0;
    edev[5] = eps[5] / 2.0;

    // Norm-squared in tensor sense: ||edev||^2 = e_ij e_ij. In Voigt
    // (with edev already storing tensor components for off-diagonals)
    // the off-diagonal terms enter twice each.
    const double norm2 = edev[0]*edev[0] + edev[1]*edev[1] + edev[2]*edev[2]
                       + 2.0 * (edev[3]*edev[3] + edev[4]*edev[4] + edev[5]*edev[5]);

    r.eps_s = std::sqrt(twothirds * norm2);

    // Lode angle: undefined for purely hydrostatic strain.
    if (r.eps_s == 0.0) {
        r.sin3t = -1.0;
        return r;
    }

    // Components of (edev_ij)(edev_jk) in Voigt with the factor of 2 on
    // off-diagonals BAKED IN. This lets `tredev3` be a plain dot product
    // rather than a Voigt-weighted one.
    double edev2[6];
    edev2[0] = edev[0]*edev[0] + edev[3]*edev[3] + edev[4]*edev[4];
    edev2[1] = edev[3]*edev[3] + edev[1]*edev[1] + edev[5]*edev[5];
    edev2[2] = edev[5]*edev[5] + edev[4]*edev[4] + edev[2]*edev[2];
    edev2[3] = 2.0 * (edev[0]*edev[3] + edev[3]*edev[1] + edev[5]*edev[4]);
    edev2[4] = 2.0 * (edev[4]*edev[0] + edev[5]*edev[3] + edev[2]*edev[4]);
    edev2[5] = 2.0 * (edev[3]*edev[4] + edev[1]*edev[5] + edev[5]*edev[2]);

    double tredev3 = 0.0;
    for (int i = 0; i < 6; ++i) {
        tredev3 += edev[i] * edev2[i];
    }

    const double numer = sqrt6 * tredev3;
    const double sqrt_norm2 = std::sqrt(norm2);
    const double denom = sqrt_norm2 * sqrt_norm2 * sqrt_norm2;
    r.sin3t = numer / denom;

    // Clamp |sin3t| to <= 1 against floating-point overshoot, sign-preserving.
    if (std::abs(r.sin3t) > 1.0) {
        r.sin3t = r.sin3t / std::abs(r.sin3t);
    }

    return r;
}

StressInvariants inv_sig(const double* sig) {
    StressInvariants r;

    const double onethird    = 1.0 / 3.0;
    const double threehalves = 3.0 / 2.0;
    const double sqrt6       = std::sqrt(6.0);
    // NOTE: this `tiny` is 1e-18, smaller than the 1e-17 used in
    // get_tan_h. Matches the Fortran inv_sig_h verbatim (.for:1136).
    const double tiny = 1.0e-18;

    // Trace and mean stress.
    r.I1 = sig[0] + sig[1] + sig[2];
    r.pp = onethird * r.I1;

    // Deviator stress. Off-diagonals carry tensor shear directly
    // (no /2 conversion, unlike inv_eps with engineering strain).
    double sdev[6];
    sdev[0] = sig[0] - r.pp;
    sdev[1] = sig[1] - r.pp;
    sdev[2] = sig[2] - r.pp;
    sdev[3] = sig[3];
    sdev[4] = sig[4];
    sdev[5] = sig[5];

    // Normalized stress eta = sig / I1. When I1 == 0, the Fortran
    // substitutes sig / tiny verbatim; this is a numerical hack, not a
    // physically meaningful normalization. See header note.
    double eta[6];
    if (r.I1 != 0.0) {
        for (int i = 0; i < 6; ++i) eta[i] = sig[i] / r.I1;
    } else {
        for (int i = 0; i < 6; ++i) eta[i] = sig[i] / tiny;
    }

    // Deviatoric eta.
    double eta_d[6];
    eta_d[0] = eta[0] - onethird;
    eta_d[1] = eta[1] - onethird;
    eta_d[2] = eta[2] - onethird;
    eta_d[3] = eta[3];
    eta_d[4] = eta[4];
    eta_d[5] = eta[5];

    // Second invariants via Voigt-stress dot product.
    const double norm2    = dot_vect(DotProductKind::StressLike, sdev,  sdev,  6);
    const double norm2sig = dot_vect(DotProductKind::StressLike, sig,   sig,   6);
    const double norm2eta = dot_vect(DotProductKind::StressLike, eta_d, eta_d, 6);

    r.qq = std::sqrt(threehalves * norm2);
    r.I2 = 0.5 * (norm2sig - r.I1 * r.I1);

    // Components of (eta_d_ij)(eta_d_jk) in Voigt WITHOUT the factor of
    // 2 on off-diagonals; the factor is supplied by `dot_vect` in the
    // tretadev3 contraction below.
    double eta_d2[6];
    eta_d2[0] = eta_d[0]*eta_d[0] + eta_d[3]*eta_d[3] + eta_d[4]*eta_d[4];
    eta_d2[1] = eta_d[3]*eta_d[3] + eta_d[1]*eta_d[1] + eta_d[5]*eta_d[5];
    eta_d2[2] = eta_d[5]*eta_d[5] + eta_d[4]*eta_d[4] + eta_d[2]*eta_d[2];
    eta_d2[3] = eta_d[0]*eta_d[3] + eta_d[3]*eta_d[1] + eta_d[5]*eta_d[4];
    eta_d2[4] = eta_d[4]*eta_d[0] + eta_d[5]*eta_d[3] + eta_d[2]*eta_d[4];
    eta_d2[5] = eta_d[3]*eta_d[4] + eta_d[1]*eta_d[5] + eta_d[5]*eta_d[2];

    // Lode angle.
    if (norm2eta < tiny) {
        r.cos3t = -1.0;
    } else {
        const double tretadev3 =
            dot_vect(DotProductKind::StressLike, eta_d, eta_d2, 6);
        const double numer = -sqrt6 * tretadev3;
        const double sqrt_norm2eta = std::sqrt(norm2eta);
        const double denom = sqrt_norm2eta * sqrt_norm2eta * sqrt_norm2eta;
        r.cos3t = numer / denom;
        // Clamp |cos3t| to <= 1, sign-preserving.
        if (std::abs(r.cos3t) > 1.0) {
            r.cos3t = r.cos3t / std::abs(r.cos3t);
        }
    }

    // Determinant I3 (cofactor expansion along row 1, in Abaqus order).
    const double xmin1 = sig[1]*sig[2] - sig[5]*sig[5];
    const double xmin2 = sig[3]*sig[2] - sig[5]*sig[4];
    const double xmin3 = sig[3]*sig[5] - sig[4]*sig[1];
    r.I3 = sig[0]*xmin1 - sig[3]*xmin2 + sig[4]*xmin3;

    return r;
}

TangentResult get_tan(const double* deps,
                      const double* sig,
                      const double* q,
                      const Parms16& parms,
                      int istrain) {
    TangentResult r{};
    r.error = kErrorOk;

    // ---- Constants (Fortran .for:666-681) -----------------------------------
    const double tiny     = 1.0e-17;        // get_tan_h's tiny (NOT inv_sig's)
    const double half     = 0.5;
    const double onethird = 1.0 / 3.0;
    const double sqrt2    = std::sqrt(2.0);
    const double sqrt3    = std::sqrt(3.0);
    const double twosqrt2 = 2.0 * sqrt2;
    const double oneeight = 1.0 / 8.0;

    // ---- Fourth-order identity tensors in Voigt notation -------------------
    // II is the symmetric Voigt identity: II[i][j] is 1 on the diagonal
    // for i<3, 1/2 on the diagonal for i>=3, zero off-diagonal.
    // IU is the plain Voigt identity (1 on every diagonal slot).
    double II[6][6] = {{0.0}};
    double IU[6][6] = {{0.0}};
    for (int i = 0; i < 3; ++i) II[i][i] = 1.0;
    for (int i = 3; i < 6; ++i) II[i][i] = half;
    for (int i = 0; i < 6; ++i) IU[i][i] = 1.0;

    // ---- Material parameters (Fortran .for:720-734) ------------------------
    const double phi    = parms[0];
    const double p_t    = parms[1];
    const double hs     = parms[2];
    const double en     = parms[3];
    const double ed0    = parms[4];
    const double ec0    = parms[5];
    const double ei0    = parms[6];
    const double alpha  = parms[7];
    const double beta   = parms[8];
    const double m_R    = parms[9];
    const double m_T    = parms[10];
    const double r_uc   = parms[11];
    const double beta_r = parms[12];
    const double chi    = parms[13];
    // parms[14] = bulk_w (unused in get_tan), parms[15] = e_0 (unused).

    // ---- Internal state (Fortran .for:742-748) -----------------------------
    double del[6];
    for (int i = 0; i < 6; ++i) del[i] = q[i];
    const double void_ratio = q[6];

    // ---- Cohesion-shifted stress (Fortran .for:752-757) --------------------
    double sig_star[6];
    sig_star[0] = sig[0] - p_t;
    sig_star[1] = sig[1] - p_t;
    sig_star[2] = sig[2] - p_t;
    sig_star[3] = sig[3];
    sig_star[4] = sig[4];
    sig_star[5] = sig[5];

    // ---- Strain-rate / IS norms and unit vectors (Fortran .for:761-787) ----
    const double norm_deps2 = dot_vect(DotProductKind::StrainLike, deps, deps, 6);
    const double norm_del2  = dot_vect(DotProductKind::StrainLike, del,  del,  6);
    const double norm_deps  = std::sqrt(norm_deps2);
    const double norm_del   = std::sqrt(norm_del2);

    double eta_del[6]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double eta_delta[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double eta_eps[6]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (norm_del >= tiny) {
        for (int i = 0; i < 6; ++i) eta_del[i] = del[i] / norm_del;
    }

    eta_delta[0] = eta_del[0];
    eta_delta[1] = eta_del[1];
    eta_delta[2] = eta_del[2];
    eta_delta[3] = half * eta_del[3];
    eta_delta[4] = half * eta_del[4];
    eta_delta[5] = half * eta_del[5];

    if (norm_deps >= tiny) {
        for (int i = 0; i < 6; ++i) eta_eps[i] = deps[i] / norm_deps;
    }

    // ---- Stress invariants of sig_star (Fortran .for:791) ------------------
    const StressInvariants si = inv_sig(sig_star);
    const double I1    = si.I1;
    const double cos3t = si.cos3t;

    // ---- Normalized stress eta and deviator eta_dev (Fortran .for:805-817)
    // Note: this divides by I1 unconditionally, matching the Fortran. The
    // commented-out tensile guard at .for:793-802 is preserved as a
    // comment in the original; we honour the post-Tamagnini umat design
    // where the caller short-circuits via inittension before reaching us.
    double eta[6];
    double eta_dev[6];
    for (int i = 0; i < 6; ++i) eta[i] = sig_star[i] / I1;
    eta_dev[0] = eta[0] - onethird;
    eta_dev[1] = eta[1] - onethird;
    eta_dev[2] = eta[2] - onethird;
    eta_dev[3] = eta[3];
    eta_dev[4] = eta[4];
    eta_dev[5] = eta[5];

    // ---- Functions a and FF (Fortran .for:821-829) -------------------------
    const double eta_dn2 = dot_vect(DotProductKind::StressLike, eta_dev, eta_dev, 6);
    const double tanpsi  = sqrt3 * std::sqrt(eta_dn2);
    const double tanpsi2 = tanpsi * tanpsi;
    const double FF_radicand =
        oneeight * tanpsi2
        + (2.0 - tanpsi2) / (2.0 + sqrt2 * tanpsi * cos3t);
    const double FF_minus     = tanpsi / twosqrt2;
    const double sinphi       = std::sin(phi);
    const double a            = sqrt3 * (3.0 - sinphi) / (twosqrt2 * sinphi);
    const double a2           = a * a;
    const double FF           = std::sqrt(FF_radicand) - FF_minus;

    // ---- Barotropy and pyknotropy (Fortran .for:833-849) -------------------
    // bauer = exp(-(-I1/hs)^en). Mirrors the Fortran exactly: for I1 >= 0
    // (tension/zero), -I1/hs is non-positive and `pow` produces NaN; the
    // caller is responsible for not reaching get_tan in that regime.
    const double bauer = std::exp(-std::pow(-I1 / hs, en));
    const double ed    = ed0 * bauer;
    const double ec    = ec0 * bauer;
    const double ei    = ei0 * bauer;

    // fb temp1 (.for:838): factor inside the fb denominator.
    const double fb_temp1 =
        3.0 + a2 - a * sqrt3
              * std::pow((ei0 - ed0) / (ec0 - ed0), alpha);
    if (fb_temp1 < 0.0) {
        // .for:839: `if(temp1.lt.zero) stop 'factor fb not defined'`.
        // We translate the STOP into an error code; M/L/N/H stay zero
        // (already initialized via {} above).
        r.error = kErrorFatal;
        return r;
    }

    const double fb = (hs / en / fb_temp1)
                    * (1.0 + ei) / ei
                    * std::pow(ei0 / ec0, beta)
                    * std::pow(-I1 / hs, 1.0 - en);
    const double fe = std::pow(ec / void_ratio, beta);
    const double fs = fb * fe;

    double fd;
    if (void_ratio >= ed) {
        fd = std::pow((void_ratio - ed) / (ec - ed), alpha);
    } else {
        fd = 0.0;
    }

    // ---- Tensor L (Fortran .for:854-860) -----------------------------------
    const double eta_n2 = dot_vect(DotProductKind::StressLike, eta, eta, 6);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            r.L[i][j] = (II[i][j] * FF * FF
                         + a2 * eta[i] * eta[j]) / eta_n2;
        }
    }

    // ---- Vector NN (Fortran .for:866-868) ----------------------------------
    for (int i = 0; i < 6; ++i) {
        r.N[i] = FF * a * (eta[i] + eta_dev[i]) / eta_n2;
    }

    // =====================================================================
    // Branch on istrain
    // =====================================================================
    if (istrain == 1) {
        // Belt-and-suspenders guard: structurally unreachable through
        // validated parms after Codex 2026-05-02 (medium follow-up)
        // hoisted "r_uc > 0 if m_R > 0.5" into `check_parms`. Kept
        // here as defense-in-depth in case a caller constructs a
        // Parms16 that bypasses check_parms (e.g. via memcpy or a
        // future test seam) -- it preserves the structural invariant
        // that the IS branch never divides by zero.
        //
        // Original rationale (Fortran deviation, Codex 2026-05-01
        // high finding): the Fortran reference computes
        // `rho = norm_del / r_uc` unconditionally; with the common
        // zero-IS initial state and an accepted r_uc == 0 (or
        // negative-passing-through-FP) this yields NaN that silently
        // propagates into M and H while error stays 0.
        if (r_uc <= 0.0) {
            r.error = kErrorFatal;
            return r;
        }

        // ---- Loading function (Fortran .for:878) ---------------------------
        const double load = dot_vect(DotProductKind::StrainLike,
                                     eta_del, eta_eps, 6);

        // ---- Intergranular-strain rho (Fortran .for:882-886) ---------------
        double rho = norm_del / r_uc;
        if (rho > 1.0) rho = 1.0;

        const double rho_chi    = std::pow(rho, chi);
        const double rho_beta_r = std::pow(rho, beta_r);

        // Leta = LL @ eta_del (Fortran .for:888) -----------------------------
        double Leta[6];
        matmul(&r.L[0][0], eta_del, Leta, 6, 6, 1);

        // ---- Tangent stiffness MM (Fortran .for:892-918) -------------------
        const double mm_temp1 =
            (rho_chi * m_T + (1.0 - rho_chi) * m_R) * fs;

        if (load > 0.0) {
            const double mm_temp2 = rho_chi * (1.0 - m_T) * fs;
            const double mm_temp3 = rho_chi * fs * fd;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double AA = mm_temp2 * Leta[i] * eta_delta[j]
                                    + mm_temp3 * r.N[i] * eta_delta[j];
                    r.M[i][j] = mm_temp1 * r.L[i][j] + AA;
                }
            }
        } else {
            const double mm_temp4 = rho_chi * (m_R - m_T) * fs;
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    const double AA = mm_temp4 * Leta[i] * eta_delta[j];
                    r.M[i][j] = mm_temp1 * r.L[i][j] + AA;
                }
            }
        }

        // ---- IS evolution H_del (Fortran .for:925-939) ---------------------
        // NOTE: Fortran reuses H_del[i][j] but only zeroes diagonal in the
        // else branch. Off-diagonals come from the zero-init at the top
        // of the routine. We mirror that with explicit zero-init.
        double H_del[6][6] = {{0.0}};
        if (load > 0.0) {
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    H_del[i][j] = IU[i][j]
                                  - rho_beta_r * eta_del[i] * eta_delta[j];
                }
            }
        } else {
            for (int i = 0; i < 6; ++i) H_del[i][i] = 1.0;
        }

        // ---- Void-ratio evolution H_e (Fortran .for:943-949) ---------------
        double H_e[6];
        for (int i = 0; i < 6; ++i) {
            H_e[i] = (i < 3) ? (1.0 + void_ratio) : 0.0;
        }

        // ---- Assemble HH (Fortran .for:953-963) ----------------------------
        for (std::size_t i = 0; i < kStateDim; ++i) {
            if (i < 6) {
                for (int j = 0; j < 6; ++j) r.H[i][j] = H_del[i][j];
            } else {
                for (int j = 0; j < 6; ++j) r.H[i][j] = H_e[j];
            }
        }
    } else if (istrain == 0) {
        // ---- No-IS branch (Fortran .for:966-985) ---------------------------
        double H_e[6];
        for (int i = 0; i < 6; ++i) {
            H_e[i] = (i < 3) ? (1.0 + void_ratio) : 0.0;
        }
        for (std::size_t i = 0; i < kStateDim; ++i) {
            if (i < 6) {
                for (int j = 0; j < 6; ++j) r.H[i][j] = 0.0;
            } else {
                for (int j = 0; j < 6; ++j) r.H[i][j] = H_e[j];
            }
        }
        // M is left zero-initialized.
    }

    // ---- Final scaling of L and N (Fortran .for:989-994) -------------------
    // Note: this happens AFTER M is built using the unscaled L and N,
    // but the temp1/temp2/temp3/temp4 inside the IS branch have already
    // multiplied by fs (and fd where appropriate), so M is correctly
    // scaled. The scaling here only affects the L/N OUTPUTS.
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            r.L[i][j] *= fs;
        }
        r.N[i] *= fs * fd;
    }

    return r;
}

RhsResult get_F_sig_q(const double* deps,
                      const double* sig,
                      const double* q,
                      const Parms16& parms) {
    RhsResult r{};

    // Fortran-pinned istrain selection (.for:602-606): m_R <= 0.5 -> 0,
    // otherwise 1. m_R lives at parms[9] (Fortran parms(10)).
    const int istrain = (parms[9] <= 0.5) ? 0 : 1;

    const TangentResult tan_r = get_tan(deps, sig, q, parms, istrain);
    if (tan_r.error != kErrorOk) {
        // Propagate the kernel-level fatal verbatim. F_sig and F_q
        // remain zero-initialized; never leak a partial rate vector.
        r.error = tan_r.error;
        return r;
    }
    r.error = kErrorOk;

    if (istrain == 1) {
        // F_sig = M * deps  (.for:614)
        matmul(&tan_r.M[0][0], deps, r.F_sig, 6, 6, 1);
    } else {
        // F_sig = L * deps + N * |deps|_strain  (.for:615-621)
        matmul(&tan_r.L[0][0], deps, r.F_sig, 6, 6, 1);
        const double norm_D2 = dot_vect(DotProductKind::StrainLike,
                                        deps, deps, 6);
        const double norm_D  = std::sqrt(norm_D2);
        for (int i = 0; i < 6; ++i) {
            r.F_sig[i] += tan_r.N[i] * norm_D;
        }
    }

    // F_q = H * deps  (.for:626)
    matmul(&tan_r.H[0][0], deps, r.F_q, kStateDim, 6, 1);

    return r;
}

YdotResult rhs(const double* y,
               const double* deps,
               const Parms16& parms) {
    YdotResult r{};  // y_dot zero-initialized

    // Unpack y -> (sig, q) (Fortran .for:1541-1547). We pass slices of
    // y directly to get_F_sig_q; both functions expect a const double*
    // 6-vector and 7-vector respectively, and y's storage is contiguous.
    const double* sig = &y[0];   // y[0..5]   = stress
    const double* q   = &y[6];   // y[6..12]  = state (del[0..5] + void)

    const RhsResult fr = get_F_sig_q(deps, sig, q, parms);
    if (fr.error != kErrorOk) {
        // Propagate fatal verbatim; y_dot stays zero-initialized.
        // The Fortran does the same via `if(error.eq.10) return`.
        r.error = fr.error;
        return r;
    }
    r.error = kErrorOk;

    // Pack F_sig into y_dot[0..5] and F_q into y_dot[6..12]
    // (Fortran .for:1554-1564).
    for (int i = 0; i < 6; ++i) {
        r.y_dot[i] = fr.F_sig[i];
    }
    for (std::size_t i = 0; i < kStateDim; ++i) {
        r.y_dot[6 + i] = fr.F_q[i];
    }

    return r;
}

PrincipalStresses principal_stresses_3(const double* sig) {
    // sig in Abaqus-Voigt order: [11, 22, 33, 12, 13, 23].
    PrincipalStresses ps{};

    // Off-diagonal magnitude squared.
    const double p1 = sig[3]*sig[3] + sig[4]*sig[4] + sig[5]*sig[5];

    if (p1 == 0.0) {
        // Diagonal matrix: eigenvalues are the diagonal entries.
        // Sort ascending with three swaps.
        double t0 = sig[0], t1 = sig[1], t2 = sig[2];
        if (t0 > t1) { const double tmp = t0; t0 = t1; t1 = tmp; }
        if (t1 > t2) { const double tmp = t1; t1 = t2; t2 = tmp; }
        if (t0 > t1) { const double tmp = t0; t0 = t1; t1 = tmp; }
        ps.s[0] = t0;
        ps.s[1] = t1;
        ps.s[2] = t2;
        return ps;
    }

    // Closed-form Cardano for symmetric 3x3.
    //
    // q  = trace(A) / 3                           (mean)
    // p2 = ||A - q I||_F^2                        (squared Frobenius of deviator)
    // p  = sqrt(p2 / 6)
    // det_dev = det(A - q I)
    // r  = det_dev / (2 * p^3)                    (in [-1, 1] mathematically)
    // phi = acos(clamp(r, -1, 1)) / 3             (in [0, pi/3])
    //
    // Eigenvalues are q + 2*p*cos(phi + 2*pi*k/3) for k = 0, 1, 2;
    // ordered: q + 2 p cos(phi)         largest  (s[2])
    //          q + 2 p cos(phi + 2pi/3) smallest (s[0])
    //          eig2 = 3q - eig1 - eig3 middle    (s[1])
    const double q  = (sig[0] + sig[1] + sig[2]) / 3.0;
    const double d0 = sig[0] - q;
    const double d1 = sig[1] - q;
    const double d2 = sig[2] - q;
    const double p2 = d0*d0 + d1*d1 + d2*d2 + 2.0 * p1;
    const double p  = std::sqrt(p2 / 6.0);

    // det of (A - q I) computed directly in Voigt.
    const double det_dev = d0 * (d1*d2 - sig[5]*sig[5])
                         - sig[3] * (sig[3]*d2 - sig[4]*sig[5])
                         + sig[4] * (sig[3]*sig[5] - sig[4]*d1);

    double r_val = (det_dev / (p * p * p)) / 2.0;
    // Clamp against floating-point overshoot of [-1, 1].
    if (r_val >  1.0) r_val =  1.0;
    if (r_val < -1.0) r_val = -1.0;

    const double pi  = 4.0 * std::atan(1.0);
    const double phi = std::acos(r_val) / 3.0;

    const double eig_max = q + 2.0 * p * std::cos(phi);
    const double eig_min = q + 2.0 * p * std::cos(phi + 2.0 * pi / 3.0);
    const double eig_mid = 3.0 * q - eig_max - eig_min;

    ps.s[0] = eig_min;
    ps.s[1] = eig_mid;
    ps.s[2] = eig_max;
    return ps;
}

CheckRkfResult check_RKF(const double* y, const Parms16& parms) {
    CheckRkfResult r{};
    r.error  = 0;
    r.failed = nullptr;

    const double p_t       = parms[1];
    const double minstress = p_t / 4.0;

    // sig from y[0..5].
    const double sig[6] = {y[0], y[1], y[2], y[3], y[4], y[5]};

    // Cohesion-shifted stress: sig_star = sig - p_t * I.
    const double sig_star[6] = {
        sig[0] - p_t,
        sig[1] - p_t,
        sig[2] - p_t,
        sig[3],
        sig[4],
        sig[5]
    };

    // Failed-label precedence (deviation from Fortran's single-bit
    // error_RKF return; Codex review 2026-05-02 P1 finding):
    //   - "pmean" / "tension" classify a STRESS-DOMAIN failure that
    //     `integrate_step` treats as legitimate inittension and routes
    //     through `calc_elasti`.
    //   - "nan" / "void" classify a NUMERIC-CORRUPTION / contract
    //     failure that `integrate_step` MUST fail-closed on (returning
    //     error=3 with state preserved); silently absorbing it as
    //     elastic would let overflow-finite or NaN state persist
    //     across step boundaries.
    // Therefore corruption labels UNCONDITIONALLY override stress-
    // domain labels. The previous "first-failure-wins" preservation
    // had a documented leak: a derived `min_principal` overflow
    // (e.g. sig = [0,0,0,1e30,1e30,1e30] -> principal eigenvalues
    // ~2e30 -> min_principal = -2e30 trips both the "tension"
    // threshold AND the umatisnan_h overflow sentinel) would set
    // `failed = "tension"` first, mask `failed = "nan"` from the
    // overflow check, and slip through integrate_step as an elastic
    // fallback. With unconditional override, the corruption class
    // wins regardless of ordering.
    //
    // Fortran reference behavior (.for:1768-1831): returns a single
    // bit, no label, all failures conflated; umat (.for:259-340)
    // routes EVERY failure through calc_elasti including overflow.
    // Path B's stricter "fail closed on corruption" behavior is the
    // logical closure of the integrate_step-side defensive deviations
    // already in place (raw-y NaN/Inf preflight, deps norm guard).
    //
    // Within each class, first-failure-wins still applies (pmean
    // takes precedence over tension; nan over void).

    // Mean stress (compression-positive). Fortran: pmean = -tr(sig_star)/3.
    r.pmean = -(sig_star[0] + sig_star[1] + sig_star[2]) / 3.0;
    if (r.pmean <= minstress) {
        r.error  = 1;
        r.failed = "pmean";
    }

    // Principal stresses; tmin = -max(S) (compression-positive), so the
    // most-tensile principal of sig_star fails the threshold.
    const PrincipalStresses ps = principal_stresses_3(sig_star);
    r.min_principal = -ps.s[2];  // s[2] is the largest (most-tensile in tension-positive)
    if (r.min_principal <= minstress) {
        if (r.error == 0) {
            r.error  = 1;
            r.failed = "tension";
        }
        // Within stress-domain class, "pmean" already wins; nothing to do.
    }

    // Fortran umatisnan_h scan over y(1..ny) and tmin. Despite the name it
    // also rejects +/-Inf and finite overflow sentinels beyond +/-1.d30.
    // CORRUPTION CLASS: unconditional override of any prior stress-domain
    // label.
    for (std::size_t i = 0; i < kYDim; ++i) {
        if (fails_umatisnan_h(y[i])) {
            r.error  = 1;
            r.failed = "nan";   // override stress-domain
            break;
        }
    }
    if (fails_umatisnan_h(r.min_principal)) {
        r.error  = 1;
        // "nan" wins over pmean/tension; if a prior raw-y scan already
        // set "nan", this re-set is idempotent.
        r.failed = "nan";
    }

    // Void-ratio admissibility (deviation from Fortran). The Fortran
    // check_RKF_h does not check void > 0; downstream get_tan computes
    // `fe = pow(ec / void_ratio, beta)`, so void <= 0 produces NaN/Inf
    // in fe -> fs -> M / L / N, which then escape via an accepted
    // y_hat whose stress was still admissible. NaN/Inf void is already
    // caught by the umatisnan_h scan above; this guard handles the
    // FINITE-but-non-positive case (e.g. y_hat[12] = -0.5 after a large
    // volumetric-compression substep).
    // Codex review 2026-05-01 (medium follow-up).
    // CORRUPTION CLASS: overrides pmean/tension; defers to a prior "nan"
    // (the existing first-failure-wins within the corruption class).
    if (y[12] <= 0.0) {
        r.error = 1;
        if (r.failed == nullptr
         || std::strcmp(r.failed, "pmean")   == 0
         || std::strcmp(r.failed, "tension") == 0) {
            r.failed = "void";
        }
        // else r.failed == "nan" -- keep it (corruption-class first wins).
    }

    return r;
}

ElasticIncrementResult calc_elasti(const double* deps,
                                   double youngel,
                                   double nuel) {
    // Both D and d_sig zero-initialized; the youngel <= 0 branch leaves
    // them at zero (defensive deviation from Fortran .for:2438-2450).
    ElasticIncrementResult r{};

    if (youngel > 0.0) {
        const double E = youngel;
        const double nu = nuel;
        const double inv_1pnu      = 1.0 / (1.0 + nu);
        const double nu_over_1m2nu = nu / (1.0 - 2.0 * nu);

        // II is the symmetric Voigt identity:
        //   II[i][i] = 1 for i in {0,1,2} (normal slots)
        //   II[i][i] = 1/2 for i in {3,4,5} (shear slots, engineering form)
        //   II[i][j] = 0 for i != j
        // krondelta = (1, 1, 1, 0, 0, 0).
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                double II_ij = 0.0;
                if (i == j) II_ij = (i < 3) ? 1.0 : 0.5;
                const double kron_i = (i < 3) ? 1.0 : 0.0;
                const double kron_j = (j < 3) ? 1.0 : 0.0;
                r.D[i][j] = E * inv_1pnu
                          * (II_ij + nu_over_1m2nu * kron_i * kron_j);
            }
        }

        // d_sig = D * deps (Fortran .for:2447).
        matmul(&r.D[0][0], deps, r.d_sig, 6, 6, 1);
    }

    return r;
}

NormResResult norm_res(const double* y_til, const double* y_hat) {
    NormResResult r{};
    r.norm_R       = 0.0;
    r.nan_detected = false;

    // Stress slice: y[0..5].
    const double* sig_til = &y_til[0];
    const double* sig_hat = &y_hat[0];

    // Internal-state slice: y[6..11] (6 IS components; void is separate).
    // Fortran loops i=1..nasv-1 = 1..6 -> q[1..6]. In C++ that's q[0..5]
    // mapping to y[6..11].
    const double* q_til = &y_til[6];
    const double* q_hat = &y_hat[6];

    // Void: y[12] (Fortran y(6+nasv) = y(13)).
    const double void_til = y_til[12];
    const double void_hat = y_hat[12];

    double del_sig[6];
    for (int i = 0; i < 6; ++i) {
        del_sig[i] = std::abs(sig_hat[i] - sig_til[i]);
    }

    double del_q[6];
    for (int i = 0; i < 6; ++i) {
        del_q[i] = std::abs(q_hat[i] - q_til[i]);
    }

    const double del_void = std::abs(void_hat - void_til);

    // Norms of the "hat" estimate.
    const double norm_sig2 = dot_vect(DotProductKind::StressLike,
                                      sig_hat, sig_hat, 6);
    const double norm_q2   = dot_vect(DotProductKind::StrainLike,
                                      q_hat,   q_hat,   6);
    const double norm_sig  = std::sqrt(norm_sig2);
    const double norm_q    = std::sqrt(norm_q2);

    // Build the relative-error vector (kYDim slots).
    double err[kYDim] = {0.0};

    if (norm_sig > 0.0) {
        for (int i = 0; i < 6; ++i) err[i] = del_sig[i] / norm_sig;
    }
    if (norm_q > 0.0) {
        for (int i = 0; i < 6; ++i) err[6 + i] = del_q[i] / norm_q;
    }

    // Void slot: unguarded division (Fortran: err(6+nasv) = del_void/void_hat).
    err[12] = del_void / void_hat;

    // Plain L2 norm of err.
    const double norm_R2 = dot_vect(DotProductKind::Plain, err, err, kYDim);
    r.norm_R = std::sqrt(norm_R2);

    // NaN/overflow sentinel on norm_sig / norm_q / void_hat (matches
    // Fortran umatisnan_h's finiteness + |v| > 1e30 check; see
    // fails_umatisnan_h in this file).
    if (fails_umatisnan_h(norm_sig)
     || fails_umatisnan_h(norm_q)
     || fails_umatisnan_h(void_hat)) {
        r.norm_R       = 1.0e20;
        r.nan_detected = true;
    }

    // Defensive post-check on the returned residual itself (deviation
    // from Fortran). The Fortran's umatisnan_h scan only inspects the
    // intermediate norms, missing two corruption paths that produce a
    // non-finite norm_R while leaving those intermediates clean:
    //   * `void_hat == 0` (finite) with `del_void > 0` makes
    //     `err[12] = del_void / 0 = +Inf`, poisoning norm_R.
    //   * NaN in y_til (rather than y_hat) makes del_sig or del_q NaN
    //     while norm_sig and norm_q (computed from y_hat) stay clean.
    // Letting either Inf or NaN escape would feed the future RKF
    // driver's step-size control with garbage. Force the sentinel.
    // See Codex review 2026-05-01 (medium finding).
    if (fails_umatisnan_h(r.norm_R)) {
        r.norm_R       = 1.0e20;
        r.nan_detected = true;
    }

    return r;
}

Rkf23Result rkf23_update(const double* y_init,
                         double dtsub_init,
                         double err_tol,
                         int    maxnint,
                         double DTmin,
                         const double* deps_np1,
                         const Parms16& parms,
                         double dtime) {
    Rkf23Result r{};
    r.error          = kErrorOk;
    r.nfev           = 0;
    r.dtsub          = dtsub_init;
    r.failure_class  = FailureClass::Ok;
    // Default y on early-return paths: input y_init (Fortran-faithful;
    // the Fortran does not write y_k back to y on error=3 / error=10).
    for (std::size_t i = 0; i < kYDim; ++i) r.y[i] = y_init[i];

    // ----- Input contract guards (deviation from Fortran) -----
    // The Fortran rkf23_update_h has no input clamp on dtsub or dtime.
    // Its umat caller is responsible for keeping both well-formed across
    // calls. Path B's kernel API is exposed independently, so naive
    // callers can hit each of these silent-failure paths individually.
    // Reject all of them up front with error = 3, no rhs work done.

    // dtime must be finite and strictly positive. Fortran with dtime <= 0
    // can silently complete (DT_k = dtsub_eff / dtime cancels signs and
    // accepts the substep) and return a negative r.dtsub; with dtime == 0
    // it produces NaN substep sizes that fail later after one wasted rhs
    // call. NaN dtime would propagate NaN throughout. Catch all three.
    // See Codex review 2026-05-01 (medium follow-up to the dtsub guard).
    if (!std::isfinite(dtime) || dtime <= 0.0) {
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        return r;
    }

    // dtsub_init <= 0: programmer error (zero/negative time step).
    // Fortran would enter an infinite zero-progress substep loop until
    // maxnint fires; Path B rejects explicitly with error = 3. Also
    // catches NaN dtsub_init since `NaN <= 0.0` is false but `dtsub_init
    // > dtime` would also be false, leaving the NaN to poison DT_k --
    // explicit finiteness check guards that case too.
    if (!std::isfinite(dtsub_init) || dtsub_init <= 0.0) {
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        return r;
    }

    // err_tol must be finite and strictly positive. Invalid err_tol
    // combined with norm_R == 0 (e.g. zero deps) is an especially nasty
    // failure: the strict-< accept check `norm_R < err_tol` is false
    // for any of (0 < 0), (0 < NaN), (0 < negative); the norm_R == 0
    // branch forces S_hull = 1; the reject path computes
    // `DT_k = max(DT_k/4, 1) = 1` -- DT_k never shrinks, T_k never
    // advances, so the loop spins until ksubst > maxnint, wasting up
    // to maxnint*3 rhs evaluations. Reject up front.
    // Codex review 2026-05-01 (medium follow-up to the dtime guard).
    if (!std::isfinite(err_tol) || err_tol <= 0.0) {
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        return r;
    }

    // DTmin must be finite and strictly positive. With DTmin <= 0 or
    // NaN, the reject path's `if (DT_k < DTmin)` exit can never fire
    // -- the loop spins until maxnint with shrunk-but-not-tiny DT_k.
    // Same family of failure as the err_tol gap above; not flagged by
    // Codex 2026-05-01 directly, but caught preemptively here so the
    // input-contract surface is symmetric.
    if (!std::isfinite(DTmin) || DTmin <= 0.0) {
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        return r;
    }

    // dtsub_init > dtime: silent clamp to dtime, mirroring the Fortran
    // accept-branch's own "min(1 - T_k, DT_k)" final-interval style.
    // This makes `r.dtsub` reusable across consecutive calls without
    // requiring an external wrapper to clamp.
    const double dtsub_effective = (dtsub_init > dtime) ? dtime : dtsub_init;

    // Working state. y_k holds the latest accepted state; loop body
    // builds y_2, y_3, y_til, y_hat as scratch.
    double y_k[kYDim];
    for (std::size_t i = 0; i < kYDim; ++i) y_k[i] = y_init[i];

    double T_k  = 0.0;
    double DT_k = dtsub_effective / dtime;   // normalized initial substep, in (0, 1]
    int ksubst  = 0;

    // Map a check_RKF "failed" reason string to FailureClass. Stress-
    // domain labels (pmean/tension) are retryable: a smaller substep may
    // keep the integrated state inside the admissible region next time.
    // Corruption labels (nan/void) are NOT retryable: the substep state
    // hit a numerical fault that smaller dtime cannot fix. Used for
    // every check_RKF gate inside the substep loop. Codex review
    // 2026-05-02 P1 follow-up: closes the rkf23-internal corruption
    // gap so `integrate_step` does not have to demote everything to
    // RkfReject.
    auto rkf_class_from_failed = [](const char* failed) -> FailureClass {
        if (failed == nullptr) return FailureClass::RkfReject;  // defensive default
        if (std::strcmp(failed, "nan")  == 0) return FailureClass::ContractViolation;
        if (std::strcmp(failed, "void") == 0) return FailureClass::ContractViolation;
        // "pmean" / "tension" / unknown future modes: stress-domain.
        return FailureClass::RkfReject;
    };

    while (T_k < 1.0) {
        ksubst += 1;

        // 1. maxnint guard (Fortran .for:1637-1642). Defensive
        //    classification: integrate_step's entry pre-flight already
        //    rejects maxnint <= 0 as ContractViolation, so reaching here
        //    means substep budget legitimately exhausted -- RkfReject.
        if (ksubst > maxnint) {
            r.error          = 3;
            r.failure_class  = FailureClass::RkfReject;
            return r;
        }

        // 2. check_RKF(y_k) before kRK_1. failure class follows the
        //    failed-reason mapping above.
        {
            const auto chk = check_RKF(y_k, parms);
            if (chk.error != 0) {
                r.error          = 3;
                r.failure_class  = rkf_class_from_failed(chk.failed);
                return r;
            }
        }
        // 3. rhs(y_k) -> kRK_1. nfev incremented BEFORE the call so
        //    error=10 from rhs still reflects in r.nfev.
        r.nfev += 1;
        const auto rhs1 = rhs(y_k, deps_np1, parms);
        if (rhs1.error != kErrorOk) {
            r.error          = rhs1.error;   // typically 10
            r.failure_class  = FailureClass::Fatal;
            return r;
        }
        const double* kRK_1 = rhs1.y_dot;

        // 4. y_2 = y_k + (DT_k / 2) * kRK_1.
        double y_2[kYDim];
        const double half_DT_k = 0.5 * DT_k;
        for (std::size_t i = 0; i < kYDim; ++i) {
            y_2[i] = y_k[i] + half_DT_k * kRK_1[i];
        }

        // 5. check_RKF(y_2).
        {
            const auto chk = check_RKF(y_2, parms);
            if (chk.error != 0) {
                r.error          = 3;
                r.failure_class  = rkf_class_from_failed(chk.failed);
                return r;
            }
        }
        // 6. rhs(y_2) -> kRK_2.
        r.nfev += 1;
        const auto rhs2 = rhs(y_2, deps_np1, parms);
        if (rhs2.error != kErrorOk) {
            r.error          = rhs2.error;
            r.failure_class  = FailureClass::Fatal;
            return r;
        }
        const double* kRK_2 = rhs2.y_dot;

        // 7. y_3 = y_k - DT_k*kRK_1 + 2*DT_k*kRK_2.
        double y_3[kYDim];
        for (std::size_t i = 0; i < kYDim; ++i) {
            y_3[i] = y_k[i] - DT_k * kRK_1[i] + 2.0 * DT_k * kRK_2[i];
        }

        // 8. check_RKF(y_3).
        {
            const auto chk = check_RKF(y_3, parms);
            if (chk.error != 0) {
                r.error          = 3;
                r.failure_class  = rkf_class_from_failed(chk.failed);
                return r;
            }
        }
        // 9. rhs(y_3) -> kRK_3.
        r.nfev += 1;
        const auto rhs3 = rhs(y_3, deps_np1, parms);
        if (rhs3.error != kErrorOk) {
            r.error          = rhs3.error;
            r.failure_class  = FailureClass::Fatal;
            return r;
        }
        const double* kRK_3 = rhs3.y_dot;

        // 10. y_til (2nd order) and y_hat (3rd order).
        double y_til[kYDim];
        double y_hat[kYDim];
        const double onesixth  = 1.0 / 6.0;
        const double twothirds = 2.0 / 3.0;
        for (std::size_t i = 0; i < kYDim; ++i) {
            y_til[i] = y_k[i] + DT_k * kRK_2[i];
            y_hat[i] = y_k[i] + DT_k * (onesixth * kRK_1[i]
                                        + twothirds * kRK_2[i]
                                        + onesixth * kRK_3[i]);
        }

        // 11. norm_R = norm_res(y_til, y_hat).
        const auto nres = norm_res(y_til, y_hat);
        const double norm_R = nres.norm_R;

        // 11b. Poisoned-residual fail-closed (Codex 2026-05-02 P2
        //      finding). norm_res sets a corruption sentinel
        //      `norm_R = 1.0e20` AND `nan_detected = true` when any
        //      intermediate norm is NaN/Inf/overflow. Fortran's
        //      rkf23_update_h discards the sentinel signal and only
        //      compares `norm_R < err_tol`; with a buggy caller-
        //      supplied `err_tol > 1.0e20` the corrupt step would be
        //      silently committed. Path B respects the explicit
        //      corruption flag instead -- same logic as check_RKF's
        //      corruption-class label precedence (nan/void override
        //      pmean/tension). Cheap, fail-closed, removes a
        //      garbage-in attack surface that the Fortran-faithful
        //      tolerance comparison cannot catch.
        if (nres.nan_detected) {
            r.error          = 3;
            r.failure_class  = FailureClass::ContractViolation;
            return r;
        }

        // 12. check_RKF(y_hat). Fortran .for:1704: final domain check
        //     before commit; preserved verbatim. Catches the case where
        //     norm_R is acceptably small but y_hat itself has drifted
        //     out of the admissible region.
        {
            const auto chk = check_RKF(y_hat, parms);
            if (chk.error != 0) {
                r.error          = 3;
                r.failure_class  = rkf_class_from_failed(chk.failed);
                return r;
            }
        }

        // 13a. Hull S (Fortran .for:1715-1719).
        double S_hull;
        if (norm_R != 0.0) {
            S_hull = 0.9 * DT_k * std::pow(err_tol / norm_R, 1.0 / 3.0);
        } else {
            S_hull = 1.0;
        }

        // 13b. Accept (norm_R < err_tol, strict) or reject.
        if (norm_R < err_tol) {
            // Commit y_hat -> y_k.
            for (std::size_t i = 0; i < kYDim; ++i) y_k[i] = y_hat[i];
            // Advance T_k, expand DT_k, write back dtsub BEFORE the
            // final-interval clamp. This ordering is locked by the
            // Fortran .for:1730-1733 sequence.
            T_k += DT_k;
            DT_k    = std::min(4.0 * DT_k, S_hull);
            r.dtsub = DT_k * dtime;
            DT_k    = std::min(1.0 - T_k, DT_k);
        } else {
            // Reject: shrink DT_k. Fail if below DTmin.
            DT_k = std::max(DT_k / 4.0, S_hull);
            if (DT_k < DTmin) {
                // Adaptive substepping bottomed out. Genuine RkfReject:
                // the input is physically valid but the outer dtime is
                // too coarse; smaller outer step has a real chance of
                // succeeding.
                r.error          = 3;
                r.failure_class  = FailureClass::RkfReject;
                return r;
            }
        }
    }

    // Loop exited with T_k >= 1.0. Write y_k back into r.y (success).
    for (std::size_t i = 0; i < kYDim; ++i) r.y[i] = y_k[i];
    return r;
}

StepResult integrate_step(const double* sig_n,
                          const double* q_n,
                          const double* deps_np1,
                          const Parms16& parms,
                          double dtime,
                          double dtsub_in,
                          double err_tol,
                          int    maxnint,
                          double DTmin) {
    StepResult r{};
    r.error          = kErrorOk;
    r.nfev           = 0;
    r.used_elastic   = false;
    r.dtsub_next     = dtsub_in;
    r.failure_class  = FailureClass::Ok;
    // r.failure_class flips to RkfReject / ContractViolation / Fatal at
    // the corresponding error exit; success path keeps Ok.

    // Default y/q on any error path: input state. Overwritten only on
    // explicit success or the elastic-fallback path.
    for (int i = 0; i < 6; ++i) r.sig[i] = sig_n[i];
    for (std::size_t i = 0; i < kStateDim; ++i) r.q[i] = q_n[i];
    // r.D is zero-initialized.

    // ---- Pre-flight on scalar tuning parameters --------------------------
    // The elastic-fallback branch below bypasses rkf23_update entirely.
    // Without this guard, an invalid dtime / err_tol / DTmin could leak
    // through the elastic path (which doesn't use them) and corrupt
    // r.dtsub_next via the bottom clamp:
    //   - dtime < 0: any positive dtsub_in satisfies `>= dtime` and
    //     gets pinned to the negative dtime.
    // Reject up front with error = 3 and zero out dtsub_next (the bottom
    // clamp cannot sanitize NaN, and a "default" suggestion of 0 is the
    // right signal for the caller -- "no progress recommendation").
    // Codex review 2026-05-01 (medium follow-up).
    //
    // dtsub_in is the EXCEPTION: it is normalized, not rejected, on
    // non-finite-AND-everything-else values. Fortran umat .for:243-244
    // does:
    //   if ((dtsub.le.0.0d0) .or. (dtsub.gt.dtime)) dtsub = dtime
    // -- treating dtsub <= 0 as the legitimate "fresh statev" first-call
    // signal (Kratos / Abaqus default-initializes statev to zero, so
    // every first call legally arrives with `dtsub_in == 0`). Rejecting
    // it would deadlock the bridge: the caller cuts back, reads zero
    // from statev again, and never makes progress (Codex review
    // 2026-05-02 high finding). NaN/Inf are still hard rejects -- they
    // cannot be normalized and almost certainly indicate a real fault
    // upstream.
    if (!std::isfinite(dtime)   || dtime   <= 0.0
     || !std::isfinite(err_tol) || err_tol <= 0.0
     || !std::isfinite(DTmin)   || DTmin   <= 0.0
     || !std::isfinite(dtsub_in)
     || maxnint <= 0) {
        // maxnint <= 0 is structurally a contract violation: a
        // smaller outer dtime cannot fix an invalid substep limit.
        // Without this guard, rkf23_update would return error == 3
        // immediately on the first ksubst > maxnint check (cpp:1040)
        // and integrate_step would classify it as RkfReject -- a
        // bridge that retried with cutback would loop forever.
        // Codex review 2026-05-02 (P2 finding: invalid maxnint
        // misclassified as retryable cutback). Logical closure of
        // the FailureClass introduction (Codex 2026-05-02 P1).
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        r.dtsub_next     = 0.0;   // explicit; bottom clamp would not handle NaN.
        return r;
    }
    // Mirror Fortran .for:243-244 normalization. After this, dtsub_in
    // is guaranteed in (0, dtime].
    if (dtsub_in <= 0.0 || dtsub_in > dtime) {
        dtsub_in     = dtime;
        r.dtsub_next = dtime;   // keep the success-path default in sync.
    }

    // Pack y_init for check_RKF.
    double y_init[kYDim];
    for (int i = 0; i < 6; ++i) y_init[i] = sig_n[i];
    for (std::size_t i = 0; i < kStateDim; ++i) y_init[6 + i] = q_n[i];

    // ---- Pre-flight input contract on y_init -------------------------
    // Inf/NaN/overflow values can DISGUISE themselves as `pmean` or
    // `tension` failures in check_RKF: e.g. sig[0] = +Inf gives
    // pmean = -Inf, which compares <= minstress and trips the pmean
    // branch -- the umatisnan_h scan that would have flagged "nan"
    // never gets to set the failed string. If integrate_step then
    // accepted "pmean" as inittension, corrupted state would silently
    // pass through the elastic fallback. Decouple from check_RKF's
    // internal failure-priority ordering by validating finiteness +
    // overflow + void here, before asking check_RKF for an opinion.
    // Codex review 2026-05-01 (high finding follow-up).
    auto value_is_corrupt = [](double v) {
        return !std::isfinite(v) || std::abs(v) > 1.0e30;
    };
    bool input_corrupt = false;
    for (std::size_t i = 0; i < kYDim; ++i) {
        if (value_is_corrupt(y_init[i])) { input_corrupt = true; break; }
    }
    if (!input_corrupt && y_init[12] <= 0.0) input_corrupt = true;
    // deps_np1 must also be finite. The elastic-fallback branch passes
    // deps_np1 directly to `calc_elasti`, which has no error channel
    // and matmul-multiplies it into d_sig; a NaN/Inf/overflow deps
    // would silently corrupt r.sig with `error == 0` and
    // `used_elastic == true`. The RKF main path catches deps NaN
    // via the chain reaction through rhs / get_tan / check_RKF, but
    // the elastic short-path bypasses those gates, so we have to
    // close the gap here. Codex review 2026-05-02 (high finding).
    if (!input_corrupt) {
        for (int i = 0; i < 6; ++i) {
            if (value_is_corrupt(deps_np1[i])) { input_corrupt = true; break; }
        }
    }
    // Norm-level overflow guard. Mirrors Fortran umat .for:201-209,
    // which computes norm_D = sqrt(dot_vect_h(2, deps_np1, deps_np1, 6))
    // and runs umatisnan_h on the SCALAR norm before any branch.
    // Per-component finiteness is necessary but not sufficient: e.g.
    // two normal components at 8e29 each pass the per-component
    // |v| <= 1e30 check, but ||deps||_strain ~ 1.13e30 still trips
    // Fortran's overflow sentinel. Without this guard, a tensile
    // initial sig would route to elastic fallback and `calc_elasti`
    // would matmul those huge-but-finite components into d_sig
    // (~1e32 magnitude), returning success with overflow-scale r.sig.
    // The kernel's contract is to match the Fortran reference; leaving
    // this gap open contradicts that. Codex review 2026-05-02
    // (high finding follow-up).
    if (!input_corrupt) {
        const double norm_D2 = dot_vect(DotProductKind::StrainLike,
                                        deps_np1, deps_np1, 6);
        const double norm_D  = std::sqrt(norm_D2);
        if (fails_umatisnan_h(norm_D)) input_corrupt = true;
    }

    if (input_corrupt) {
        // r.sig / r.q already hold sig_n / q_n (caller-side corrupted
        // inputs propagate back unchanged -- we don't sanitize). r.D
        // stays zero; r.nfev = 0; r.dtsub_next stays at dtsub_in and
        // gets clamped at the end of this routine.
        r.error          = 3;
        r.failure_class  = FailureClass::ContractViolation;
        r.used_elastic   = false;
    } else {

    // ---- Initial RKF gate: detect inittension ----
    const auto chk = check_RKF(y_init, parms);

    if (chk.error != 0) {
        // After pre-flight, the only failure modes check_RKF can fire
        // are "pmean" or "tension" -- the umatisnan_h and void rules
        // are ruled out by the pre-flight. The defensive failed-name
        // split below preserves correctness if a future check_RKF
        // adds new failure modes that the pre-flight does not cover.
        const bool is_inittension = chk.failed != nullptr
            && (std::strcmp(chk.failed, "pmean")   == 0
             || std::strcmp(chk.failed, "tension") == 0);

        if (is_inittension) {
            // ===== Elastic fallback (Fortran inittension) =====
            // Hardcoded constants from .for:329 area: youngel = 100, nuel = 0.48.
            const double youngel = 100.0;
            const double nuel    = 0.48;
            const auto el = calc_elasti(deps_np1, youngel, nuel);
            for (int i = 0; i < 6; ++i) r.sig[i] = sig_n[i] + el.d_sig[i];
            // r.q stays at q_n (do not advance hypoplastic state in the
            // elastic regime; this matches the Fortran behavior where
            // statev(1..7) are not updated under inittension).
            for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                    r.D[i][j] = el.D[i][j];
                }
            }
            r.error        = kErrorOk;
            r.used_elastic = true;
            r.nfev         = 0;
            // r.dtsub_next stays at dtsub_in -- will be clamped at end.
        } else {
            // Numerical / void / other contract violation. Fail closed.
            // r.sig and r.q already at sig_n / q_n (set above before the
            // initial gate). r.D stays zero. r.nfev stays at 0.
            // r.dtsub_next stays at dtsub_in -- will be clamped at end.
            r.error          = 3;
            r.failure_class  = FailureClass::ContractViolation;
            r.used_elastic   = false;
        }
    } else {
        // ===== Normal RKF path or zero-deps short circuit =====
        //
        // Zero-strain-rate short circuit (Fortran umat .for:254
        // `if(norm_D.eq.0) testing=2`, then .for:273-276 just copies
        // y_n into y without calling rkf23_update_h). Skipping the
        // integrator avoids three unnecessary rhs evaluations and
        // keeps the persisted dtsub untouched on tangent-only queries
        // (e.g. Kratos `CalculateMaterialResponse` with zero strain
        // increment to fetch stiffness). Fortran's perturbate_h is
        // still invoked afterwards (.for:328) to compute the final
        // tangent; we mirror that with the shared get_tan call below.
        // Codex review 2026-05-02 (P2 finding: zero-deps testing=2).
        const double norm_D2_check = dot_vect(DotProductKind::StrainLike,
                                              deps_np1, deps_np1, 6);
        const bool zero_deps = (norm_D2_check == 0.0);

        if (zero_deps) {
            // y unchanged: r.sig and r.q already hold sig_n / q_n.
            r.error = kErrorOk;
            r.nfev  = 0;
            // r.dtsub_next stays at the normalized dtsub_in (set at
            // the top of this routine). Fortran .for:343-348 also
            // leaves dtsub untouched on the testing=2 path -- the
            // bottom clamp at the end of integrate_step keeps it
            // semantically valid for persistence.
        } else {
            const auto rkf = rkf23_update(y_init, dtsub_in, err_tol,
                                          maxnint, DTmin, deps_np1, parms,
                                          dtime);
            r.error      = rkf.error;
            r.nfev       = rkf.nfev;
            r.dtsub_next = rkf.dtsub;

            if (rkf.error == kErrorOk) {
                // Commit y_final -> sig and q.
                for (int i = 0; i < 6; ++i) r.sig[i] = rkf.y[i];
                for (std::size_t i = 0; i < kStateDim; ++i) r.q[i] = rkf.y[6 + i];
            } else {
                // Path B closure (Codex 2026-05-02 P1 follow-up):
                // rkf23_update now sets failure_class precisely at every
                // error exit (mid-substep check_RKF labels mapped via
                // corruption-class precedence, nan_detected ->
                // ContractViolation, DTmin -> RkfReject, rhs error=10
                // -> Fatal, scalar contract guards -> ContractViolation).
                // Pass it through verbatim. No more demotion to
                // RkfReject for rkf23-internal corruption.
                r.failure_class = rkf.failure_class;
            }
        }

        // Shared final-tangent computation. Runs whenever the path so
        // far succeeded (zero-deps short circuit OR successful RKF).
        // The Fortran reference is perturbate_h (.for:1424-1493),
        // which extracts sig/q from y_n (NOT y_np1; see .for:1466-1469)
        // before calling get_tan_h. Tangent at y_n is the "consistent
        // linearization" around the start of the increment -- it is
        // what an outer-loop Newton iteration uses to predict the next
        // trial state. Tangent at y_np1 would be self-consistent with
        // the converged stress but is NOT what the Fortran returns to
        // its caller, so Path A and Path B would diverge here.
        //
        // Fortran perturbate_h selects:
        //   istrain = 0 -> D = L          (.for:1478-1483)
        //   istrain = 1 -> D = m_R * L    (.for:1484-1490)
        // (NOT the IS-modified M.)
        // Codex review 2026-05-02 (high finding).
        if (r.error == kErrorOk) {
            const int istrain = (parms[9] <= 0.5) ? 0 : 1;
            const auto tan_final = get_tan(deps_np1, sig_n, q_n,
                                           parms, istrain);
            if (tan_final.error != kErrorOk) {
                // Tangent failed at the begin-of-step state. Promote
                // to the same fatal-style passthrough as rkf23
                // error=10: revert y to input and zero out D so the
                // caller cannot act on a half-baked result.
                //
                // Reachable through validated parms via the fb_temp1
                // < 0 path (parms-only failure mode that check_parms
                // does not validate). The rkf23-side error=10
                // rollback shares this output contract; both are
                // covered directly by test_integrate_step case 13
                // (Codex 2026-05-02 P2 finding -- error=10 rollback
                // had no direct test through validated parms).
                r.error = tan_final.error;   // kErrorFatal (10)
                r.failure_class = FailureClass::Fatal;
                for (int i = 0; i < 6; ++i) r.sig[i] = sig_n[i];
                for (std::size_t i = 0; i < kStateDim; ++i) r.q[i] = q_n[i];
                // r.D stays zero (already initialized).
            } else {
                const double m_R   = parms[9];
                const double scale = (istrain == 0) ? 1.0 : m_R;
                for (int i = 0; i < 6; ++i) {
                    for (int j = 0; j < 6; ++j) {
                        r.D[i][j] = scale * tan_final.L[i][j];
                    }
                }
            }
        }
    }
    }  // end "if (input_corrupt) { ... } else { ... }"

    // ---- dtsub write-back clamp (output-side guard) ----
    // The kernel rkf23_update's `r.dtsub` can exceed dtime as an
    // adaptive "try larger next step" suggestion. A step-level caller
    // persists dtsub across invocations, so clamp the persisted value
    // into the safe range here. Negative or zero is mapped to 0
    // (caller will treat as "no progress recommendation").
    if (r.dtsub_next <= 0.0) {
        r.dtsub_next = 0.0;
    } else if (r.dtsub_next >= dtime) {
        r.dtsub_next = dtime;
    }

    return r;
}

}  // namespace SandHypoCpp
}  // namespace Kratos

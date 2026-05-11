// Path B: pure C++ inline port of the sand-hypoplastic Fortran kernel
// (`00_hypo_is_original.for`).
//
// This header declares the path-private kernel namespace
// `Kratos::SandHypoCpp`. Each translated routine mirrors a Fortran
// subroutine of the same root name (`check_parms_h` -> `check_parms`,
// `dot_vect_h` -> `dot_vect`, `matmul_h` -> `matmul`).
//
// Path B must not include or call into Path A or Path C code.
//
// =====================================================================
// Voigt convention (single source of truth for the kernel)
// =====================================================================
//
// All 6-component "Voigt" vectors and 6x6 tangent matrices INSIDE this
// kernel use the Fortran/Abaqus order from `00_hypo_is_original.for`:
//
//   1-based Fortran slot   |  0-based C++ index  |  tensor component
//   ---------------------- | ------------------- | -----------------
//                       1  |                  0  |         (1,1)
//                       2  |                  1  |         (2,2)
//                       3  |                  2  |         (3,3)
//                       4  |                  3  |         (1,2) = (2,1)
//                       5  |                  4  |         (1,3) = (3,1)
//                       6  |                  5  |         (2,3) = (3,2)
//
// Equivalently, the 0-based C++ layout is `[11, 22, 33, 12, 13, 23]`.
//
// The Kratos convention is `[11, 22, 33, 12, 23, 13]`, which differs
// from this kernel by swapping C++ indices 4 and 5 (i.e., the 1-based
// Fortran slots 5 and 6). The Path B bridge layer
// (`SandHypoplasticCppLaw`, when written) is responsible for swapping
// C++ indices 4 <-> 5 on:
//   - input stress (`sig_n`) coming from Kratos (e.g.
//     `STRESS_VECTOR`, `INITIAL_STRESS_VECTOR`) before passing it
//     to `integrate_step`,
//   - input strain increments coming from Kratos,
//   - output stress and tangent vectors going back to Kratos,
//   - rows AND columns of the 6x6 tangent matrix going back to Kratos,
//   - the intergranular-strain state when written to the shared CSV
//     in Kratos order.
//
// "Swap C++ indices 4 and 5", not "swap 5 and 6": there is no index 6
// in a 6-element array, so 0-based vs 1-based confusion is the most
// common bug here. Inside this kernel namespace, every routine assumes
// the Abaqus 0-based layout above and never sees the Kratos layout.
//
// State vector (`q_n` / `StepResult::q`) Voigt order is also Abaqus.
// Slots q[0..5] are the intergranular-strain components `del`; slots
// 0..2 are normal, slots 3..5 are shear in `[del_12, del_13, del_23]`
// order. Slot q[6] is the scalar void ratio (no Voigt issue).
//
// Runtime persistence rule (binding on the bridge): Path B's persisted
// statev MUST be stored in kernel/Abaqus order between calls. Bridges
// that store statev in Kratos order would have to swap q[4] <-> q[5]
// on every entry to AND every exit from `integrate_step`; missing
// either swap silently feeds wrong IS components into the next
// step's get_tan and accumulates undetected error. The Kratos slot-4/5
// swap therefore happens ONLY at Kratos exposure boundaries (e.g.
// CalculateValue / SetValue overrides for IS-related variables, CSV
// writers, MaterialResponse setters/getters), never on the persisted
// statev itself. The kernel's Abaqus order is the storage canonical.

#pragma once

#include <array>
#include <cstddef>

namespace Kratos {
namespace SandHypoCpp {

// =====================================================================
// Material parameters: `check_parms`
// =====================================================================

constexpr std::size_t kNumProps = 16;
// State vector dimension matching the Fortran `nasv = 7`:
// 6 intergranular-strain components (`del[1..6]`) plus 1 void ratio.
constexpr std::size_t kStateDim = 7;

using Props16 = std::array<double, kNumProps>;
using Parms16 = std::array<double, kNumProps>;

// Mirrors the Fortran `error` codes: 0 OK, 10 fatal.
constexpr int kErrorOk = 0;
constexpr int kErrorFatal = 10;

// Result of `check_parms`. On success (`error == kErrorOk`), `parms`
// holds the kernel-internal parameter array (with `phi` converted from
// degrees to radians); `failed` is null. On failure, `parms` is left in
// the partial state at which validation aborted, `failed` names the
// offending parameter, and `bad_value` holds its rejected value.
struct CheckParmsResult {
    Parms16 parms;
    int error;
    const char* failed;
    double bad_value;
};

// Validates the 16-slot Fortran-ordered material parameter vector and
// produces the kernel-internal `parms` array.
//
// Index layout (1-based Fortran -> 0-based C++ shown for clarity):
//   props(1)  = parms[0]  phi (deg on input, radians on output)
//   props(2)  = parms[1]  p_t
//   props(3)  = parms[2]  hs
//   props(4)  = parms[3]  en
//   props(5)  = parms[4]  ed0
//   props(6)  = parms[5]  ec0
//   props(7)  = parms[6]  ei0
//   props(8)  = parms[7]  alpha
//   props(9)  = parms[8]  beta
//   props(10) = parms[9]  m_R
//   props(11) = parms[10] m_T
//   props(12) = parms[11] r_uc
//   props(13) = parms[12] beta_r
//   props(14) = parms[13] chi
//   props(15) = parms[14] bulk_w
//   props(16) = parms[15] e_0 (raw encoding; not validated here)
//
// Note: this routine mirrors the Fortran's IEEE-blind comparisons. NaN
// or +/-Inf inputs will pass validation. Callers (in particular the
// Path B host wrapper's `Check()`) MUST guarantee finiteness before
// passing `props` here.
CheckParmsResult check_parms(const Props16& props);

// =====================================================================
// Linear-algebra helpers: `dot_vect`, `matmul`
// =====================================================================

// Three contraction modes corresponding to the `flag` argument of
// `dot_vect_h` (.for:532-579). `StressLike` and `StrainLike` apply
// Voigt-tensor contraction weights to off-diagonal slots (C++ indices
// >= 3); `Plain` is an ordinary R^n dot product with no weighting.
//
// The numeric values of these enumerators are deliberately unspecified:
// callers must use the named constants. Do not `static_cast` an integer
// (least of all a Fortran flag value) into this enum.
enum class DotProductKind {
    StressLike,  // off-diagonal weight = 2
    StrainLike,  // off-diagonal weight = 1/2
    Plain        // off-diagonal weight = 1
};

// Voigt-aware dot product.
//
// Returns sum over i in [0, n):
//     a[i]*b[i]                              if i <  3
//     w * a[i]*b[i]   where w = weight(kind) if i >= 3
//
// For n <= 3 the kind has no effect (no off-diagonal slots are
// reached). For Voigt-stress / Voigt-strain semantics callers must
// pass n = 6 with the Abaqus-ordered layout documented at the top
// of this file; the routine itself does not depend on which physical
// pair (1,3) vs (2,3) sits in slot 4 vs slot 5, only on i < 3 vs i >= 3.
double dot_vect(DotProductKind kind,
                const double* a,
                const double* b,
                std::size_t n);

// Row-major matrix multiplication: c (l x n) = a (l x m) * b (m x n).
//
// Indexing:
//   a[i, k] = a[i*m + k]
//   b[k, j] = b[k*n + j]
//   c[i, j] = c[i*n + j]
//
// Aliasing: `c` MUST NOT alias `a` or `b` (whether equal-pointer or
// partial overlap). The implementation writes c[i,j] in place after
// reading all required a/b entries, so aliasing yields undefined
// values. No runtime alias check is performed.
//
// This is the literal port of `matmul_h` (.for:1233-1253). It will
// likely be replaced by Eigen at the call sites once `get_tan_h` is
// translated and the kernel becomes Kratos/Eigen-aware; until then,
// keep raw-pointer flat arrays so the line-by-line correspondence to
// the Fortran source is preserved.
void matmul(const double* a,
            const double* b,
            double* c,
            std::size_t l,
            std::size_t m,
            std::size_t n);

// =====================================================================
// Strain and stress invariants: `inv_eps`, `inv_sig`
// =====================================================================
//
// Both routines take a 6-component Voigt vector in the kernel-internal
// Abaqus order (see Voigt convention at the top of this file) and
// return scalar invariants. They mirror `inv_eps_h` (.for:1026-1104) and
// `inv_sig_h` (.for:1106-1231) verbatim, including the latter's I1==0
// numerical hack (divide by tiny=1e-18) and the Lode-angle clamps.

struct StrainInvariants {
    double eps_v;   // volumetric strain: trace(eps) = eps_11 + eps_22 + eps_33
    double eps_s;   // equivalent shear strain: sqrt(2/3 * ||edev||^2)
    double sin3t;   // Lode angle sin(3*theta), in [-1, 1] (post-clamp)
};

// Strain invariants. `eps` must be a length-6 Voigt-strain (Abaqus
// order, off-diagonals are engineering shear strain gamma_ij = 2*eps_ij).
//
// Edge case (mirrors Fortran early return): if the equivalent shear
// strain `eps_s` is exactly zero (purely hydrostatic strain), the Lode
// angle is undefined and the routine returns `sin3t = -1`.
StrainInvariants inv_eps(const double* eps);

struct StressInvariants {
    double pp;      // mean stress: I1/3
    double qq;      // von Mises equivalent: sqrt(3/2 * ||sdev||^2)
    double cos3t;   // Lode angle cos(3*theta), in [-1, 1] (post-clamp)
    double I1;      // first principal invariant: trace(sig)
    // Second invariant in the Fortran convention: I2 = 0.5*(sig:sig - I1^2).
    // This is the OPPOSITE sign of the characteristic-polynomial convention
    // (`det(sig - lambda I) = -lambda^3 + I1 lambda^2 - I2 lambda + I3`)
    // common in many textbooks. Hydrostatic compression therefore yields
    // *negative* I2 here. Do not "correct" the sign without matching
    // changes in Fortran and Path A.
    double I2;
    double I3;      // third principal invariant: det(sig)
};

// Stress invariants. `sig` must be a length-6 Voigt-stress (Abaqus
// order, off-diagonals are tensor shear sig_ij directly).
//
// Numerical caveat preserved verbatim from Fortran (`inv_sig_h`):
// when `I1 == 0`, the normalized stress eta=sig/I1 is replaced by
// eta=sig/tiny with `tiny = 1.0e-18`. This hack lets the eta-based
// Lode-angle formula run without dividing by zero. The cos3t result
// in this branch is dominated by the tiny-divide rescaling and the
// post-clamp behavior; it should be considered a numerical placeholder,
// not a physically accurate Lode angle. Do not "fix" this without
// matching changes in Path A and the Fortran reference.
StressInvariants inv_sig(const double* sig);

// =====================================================================
// Constitutive tangent: `get_tan`
// =====================================================================
//
// Computes the von Wolffersdorff sand-hypoplastic + Niemunis-Herle
// intergranular-strain tangent at the current step. Mirrors `get_tan_h`
// in `00_hypo_is_original.for` (.for:631-997) verbatim. This is the
// physics core of Path B; it is intentionally a literal port and will
// be optimized only after Path A and Path B agree on numerical outputs.

struct TangentResult {
    // Full tangent stiffness matrix with intergranular-strain
    // modifications applied. Used by callers when istrain == 1.
    // When istrain == 0 the Fortran leaves M zero and the caller is
    // expected to use the post-scaling L instead.
    double M[6][6];

    // Hardening matrix (state evolution operator). Rows 0..5 are the
    // intergranular-strain evolution H_del; row 6 is the void-ratio
    // evolution H_e. Layout matches Fortran HH(nasv,6) with nasv = 7.
    double H[kStateDim][6];

    // Basic linear tangent L (in the Fortran sense), POST fs scaling.
    double L[6][6];

    // Basic non-linear N vector, POST fs*fd scaling.
    double N[6];

    // Mirrors Fortran error: 0 OK, 10 fatal. Two fatal paths exist:
    //   - "factor fb not defined" guard (.for:839, fb_temp1 < 0):
    //     fires BEFORE L/N are computed, so M/H/L/N are all zero on
    //     this path.
    //   - r_uc <= 0 belt-and-suspenders guard (Path B deviation,
    //     .for-equivalent absent): fires AFTER L/N are computed,
    //     so L and N may carry partial values; only M and H stay
    //     at the zero-init default.
    // For the routine-wide failure-output contract see the docstring
    // above `get_tan` -- on any non-zero `error`, the caller MUST
    // check `error` first and treat all output fields as undefined.
    // Don't depend on either zero pattern for routing.
    int error;
};

// Compute the constitutive tangent.
//
// Inputs:
//   deps    [6]  strain-rate vector (Voigt-strain, Abaqus order;
//                off-diagonals are engineering shear)
//   sig     [6]  current effective Cauchy stress (Voigt-stress)
//   q       [7]  internal state: q[0..5] = del (intergranular strain),
//                                q[6]    = void ratio
//   parms        validated material parameters from check_parms
//   istrain      0 = no-IS branch (HH carries only void evolution; M
//                    is zero, caller must use L), 1 = IS branch
//                    (full HH and M)
//
// The selection of istrain in the Fortran callers is governed by the
// rule `parms(10) <= 0.5 -> istrain = 0` (parms(10) == m_R). That
// decision is the caller's job in this port too; get_tan trusts the
// flag it is given.
//
// Caveats preserved verbatim from Fortran:
//   - Tensile/zero mean stress is NOT guarded here; the caller is
//     expected to short-circuit via the inittension/calc_elasti
//     fallback (see `00_hypo_is_original.for:766+` and the project
//     UMAT-call-graph memory). With I1 >= 0 this routine produces
//     NaN via `pow(-I1/hs, en)`.
//   - The dependence on the Lode angle `cos3t` is taken from
//     inv_sig(sig_star); when sig_star is hydrostatic the cos3t early
//     return path is used (see inv_sig docstring).
//
// Failure-output contract: on any non-zero error code, all output
// matrices/vectors (M, H, L, N) are UNDEFINED -- the caller MUST check
// `error` first and ignore every output field on failure. Different
// failure paths leave different intermediate values:
//   - The fb_temp1 < 0 guard fires before L and N are computed, so
//     those happen to be zero on that path.
//   - The r_uc guard below fires AFTER L and N have been computed for
//     the no-IS branch, so L and N may carry partial-tangent values
//     even though M and H stay at the zero-init default.
// Don't depend on either pattern. test_get_tan case 5's M/H-only
// assertion accommodates this; if you write new failure-path tests,
// match that scope.
//
// Deliberate deviation from Fortran (belt-and-suspenders after the
// check_parms hoist):
//   - When istrain == 1 AND r_uc <= 0, this routine returns
//     kErrorFatal. The Fortran original
//     unconditionally computes rho = norm_del / r_uc, which is NaN
//     under the common zero-IS initial state when r_uc == 0.
//   - As of Codex review 2026-05-02 (medium follow-up: path-
//     dependence), check_parms now rejects the cross-parameter combo
//     `m_R > 0.5 && r_uc <= 0` at validation time (failed = "r_uc"),
//     so this in-kernel guard is structurally unreachable through
//     `Parms16` produced by check_parms.
//   - The guard is retained as defense-in-depth for callers that
//     construct a `Parms16` without running check_parms (e.g.
//     deliberately invalid parms in lower-layer tests, memcpy from
//     a serialized blob, future test seams). See Codex review
//     2026-05-01 (high-severity finding) for the original rationale.
TangentResult get_tan(const double* deps,
                      const double* sig,
                      const double* q,
                      const Parms16& parms,
                      int istrain);

// =====================================================================
// Right-hand-side rates: `get_F_sig_q`
// =====================================================================
//
// Mirrors `get_F_sig_q_h` in `00_hypo_is_original.for` (.for:581-629).
// Composes get_tan + matmul to produce the time-rate vectors F_sig
// (stress rate) and F_q (state rate) for the ODE integrator.
//
// istrain selection: the Fortran callers (`get_F_sig_q_h`,
// `rkf23_update_h`, `calc_statev_h`) uniformly use `parms(10) <= 0.5`
// (m_R <= 0.5) to choose istrain. This routine bakes that rule in;
// callers do NOT pass istrain.
//
// Branch structure for F_sig:
//   istrain == 0 (no-IS): F_sig = L * deps + N * |deps|_strain
//                         where |deps|_strain^2 = StrainLike-dot(deps, deps)
//   istrain == 1 (IS):    F_sig = M * deps
//
// F_q is uniformly H * deps in both branches. Note that in the no-IS
// branch H[0..5] are zero and only the void-ratio row (H[6]) carries
// non-zero entries; F_q[0..5] are therefore zero in that branch.
//
// Error transparency: when get_tan returns a non-zero error code
// (today: r_uc <= 0 in the IS branch, or fb_temp1 < 0 anywhere), this
// routine propagates it verbatim and leaves F_sig and F_q at their
// zero-initialized state. No partial/garbage rates ever escape on a
// fatal-error path.
//
// Precondition (caller's responsibility, NOT enforced here):
//   sig must be in the admissible domain for get_tan, i.e. the shifted
//   trace must be compressive: sig[0] + sig[1] + sig[2] < 3 * p_t.
//   Outside that domain get_tan inherits NaN through `pow(-I1/hs, en)`
//   with negative base (.for:833 / kernel.cpp), and this routine
//   currently propagates the NaN into F_sig and F_q while reporting
//   error == 0. The Fortran architecture handles this via a
//   `check_RKF_h` gate at the RKF substep boundary (.for:1768-1831),
//   which Path B will translate in a later round; once that gate is
//   in place, this contract becomes structurally enforceable.
//   See Codex review 2026-05-01 (medium finding, deferred).

struct RhsResult {
    double F_sig[6];
    double F_q[kStateDim];
    int error;
};

RhsResult get_F_sig_q(const double* deps,
                      const double* sig,
                      const double* q,
                      const Parms16& parms);

// =====================================================================
// ODE right-hand side: `rhs`
// =====================================================================
//
// Mirrors `rhs_h` in `00_hypo_is_original.for` (.for:1511-1567). Wraps
// get_F_sig_q with the y -> (sig, q) unpacking and the F_sig/F_q ->
// y_dot packing expected by an RKF integrator.
//
// State vector layout (Fortran y(1..13) -> C++ y[0..12]):
//   y[0..5]   stress sigma (Voigt-stress, Abaqus order)
//   y[6..11]  intergranular strain del (state q[0..5])
//   y[12]     void ratio (state q[6])
//
// Deviation from Fortran (single, narrow): the Fortran takes nfev
// in/out and increments it inside; this port leaves the
// function-evaluation counter to the future RKF driver. `rhs` itself
// is purely functional in its result.
//
// Error transparency: any non-zero error from get_F_sig_q (today: the
// r_uc <= 0 IS-branch guard or the fb_temp1 < 0 guard) propagates
// verbatim and leaves y_dot at its zero-initialized state. No partial
// rate ever escapes on the fatal path.

constexpr std::size_t kYDim = 6 + kStateDim;  // 13: stress + state

struct YdotResult {
    double y_dot[kYDim];
    int error;
};

YdotResult rhs(const double* y,
               const double* deps,
               const Parms16& parms);

// =====================================================================
// Principal stresses (3x3 symmetric eigenvalues, closed-form)
// =====================================================================
//
// Closed-form Cardano solution for the eigenvalues of a 3x3 symmetric
// tensor in Abaqus-Voigt order. Replaces the Fortran's `Eig_3a_h`
// (Jacobi rotation iteration). The closed form is well-conditioned
// for the stress regimes we exercise (compressive, mildly anisotropic).
//
// Voigt convention: this routine takes an Abaqus-Voigt stress vector
// directly. The Fortran source's `check_RKF_h` does a slot 5<->6 swap
// before calling PrnSig_h because PrnSig uses a different convention;
// in Path B we drop that swap entirely and keep the kernel consistent
// with the file-level "all internal arrays are Abaqus-Voigt" rule.

struct PrincipalStresses {
    // Sorted ascending: s[0] <= s[1] <= s[2].
    // In tension-positive convention, s[2] is the most-tensile principal.
    double s[3];
};

PrincipalStresses principal_stresses_3(const double* sig);

// =====================================================================
// RKF substep validity gate: `check_RKF`
// =====================================================================
//
// Mirrors `check_RKF_h` in `00_hypo_is_original.for` (.for:1768-1831).
// Validates that the current RKF substep state y is admissible for
// the constitutive routines (in particular get_tan, which produces
// NaN outside its compressive-stress domain).
//
// Fails (error = 1) if any of:
//   - shifted mean stress (compression-positive) <= p_t / 4,
//   - the most-tensile principal of sig_star (compression-positive)
//     <= p_t / 4 (i.e., any principal stress is insufficiently
//     compressive after the cohesion shift),
//   - any of y[0..12] or computed tmin/min_principal fails the Fortran
//     `umatisnan_h` predicate: NaN, +/-Inf, or outside +/-1.0e30,
//   - the void ratio y[12] is finite but <= 0 (deliberate deviation
//     from Fortran). The Fortran's check_RKF_h does not enforce
//     void > 0, but get_tan's `fe = pow(ec / void_ratio, beta)`
//     produces NaN/Inf when void_ratio <= 0; the silent-corruption
//     path is an accepted y_hat that has admissible stress but
//     non-positive void after a large volumetric-compression
//     substep. NaN/Inf void is already caught by the umatisnan_h
//     scan; this rule handles the finite-but-non-positive case.
//     See Codex review 2026-05-01 (medium follow-up).
//
// On a successful call (error = 0) the diagnostic fields `pmean` and
// `min_principal` carry the computed values, useful for post-hoc
// inspection. On the fatal path, `failed` names the most-severe
// check that fired, with class precedence (Codex 2026-05-02 P1):
//   - CORRUPTION class:   "nan", "void"
//   - STRESS-DOMAIN class: "pmean", "tension"
// Corruption-class labels UNCONDITIONALLY override stress-domain
// labels, so that `integrate_step` (the only consumer of the label)
// can route stress-domain failures to the elastic inittension
// fallback while still failing-closed on numeric corruption.
//
// Within each class first-failure-wins still applies: "pmean" beats
// "tension"; an early "nan" from the raw-y umatisnan_h scan is kept
// even if a later "nan" from the min_principal overflow check would
// have set the same label. Across classes, however, even a late
// "nan" / "void" overwrites an earlier "pmean" / "tension".
//
// The "nan" label is retained for all `umatisnan_h`-style numeric
// failures, including infinities and finite overflow sentinels --
// including the DERIVED `min_principal` overflow path that is the
// motivating example for the precedence rule (a sig with all
// shears at +/-1e30 lifts principal eigenvalues to ~2e30, which
// trips the tension threshold first but is morphologically a
// numeric overflow and must not be silently absorbed as elastic).
//
// This is the gate referenced by the deferred Codex finding in the
// `get_F_sig_q` docstring: when an RKF driver invokes check_RKF
// before each rhs evaluation, the get_F_sig_q precondition becomes
// structurally enforced.

struct CheckRkfResult {
    int error;             // 0 OK, 1 fatal (matches Fortran error_RKF).
    const char* failed;    // First failure: "pmean", "tension", "nan", or nullptr.
    double pmean;          // Computed mean stress (compression-positive) of sig_star.
    double min_principal;  // -max(principal stresses) of sig_star (compression-positive).
};

CheckRkfResult check_RKF(const double* y, const Parms16& parms);

// =====================================================================
// Linear-elastic fallback: `calc_elasti`
// =====================================================================
//
// Mirrors `calc_elasti_h` in `00_hypo_is_original.for` (.for:2372-2453).
// Used in the umat driver's `inittension` short-circuit path: when the
// incoming stress state is tensile/invalid for hypoplasticity, the
// driver substitutes a linear-elastic response (E = 100 kPa, nu = 0.48
// in the original) instead of integrating the rate equations.
//
// Builds the symmetric Voigt elastic stiffness D from young's modulus
// E and Poisson's ratio nu using the standard
//   D_ijkl = (E / (1 + nu)) * (II_ijkl + (nu / (1 - 2 nu)) * delta_ij delta_kl)
// formulation, where II is the symmetric identity (1 on the normal
// diagonal, 1/2 on the shear diagonal — this is the engineering-strain
// form, so deps off-diagonals carry gamma_ij = 2 eps_ij). Returns
// d_sig = D * deps as the stress increment for the elastic step.
//
// Defensive deviation from Fortran: when `youngel <= 0`, both D and
// d_sig are returned as zero. The Fortran original calls matmul on an
// uninitialized DDtan in that branch, which is a latent bug; mirroring
// it would expose Path B callers to undefined values whenever the
// driver passes the test=3 sentinel `youngel = -100`. Path B's
// behaviour here is "no elastic update" and is what the umat driver
// would have wanted in that branch.

struct ElasticIncrementResult {
    double D[6][6];     // elastic tangent stiffness, Voigt order
    double d_sig[6];    // stress increment = D * deps
};

ElasticIncrementResult calc_elasti(const double* deps,
                                   double youngel,
                                   double nuel);

// =====================================================================
// RKF residual norm: `norm_res`
// =====================================================================
//
// Mirrors `norm_res_h` in `00_hypo_is_original.for` (.for:1342-1421).
// Used by the RKF23 substepping driver as the Hull-style error
// estimate that decides whether to accept or reject a substep.
//
// Computes a global relative-error norm between two ODE state vectors
// y_til (lower-order RKF estimate) and y_hat (higher-order):
//   err[0..5]  = |sig_hat - sig_til| / ||sig_hat||_stress
//   err[6..11] = |q_hat   - q_til|   / ||q_hat||_strain     (IS components)
//   err[12]    = |void_hat - void_til| / void_hat           (void ratio)
//   norm_R     = sqrt(sum_i err[i]^2)
// Components with norm_sig == 0 or norm_q == 0 contribute zero (the
// Fortran skips the divide). The void slot is divided by void_hat
// unconditionally; void_hat == 0 propagates Inf/NaN into norm_R.
//
// NaN sentinel (matches Fortran .for:1412-1418): if any of norm_sig /
// norm_q / void_hat fails the umatisnan_h finiteness/overflow guard
// (NaN, +/-Inf, or |v| > 1e30), norm_R is forced to 1.0e20 so the RKF
// driver rejects the substep.
//
// Defensive deviation from Fortran: this routine ALSO post-checks the
// returned norm_R itself with the same guard. The Fortran's scan
// inspects only the intermediate norms, missing two corruption paths
// that yield a non-finite norm_R while the intermediates stay clean
// (finite zero void_hat with non-zero del_void; NaN seeded in y_til
// rather than y_hat). Catching them at the source rather than relying
// on every caller to recheck. See Codex review 2026-05-01.

struct NormResResult {
    double norm_R;        // global relative error norm
    bool   nan_detected;  // true iff the 1e20 sentinel fired
};

NormResResult norm_res(const double* y_til, const double* y_hat);

// =====================================================================
// Adaptive RKF23 substepping driver: `rkf23_update`
// =====================================================================
//
// Mirrors `rkf23_update_h` in `00_hypo_is_original.for` (.for:1569-1764).
// Integrates y'(t) = f(y) over the unit normalized interval [0, 1]
// using the explicit Tamagnini-style RKF(2,3) scheme with Hull-error
// step-size control.
//
// Time-step parameterization (Fortran-faithful):
//   - DT_k is the NORMALIZED substep size, T_k in [0, 1].
//   - dtsub is the ABSOLUTE-time suggestion. DT_k = dtsub / dtime.
//   - DTmin is also in normalized space (compared against DT_k).
//
// Substep body (every iteration of the do-while loop runs THIS sequence):
//   1.  ksubst > maxnint guard.
//   2.  check_RKF(y_k)           -> if fail: error = 3 return
//   3.  rhs(y_k)     -> kRK_1    -> if error = 10: propagate
//   4.  y_2 = y_k + (DT_k/2) kRK_1
//   5.  check_RKF(y_2)           -> if fail: error = 3
//   6.  rhs(y_2)     -> kRK_2    -> if error = 10: propagate
//   7.  y_3 = y_k - DT_k kRK_1 + 2 DT_k kRK_2
//   8.  check_RKF(y_3)           -> if fail: error = 3
//   9.  rhs(y_3)     -> kRK_3    -> if error = 10: propagate
//   10. y_til = y_k + DT_k kRK_2
//       y_hat = y_k + DT_k (kRK_1/6 + 2 kRK_2/3 + kRK_3/6)
//   11. norm_R = norm_res(y_til, y_hat)
//       11b. If `norm_res` reports `nan_detected` (the residual
//            computation hit NaN/Inf/overflow), reject the substep
//            with error = 3 + failure_class = ContractViolation
//            BEFORE the tolerance compare. Without this gate, a
//            caller passing `err_tol > 1.0e20` would satisfy
//            `norm_R = sentinel(1.0e20) < err_tol` and commit a step
//            whose error estimate was explicitly poisoned. Same
//            logic as check_RKF's corruption-class precedence:
//            explicit corruption signals win over tolerance compare.
//            Codex review 2026-05-02 (P2 finding); failure_class
//            propagation through Rkf23Result added in 2026-05-02
//            P1 follow-up so integrate_step does not demote this
//            to RkfReject.
//   12. check_RKF(y_hat)         -> if fail: error = 3
//       (preserved verbatim from Fortran .for:1704; this final domain
//       check prevents accepting a substep whose error norm is small
//       but whose state has already left the admissible region.)
//   13. Hull S = 0.9 DT_k (err_tol/norm_R)^(1/3) if norm_R != 0 else 1.
//       Accept (norm_R < err_tol, strict): commit y_hat to y_k; advance
//       T_k by DT_k; expand DT_k via min(4 DT_k, S_hull); WRITE BACK
//       dtsub = DT_k * dtime; THEN clamp DT_k by min(1 - T_k, DT_k).
//       (The dtsub write happens before the final-interval clamp; the
//       caller's "next step suggestion" thus reflects the post-expand
//       value, not the trimmed value.)
//       Reject (norm_R >= err_tol): DT_k = max(DT_k/4, S_hull); if
//       DT_k < DTmin, error = 3.
//
// nfev semantics: incremented BEFORE each rhs call (mirrors Fortran's
// `nfev = nfev + 1` as the first line of rhs_h). On error = 10 from
// rhs, nfev has already been incremented for that call.
//
// Output y semantics on error paths: when error != 0, r.y is left at
// the function's input y_init (the Fortran does not write y_k back to
// y on early-return paths). On error = 0, r.y holds the final accepted
// y_k (fully integrated state).
//
// Input contract (deliberate deviation from Fortran). All scalar tuning
// parameters must be finite and strictly positive; bad values are
// rejected up front with error = 3, nfev = 0, no rhs work performed.
//
//   - dtime: dtime <= 0 with dtsub_init > 0 silently completes in
//     Fortran with a sign-cancelling DT_k = 1 and a negative r.dtsub;
//     dtime == 0 produces NaN DT_k.
//   - dtsub_init: dtsub_init <= 0 in Fortran starts a zero-progress
//     substep loop until maxnint fires.
//   - err_tol: with norm_R == 0 (e.g. zero deps) and invalid err_tol
//     (0/NaN/negative), the strict-< accept check fails AND the
//     S_hull = 1 (norm_R==0) branch makes the reject path's
//     `max(DT_k/4, S_hull)` keep DT_k = 1 forever -- the loop spins
//     until maxnint, wasting up to maxnint*3 rhs evaluations.
//   - DTmin: DTmin <= 0 or NaN disables the reject-shrunk-DT_k early
//     exit, same spin-until-maxnint family of failure.
//
// Plus one silent-clamp rule (also a deviation from Fortran):
//   - dtsub_init > dtime -> silent clamp to dtime. The Fortran has no
//     clamp here; a naive caller that reuses `r.dtsub` (which can
//     carry up to a 4x expansion from the accept branch) would set
//     DT_k > 1 and over-integrate past the unit interval. The clamp
//     mirrors the Fortran accept-branch's own `min(1 - T_k, DT_k)`
//     style, just at the entry point. The OUTPUT `r.dtsub` is NOT
//     clamped (a caller using adaptive wrapper-side dtime can read the
//     unclamped expansion as a "try larger next step" suggestion).
//
// `maxnint` is naturally robust (maxnint <= 0 fails fast on iter 1)
// and is not separately validated. Array inputs (y_init, deps_np1,
// parms) are caller-trusted; parms in particular is expected to come
// from check_parms.
//
// See Codex review 2026-05-01 (high+medium findings).

// Disambiguates the meaning of a non-zero `error` for callers that need
// to route on retryability. Codex review 2026-05-02 (P1 high finding):
// `error` alone is Fortran-faithful (values 0/3/10 mirror umat .for) but
// conflates a retryable RKF reject with a contract violation that looks
// identical (both report `error = 3`). The bridge must escalate
// contract violations to KRATOS_ERROR rather than cutback the outer
// step (a smaller dtime cannot fix corrupt input).
//
// `error` field values stay Fortran-faithful (0/3/10); this enum is a
// parallel signal for routing decisions only. Both `Rkf23Result` and
// `StepResult` carry the field; `integrate_step` passes through the
// value rkf23_update reports.
//
// Path B closure (Codex 2026-05-02 follow-up): `failure_class` is set
// precisely at EVERY error exit -- `integrate_step`'s OWN entry, the
// initial `check_RKF` gate, AND `rkf23_update`'s mid-substep gates
// (including `check_RKF` corruption-class labels in k1/k2/k3/y_hat
// stages and `nres.nan_detected` for poisoned residuals). No
// rkf23-internal corruption is silently demoted to `RkfReject`. This
// is one more deviation from Fortran's int-only error model
// (.for:302-316 conflates all RKF-side error == 3 unconditionally),
// consistent with Path B's broader "explicit fail-closed on every
// detected corruption signal" framework.
enum class FailureClass {
    Ok,                  // error == 0
    RkfReject,           // error == 3, retryable: caller cuts back outer dtime.
                         //   Reserved for adaptive substepping bottoming out
                         //   (DT_k < DTmin) and for substep `check_RKF`
                         //   stress-domain failures (failed = "pmean" or
                         //   "tension"; smaller substep may keep state in
                         //   the admissible region next time).
    ContractViolation,   // error == 3, NOT retryable: input or derived
                         //   state was corrupt (entry pre-flight failure,
                         //   `check_RKF` corruption-class label "nan"/"void"
                         //   at any gate, scalar contract failure inside
                         //   `rkf23_update`, or `nres.nan_detected` poisoned
                         //   residual). Bridge must escalate, not retry.
    Fatal,               // error == 10, get_tan / get_F_sig_q fatal
                         //   (parms-only: `r_uc <= 0`, `fb_temp1 < 0`).
};

struct Rkf23Result {
    double y[kYDim];   // final integrated state (= y_init on error paths)
    double dtsub;      // updated suggested next absolute step
    int    nfev;       // total number of rhs evaluations performed
    // 0 = success (T_k reached 1.0 with all substeps accepted)
    // 3 = RKF reject (maxnint exceeded, check_RKF failed, DT_k < DTmin,
    //     poisoned residual, or scalar contract guard)
    // 10 = fatal from get_tan / get_F_sig_q (propagated through rhs)
    int    error;
    // Mirrors StepResult::failure_class. Set at every error exit so
    // integrate_step can pass it through unchanged. See FailureClass
    // enum docstring above for the value -> meaning mapping.
    FailureClass failure_class;
};

Rkf23Result rkf23_update(const double* y_init,
                         double dtsub_init,
                         double err_tol,
                         int    maxnint,
                         double DTmin,
                         const double* deps_np1,
                         const Parms16& parms,
                         double dtime);

// =====================================================================
// Step wrapper: `integrate_step`
// =====================================================================
//
// Path B's pure-C++ analog of the Fortran `umat` driver: a step-level
// composer that wires together the kernel pieces (check_RKF +
// rkf23_update + calc_elasti + get_tan) into one call that a
// constitutive-law caller can invoke per time increment. It does NOT
// touch Kratos types, statev re-encoding, or pore pressure -- those
// belong in the bridge layer that wraps this routine.
//
// Internal flow (mirrors Fortran umat .for:19-386 minus the testing
// hooks, statev re-encoding, and pore-pressure plumbing -- the
// inittension routing IS reproduced below in step 3a):
//   1. Pack y_init = (sig_n, q_n) -> 13-slot vector.
//   2. check_RKF(y_init, parms) -> detect inittension.
//   3a. inittension fallback. Triggered ONLY when check_RKF fails with
//       failed == "pmean" or failed == "tension" (the stress-domain
//       failure modes that match the Fortran-original inittension
//       semantics). Note: check_RKF's class-precedence rule (see its
//       docstring) guarantees that corruption-class labels ("nan",
//       "void") override stress-domain labels, so a "tension"
//       reaching this branch is provably NOT masking an underlying
//       overflow / NaN / non-positive-void condition. Codex review
//       2026-05-02 (P1: derived overflow leak through tension label).
//        - calc_elasti(deps_np1, youngel = 100, nuel = 0.48)
//        - r.sig = sig_n + d_sig
//        - r.q   = q_n (unchanged)
//        - r.D   = D_elastic
//        - r.error = 0 (the elastic path is intentional, NOT a reject)
//        - r.used_elastic = true
//        - r.nfev = 0
//       Any other check_RKF failure (failed == "nan" for NaN/Inf/
//       overflow in y_init OR derived-quantity overflow surfaced by
//       the corruption-class precedence rule, failed == "void" for
//       non-positive void ratio, or any future mode) is a caller-side
//       CONTRACT VIOLATION and must NOT be silently absorbed into the
//       elastic path (silent corruption would persist across step
//       boundaries). For those, integrate_step fails closed:
//        - r.error = 3 (the Fortran-faithful int code; see also
//          r.failure_class)
//        - r.failure_class = FailureClass::ContractViolation
//          (DISTINCT from RkfReject; bridge MUST escalate as
//          KRATOS_ERROR, MUST NOT cutback the outer dtime -- a
//          smaller dtime cannot fix corrupt input)
//        - r.sig = sig_n, r.q = q_n (input state preserved)
//        - r.D = 0, r.nfev = 0
//        - r.used_elastic = false
//       See Codex review 2026-05-01 (high finding) for the original
//       fail-closed deviation, and 2026-05-02 (P1 high follow-up)
//       for the FailureClass disambiguation that makes the bridge
//       routing unambiguous.
//   3b. Normal RKF path (check_RKF.error == 0):
//        - rkf23_update(...) -> result with error = 0 / 3 / 10
//        - error = 0: commit y_final into r.sig / r.q. The output
//          tangent r.D is then computed at the BEGIN-OF-STEP state
//          (sig_n, q_n) -- NOT at y_final -- to match the Fortran
//          reference perturbate_h, which extracts sig/q from y_n
//          before calling get_tan_h (.for:1466-1475). Path A vs Path B
//          tangent agreement requires this convention exactly.
//          Fortran perturbate_h selects (.for:1478-1490):
//            istrain = 0 -> D = L
//            istrain = 1 -> D = m_R * L
//          (NOT the full IS-modified M.)
//        - error = 3: RKF rejected substep; r.sig = sig_n,
//          r.q = q_n (unchanged); r.D = 0. Caller should cut back
//          its outer time step.
//        - error = 10: fatal from get_tan / get_F_sig_q (today: r_uc
//          guard or fb_temp1 guard); same y revert + D = 0; the bridge
//          layer maps this to KRATOS_ERROR.
//        - r.used_elastic = false in all three sub-cases.
//   4. dtsub write-back clamp on the OUTPUT side (independent of which
//      branch ran):
//        if (dtsub_next <= 0) dtsub_next = 0;
//        else if (dtsub_next >= dtime) dtsub_next = dtime;
//      The kernel `rkf23_update` allows `r.dtsub` up to 4 * dtime
//      as a "try larger" suggestion; a step-level caller persists
//      dtsub across calls, and an unclamped value would trigger
//      rkf23_update's input clamp on the next invocation. Clamping
//      here keeps the persisted value semantically valid.
//
// Voigt order: this routine takes and returns Abaqus-Voigt sig
// (`[11, 22, 33, 12, 13, 23]`) and the kernel-internal q. The bridge
// layer is responsible for the Kratos slot-4/5 swap.
//
// Default tuning parameters match the Fortran umat constants:
//   err_tol = 1.0e-3, maxnint = 10000, DTmin = 1.0e-17 (normalized).
//
// Input contract on scalar tuning parameters:
//   - dtime, err_tol, DTmin must each be finite and strictly positive.
//     Bad values reject up front with error = 3 and r.dtsub_next = 0.
//     The elastic-fallback branch bypasses rkf23_update, so without
//     this entry-side check an invalid scalar would leak past rkf23's
//     own guard and corrupt the wrapper's output (e.g. negative dtime
//     pinning r.dtsub_next negative). See Codex review 2026-05-01
//     (medium follow-up).
//   - maxnint must be > 0 (integer count of allowed substeps).
//     `maxnint <= 0` is a contract violation: rkf23_update would
//     return error = 3 on its first ksubst > maxnint check, and
//     without this guard the wrapper would route it as a retryable
//     RkfReject -- a bridge that cuts back the outer dtime would
//     loop forever (smaller outer step cannot fix an invalid
//     substep budget). See Codex review 2026-05-02 (P2 finding).
//   - dtsub_in must be finite. Non-finite (NaN/Inf) is a hard reject
//     with error = 3 and r.dtsub_next = 0.
//   - dtsub_in is otherwise NORMALIZED, not rejected. Mirroring
//     Fortran umat .for:243-244, any finite `dtsub_in <= 0` or
//     `dtsub_in > dtime` is silently set to `dtime` before integration.
//     This is the legitimate first-call default: Kratos / Abaqus
//     zero-initialize persisted statev, so the very first call to
//     `integrate_step` arrives with `dtsub_in == 0` -- treating that
//     as a contract violation deadlocks the bridge (caller cuts back,
//     reads zero again, never makes progress). See Codex review
//     2026-05-02 (high finding).
//
// All scalar contract failures classify as
// `FailureClass::ContractViolation` -- the bridge MUST escalate as
// KRATOS_ERROR, NOT cutback. Sets r.error = 3 + r.dtsub_next = 0.
//
// Input contract on array inputs: every entry of sig_n (6), q_n (7),
// and deps_np1 (6) must be finite with |v| <= 1e30. q_n[6] (void
// ratio) must additionally be > 0. The deps_np1 vector must also
// satisfy ||deps_np1||_strain <= 1e30 (StrainLike norm), matching
// Fortran umat's `umatisnan_h(norm_D)` sentinel at .for:201-209;
// per-component finiteness is necessary but not sufficient because
// pairs of large finite components (~8e29 each) can still overflow
// the norm. The elastic-fallback branch passes deps_np1 directly to
// `calc_elasti`, which has no error channel, so without these guards
// a NaN/Inf strain increment OR an overflow-norm strain increment
// in a tensile initial state would silently corrupt r.sig with
// error == 0 and used_elastic == true. See Codex review 2026-05-02
// (high finding) and 2026-05-02 follow-up (norm-level overflow).

struct StepResult {
    double sig[6];          // updated EFFECTIVE Cauchy stress
                            // (tension-positive, Voigt-Abaqus). Bridge
                            // converts to total via
                            //   sig_total_next = sig_eff_next - pore_next * I_normal
                            // (normal slots only; see memory bridge
                            // contract pore equations). The kernel
                            // never sees total stress or pore.
    double q[kStateDim];    // updated kernel state: q[0..5] = del
                            // (intergranular strain, Voigt-Abaqus,
                            // continuum sign -- persist directly to
                            // external statev with NO sign flip),
                            // q[6] = void ratio. Equal to q_n on
                            // inittension / error paths.
    double D[6][6];         // ddsdde-style tangent stiffness. Zero on
                            // error = 3 / error = 10 paths.
    double dtsub_next;      // suggested next absolute step. Already
                            // clamped to [0, dtime] -- safe to persist.
    int    nfev;            // total rhs evaluations performed.
    int    error;           // 0 OK, 3 RKF reject OR contract violation, 10 fatal.
                            // For routing, prefer `failure_class` -- `error`
                            // alone cannot disambiguate the two error == 3
                            // sub-cases. Field kept for Fortran-faithful
                            // value-mapping (umat .for uses 0/3/10).
    bool   used_elastic;    // true iff inittension fallback fired.
    FailureClass failure_class;  // Bridge routing signal; see enum docstring.
};

StepResult integrate_step(const double* sig_n,
                          const double* q_n,
                          const double* deps_np1,
                          const Parms16& parms,
                          double dtime,
                          double dtsub_in,
                          double err_tol = 1.0e-3,
                          int    maxnint = 10000,
                          double DTmin   = 1.0e-17);

}  // namespace SandHypoCpp
}  // namespace Kratos

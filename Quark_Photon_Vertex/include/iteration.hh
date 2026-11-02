#pragma once

#include <complex>
#include <chrono>
#include <type_traits>
#include "omp.h"

#include "Utils.hh"
#include "types.hh"
#include "QuadratureIntegral.hh"
#include "ChebyshevPolynomial2.hh"
#include "LinearInterpolate.hh"
#include "complex_spline.hh"
#include "hdf5IO.hh"

#include "parameters.hh"
#include "flavor.hh"
#include "Kernels_G.hh"
#include "Kernels_K.hh"
#include "momentumtransform.hh"
#include "basistransform.hh"
#include "maris_tandy.hh"
#include "WTI.hh"

// b-slices needed by the HVP loop (b1, b7, b10 ↔ struct indices 0, 6, 9),
// returned in memory from iterate_a_and_b so hvp.hh no longer re-reads files.
// Each is indexed [q_iter][k_idx][z_idx].
struct HvpBInput {
  tens_cmplx b1, b7, b10;
  double z2;  // quark wave-function renormalisation (bare outer HVP vertex)
};

// The k/y/z quadratures are all open-coded against hoisted node/weight tables
// (see iterate_a_and_b and precalculate_K_kernel), so no qIntegral object is
// instantiated here any more.
//
// z dimension uses Gauss–Chebyshev type 2: the √(1−z²) Jacobian of the
// 4D loop measure is absorbed into the quadrature weights, so it must
// NOT appear in the integrand and the z-bounds must be (-1, 1).

double update_accuracy(const unsigned z_0, const tens_cmplx &a, const mat_cmplx &a_old)
{
  double current_acc = 0.;

  using namespace parameters::numerical;
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    {
      std::complex<double> current_a = 0.;
      for (unsigned z_idx = 0; z_idx < a[i][k_idx].size(); ++z_idx)
        current_a += std::abs(a[i][k_idx][z_idx]) / double(a[i][k_idx].size());
      const auto diff = current_a - a_old[i][k_idx];
      const auto sum = current_a + a_old[i][k_idx];
      if(!isEqual(std::abs(sum), 0.))
        current_acc = std::max(current_acc, abs(diff) / std::abs(sum));
    }
  return current_acc;
}

mat_cmplx average_array_full(const tens_cmplx &a)
{
  using namespace parameters::numerical;
  mat_cmplx a_z0(a.size(), vec_cmplx(a[0].size(), 0.0));

  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k_idx = 0; k_idx < a[i].size(); ++k_idx)
      for (unsigned z_idx = 0; z_idx < a[i][k_idx].size(); ++z_idx)
        a_z0[i][k_idx] += std::abs(a[i][k_idx][z_idx]) / double(a[i][k_idx].size());
  return a_z0;
}

mat_cmplx average_array_z0(const tens_cmplx &a, const unsigned z_0)
{
  using namespace parameters::numerical;
  mat_cmplx a_z0(a.size(), vec_cmplx(a[0].size(), 0.0));

  for (unsigned i = 0; i < a.size(); ++i)
    for (unsigned k_idx = 0; k_idx < a[i].size(); ++k_idx)
      a_z0[i][k_idx] = 0.5 * (a[i][k_idx][z_0-1] + a[i][k_idx][z_0]);
  return a_z0;
}

double a0 (const unsigned& i)
{
  if(i == 0)
    return std::sqrt(2.);
  else if(i == 6 || i == 9)
    return 1.;
  return 0.;
}

// Initializes a_i with the bare vertex (inhomogeneous term a_i^0).
// See Eq. (48) in the project description.
template<typename Quark>
void a_initialize(tens_cmplx &a, const Quark& quark)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(2)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        a[i][k_idx][z_idx] = quark.z2() * a0(i);
}

// Build cubic splines of b[j][:, zp_idx] over log(k'²) for each (j, zp_idx).
// Done once per BSE iteration after b_iteration_step refreshes b. The 384
// splines (n_structs × z_steps at default sizes) are cheap: each is an
// O(k_steps) tridiagonal solve. Storage ≈ 4·k_steps doubles per spline,
// total ~1.6 MB at the default grid.
std::vector<std::vector<ComplexSpline>> build_b_splines(
    const tens_cmplx& b, const vec_double& k_grid)
{
  using namespace parameters::numerical;
  std::vector<std::vector<ComplexSpline>> spline_b(
      n_structs, std::vector<ComplexSpline>(z_steps));
  for (unsigned j = 0; j < n_structs; ++j) {
    for (unsigned zp_idx = 0; zp_idx < z_steps; ++zp_idx) {
      vec_cmplx b_slice(k_steps);
      for (unsigned k = 0; k < k_steps; ++k)
        b_slice[k] = b[j][k][zp_idx];
      spline_b[j][zp_idx].set_points(k_grid, b_slice);
    }
  }
  return spline_b;
}

template<typename Quark>
void a_iteration_step(const tens_cmplx & /*b*/,
    const ijtens2_double &K_prime, const vec_double &z_grid, const vec_double &k_grid, tens_cmplx &a,
    const Quark& quark,
    const std::vector<std::vector<ComplexSpline>>& spline_b,
    const vec_double &k_quad_log, const vec_double &k_sq_quad_table)
{
  using namespace parameters::numerical;

  // Cheb2 z'-quadrature: nodes coincide with z_grid by construction (both this
  // and Simulation.cpp build z from ChebyshevPolynomial2<z_steps> on (−1, 1)),
  // so no z-interpolation of b/K' is needed — we just sample the stored grid
  // value at zp_idx.
  static const ChebyshevPolynomial2<z_steps> cp_z;
  static const auto z_weights = cp_z.weights();

  // Legendre k-quadrature weights (nodes/log were built once in
  // iterate_a_and_b and threaded through here). K_prime is tabulated
  // directly at these nodes (see precalculate_K_kernel) so we never
  // interpolate K' along k' — the per-(k_idx,z_idx,j) inner sum is just a
  // weighted dot product over k'_quad and zp.
  static const LegendrePolynomial<k_steps> lp_k;
  static const auto& k_unit_weights = lp_k.weights();

  const double k_a = k_grid[0];
  const double k_b = k_grid[k_steps - 1];
  const double dk = 0.5 * (k_b - k_a);

  std::vector<double> k_quad_powr2(k_steps);
  for (unsigned q = 0; q < k_steps; ++q)
    k_quad_powr2[q] = k_sq_quad_table[q] * k_sq_quad_table[q];

  // Tabulate b at the k'-quadrature nodes once per BSE step.
  //
  // spline_b[j][zp](k_quad_log[q]) depends only on (j, zp, q) — never on
  // (i, k_idx, z_idx) — but it used to be evaluated inside the collapse(3) over
  // those, so every value was recomputed n_structs*k_steps*z_steps = 27,648
  // times (~4,992x redundancy: 1.38e8 spline evaluations per step where 27,648
  // distinct ones exist, each costing two binary searches plus two Horner
  // evaluations). The table is n_structs*z_steps*k_steps complex = 432 KiB at
  // the production grid, and turns the inner q-loop into a contiguous read.
  vec_cmplx b_at_quad(std::size_t(n_structs) * z_steps * k_steps);
  #pragma omp parallel for collapse(2)
  for (unsigned j = 0; j < n_structs; ++j)
    for (unsigned zp_idx = 0; zp_idx < z_steps; ++zp_idx)
    {
      const auto& spline_bj = spline_b[j][zp_idx];
      const std::size_t base = (std::size_t(j) * z_steps + zp_idx) * k_steps;
      for (unsigned q = 0; q < k_steps; ++q)
        b_at_quad[base + q] = spline_bj(k_quad_log[q]);
    }

  #pragma omp parallel for collapse(3)
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
      {
        // Inhomogeneous term: a^0_i (Eq. 48), scaled by Z_2.
        a[i][k_idx][z_idx] = quark.z2() * a0(i);

        for (unsigned j = 0; j < n_structs; ++j)
        {
          if (K::isZeroIndex(i, j))
            continue;

          std::complex<double> integral{0., 0.};
          for (unsigned zp_idx = 0; zp_idx < z_steps; ++zp_idx)
          {
            // Pre-tabulated b(j, zp, k'_q) — same values the spline returned,
            // so the accumulation below is unchanged term for term.
            const std::complex<double>* b_row =
                &b_at_quad[(std::size_t(j) * z_steps + zp_idx) * k_steps];

            std::complex<double> partial{0., 0.};
            for (unsigned q = 0; q < k_steps; ++q)
            {
              const double K_at = K_prime(i, k_idx, z_idx, j, q, zp_idx);
              const auto b_at = b_row[q];
              partial += k_unit_weights[q] * K_at * b_at * k_quad_powr2[q];
            }
            partial *= dk;
            // Cheb2 weight already absorbs √(1−z'²); no extra factor here.
            integral += z_weights[zp_idx] * partial;
          }

          a[i][k_idx][z_idx] += integral * 2.0 * M_PI * int_factors;
        }
      }
}

template<typename Quark>
void b_iteration_step(const tens_cmplx &a, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid, tens_cmplx &b, const Quark& quark)
{
  using namespace parameters::numerical;
  #pragma omp parallel for collapse(2)
  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) 
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      const double &z = z_grid[z_idx];

      // Evaluate Gij
      const G<Quark> g_kernel(k_sq, z, q_sq, quark);
      for (unsigned i = 0; i < n_structs; ++i)
      {
        // Initialize the b's to 0
        b[i][k_idx][z_idx] = 0.0;

        // Add stuff to the b's
        for (unsigned j = 0; j < n_structs; ++j)
          b[i][k_idx][z_idx] += g_kernel.get(i, j) * a[j][k_idx][z_idx];
      }
    }
}

template<typename Quark>
void precalculate_K_kernel(const vec_double &y_nodes,
    const vec_double &y_weights, const double y_dx, const double &q_sq,
    const vec_double &z_grid, const vec_double &k_grid,
    const vec_double &k_sq_table, const vec_double &k_table,
    const vec_double &k_sq_quad_table, const vec_double &k_quad_table,
    const vec_double &s_z_table,
    ijtens2_double &K_prime, const Quark& quark, const bool use_PauliVillars)
{
  using namespace parameters::numerical;

  // --- Set up for the compute region (shared by the OpenACC and OpenMP paths).
  // Extract the three scalar quark quantities so the device region never needs
  // the quark object itself (in the DSE case it holds host-side splines).
  const double q_z2  = quark.z2();
  const double q_eta = quark.eta_mt();
  const double q_lam = quark.lambda_mt();

  // Raw pointers so the OpenACC region maps flat arrays, not std::vector/class
  // objects. All are read-only inside the loop. __restrict__ is essential: it
  // tells nvc the output buffer Kp does not alias any read-only table, so the
  // Kp[flat] store carries no cross-iteration dependence. Without it nvc reports
  // "loop carried dependence of Kp-> prevents parallelization" and downgrades
  // the whole collapse(5) to `loop seq` (single-threaded → ~100× slower).
  const double* __restrict__ k_sq_table_p  = k_sq_table.data();
  const double* __restrict__ k_table_p     = k_table.data();
  const double* __restrict__ k_sq_quad_p   = k_sq_quad_table.data();
  const double* __restrict__ k_quad_p      = k_quad_table.data();
  const double* __restrict__ s_z_table_p   = s_z_table.data();
  const double* __restrict__ z_grid_p      = z_grid.data();
  const double* __restrict__ y_nodes_p     = y_nodes.data();
  const double* __restrict__ y_weights_p   = y_weights.data();

  // K' flat storage. Writing by explicit flat index (same layout as
  // SparseTensor6::flat) keeps the whole write device-friendly.
  double* __restrict__ Kp = K_prime.data();
  [[maybe_unused]] const std::size_t Kp_n = K_prime.size();  // used in the ACC copyout clause
  const std::size_t Kd1 = K_prime.d1();
  const std::size_t Kd2 = K_prime.d2();
  const std::size_t Kd4 = K_prime.d4();
  const std::size_t Kd5 = K_prime.d5();

  // Dense list of the non-zero (i,j) kernel pairs, indexed by slot. We iterate
  // over `slot` directly — a genuine loop induction variable — instead of over
  // (i,j) with an indirect slot_map[i*12+j] lookup. That makes the leading
  // factor of the Kp[flat] write index (slot) provably distinct per iteration,
  // so nvc can prove the loop nest is parallel. With the indirection nvc could
  // not disprove aliasing on the write and reported "loop carried dependence of
  // Kp-> prevents parallelization", generating `loop seq` (single GPU thread,
  // ~100× slower than the CPU). slot == SparseTensor6's ij_to_slot_ value, so
  // the flat layout is unchanged.
  const unsigned n_slots = static_cast<unsigned>(K_prime.n_slots());
  unsigned slot_i[144];
  unsigned slot_j[144];
  for (unsigned i = 0; i < n_structs; ++i)
    for (unsigned j = 0; j < n_structs; ++j) {
      const int s = K_prime.slot_of(i, j);
      if (s >= 0) { slot_i[s] = i; slot_j[s] = j; }
    }

  // Parallelise over ALL of (slot, k, z, k', z'); only the innermost
  // y-quadrature stays sequential, as a genuine per-thread reduction into
  // integral_y. This is the mapping nvc handles correctly. Two earlier variants
  // both silently produced zero on the GPU: leaving the z' loop sequential
  // inside a collapse(4) region, and hoisting the y-loop into an `acc routine`
  // — in each case nvc mis-generated a cross-thread reduction on integral_y.
  // k_prime_idx indexes the Legendre k-quadrature nodes used downstream by
  // a_iteration_step (K_prime is tabulated *at* the BSE quadrature, so no
  // k-interpolation is needed later).
  //
  // OpenACC (nvc++ -acc) offloads this to the GPU; every other compiler falls
  // back to the identical OpenMP CPU loop. Data clauses copy the small read-only
  // tables in and the K' buffer out once per q-point.
#ifdef _OPENACC
  #pragma acc parallel loop independent gang vector collapse(5) \
      copyin(k_sq_table_p[0:k_steps], k_table_p[0:k_steps], \
             k_sq_quad_p[0:k_steps], k_quad_p[0:k_steps], \
             s_z_table_p[0:z_steps], z_grid_p[0:z_steps], \
             y_nodes_p[0:y_steps], y_weights_p[0:y_steps], \
             slot_i[0:n_slots], slot_j[0:n_slots]) \
      copyout(Kp[0:Kp_n])
#else
  #pragma omp parallel for collapse(5)
#endif
  for (unsigned slot = 0; slot < n_slots; ++slot)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        for (unsigned k_prime_idx = 0; k_prime_idx < k_steps; ++k_prime_idx)
          for (unsigned z_prime_idx = 0; z_prime_idx < z_steps; ++z_prime_idx)
          {
            const unsigned i = slot_i[slot];
            const unsigned j = slot_j[slot];

            const double k_sq = k_sq_table_p[k_idx];
            const double k_prime_sq = k_sq_quad_p[k_prime_idx];
            const double k_v = k_table_p[k_idx];
            const double k_prime_v = k_quad_p[k_prime_idx];
            const double k_kp_sqrt = k_v * k_prime_v;
            const double z = z_grid_p[z_idx];
            const double s_z = s_z_table_p[z_idx];
            const double z_prime = z_grid_p[z_prime_idx];
            const double s_z_prime = s_z_table_p[z_prime_idx];

            // y-quadrature over the pre-mapped Legendre nodes. Bounds are fixed
            // for the whole run, so y_nodes/y_weights/y_dx are hoisted once in
            // iterate_a_and_b (removing the per-call std::vector alloc qIntegral
            // did). Accumulation order matches qint1d exactly (dx * w[i]·f(y_i)).
            double integral_y = 0.;
            for (unsigned y_idx = 0; y_idx < y_steps; ++y_idx)
            {
              const double y = y_nodes_p[y_idx];

              // Inlined momentumtransform::l2 with the precomputed sqrt(k²·k'²)
              // and sin terms; matches the old formula modulo associativity.
              const double l_sq = k_sq + k_prime_sq
                  - 2. * k_kp_sqrt * (z * z_prime + y * s_z * s_z_prime);

              // Scalar (device-callable) gluon overloads — no quark object.
              const double gl = use_PauliVillars
                  ? pauli_villars_g(l_sq, q_z2, q_eta, q_lam)
                  : maris_tandy_g(l_sq, q_z2, q_eta, q_lam);

              // K-kernel quantities (identical to the K constructor), then the
              // free-function k_component. The K *class* method dispatch
              // (get()->switch->private method) miscompiles to zero on the nvc
              // device path, while this inline+free-function form is correct on
              // both host and device.
              const double u = k_v * s_z;
              const double uprime = k_prime_v * s_z_prime;
              const double inv_l2 = 1.0 / l_sq;
              const double V = (k_v * z - k_prime_v * z_prime) * inv_l2;
              const double w = u * u * inv_l2;
              const double wprime = uprime * uprime * inv_l2;
              const double X = u * uprime * inv_l2;
              const double kij = k_component(i, j, y, l_sq, u, uprime, V, w, wprime, X);
              integral_y += y_dx * y_weights_p[y_idx] * (gl * kij);
            }

            // Flat write, matching SparseTensor6::flat with slot as the leading
            // (i,j) coordinate.
            const std::size_t flat =
                ((((std::size_t(slot) * Kd1 + k_idx) * Kd2 + z_idx) * Kd4
                    + k_prime_idx) * Kd5 + z_prime_idx);
            Kp[flat] = integral_y;
          }
}

void transform_a_to_fg(tens_cmplx &a, const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  using namespace parameters::numerical;
  using namespace basistransform;

  for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
    for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
    {
      vec_cmplx a_copy(n_structs);
      for(unsigned i = 0; i < n_structs; ++i)
        a_copy[i] = a[i][k_idx][z_idx];

      const double Q = std::sqrt(q_sq);
      const double k_sq = std::exp(k_grid[k_idx]);
      const double k = std::sqrt(k_sq);
      const double& z = z_grid[z_idx];
      const double s = std::sqrt(1. - powr<2>(z));

      a[0][k_idx][z_idx] = f1(Q, s, z, k, a_copy);
      a[1][k_idx][z_idx] = f2(Q, s, z, k, a_copy);
      a[2][k_idx][z_idx] = f3(Q, s, z, k, a_copy);
      a[3][k_idx][z_idx] = f4(Q, s, z, k, a_copy);
      a[4][k_idx][z_idx] = f5(Q, s, z, k, a_copy);
      a[5][k_idx][z_idx] = f6(Q, s, z, k, a_copy);
      a[6][k_idx][z_idx] = f7(Q, s, z, k, a_copy);
      a[7][k_idx][z_idx] = f8(Q, s, z, k, a_copy);

      a[8][k_idx][z_idx] = g1(Q, s, z, k, a_copy);
      a[9][k_idx][z_idx] = g2(Q, s, z, k, a_copy);
      a[10][k_idx][z_idx] = g3(Q, s, z, k, a_copy);
      a[11][k_idx][z_idx] = g4(Q, s, z, k, a_copy);
    }
}

template<typename Quark>
HvpBInput iterate_a_and_b(const vec_double &q_grid, const vec_double &z_grid, const vec_double &k_grid, const vec_double &y_grid, const bool use_PauliVillars, const bool debug,
    const FlavorParams &fp, const bool use_quark_DSE, const std::string &h5_path)
{
  using namespace std::chrono;
  const auto start_time = steady_clock::now();

  using namespace parameters::numerical;
  const unsigned z_0 = z_grid.size() / 2;

  // Flat, row-major result buffers accumulated over the (serial) q-loop and
  // written once at the end via qpv_hdf5::write_run. Layout has q as the
  // leading dimension: [q, struct, k, z] (or [q, struct, k] for z0 arrays).
  // The q-loop and WTI loop are serial, so per-q slices never race.
  const std::size_t sz_qskz = std::size_t(q_steps) * n_structs * k_steps * z_steps;
  const std::size_t sz_qsk  = std::size_t(q_steps) * n_structs * k_steps;
  const std::size_t sz_qwkz = std::size_t(q_steps) * 3 * k_steps * z_steps;
  const std::size_t sz_qwk  = std::size_t(q_steps) * 3 * k_steps;
  qpv_hdf5::cvec fg_buf(sz_qskz), b_buf(sz_qskz), fg_z0_buf(sz_qsk);
  qpv_hdf5::cvec w_buf(sz_qwkz), w_z0_buf(sz_qwk);

  // Flat-index helpers into the [q, n, k, z] and [q, n, k] buffers.
  auto idx_kz = [&](unsigned q_it, unsigned n, unsigned nn, unsigned ki, unsigned zi) {
    return ((std::size_t(q_it) * nn + n) * k_steps + ki) * z_steps + zi;
  };
  auto idx_k = [&](unsigned q_it, unsigned n, unsigned nn, unsigned ki) {
    return (std::size_t(q_it) * nn + n) * k_steps + ki;
  };

  // Do some Legendre Magic
  const Quark quark(fp);

  // Pre-map the y-quadrature once; hoisting it lets precalculate_K_kernel
  // open-code the y-sum. y (cosine of the gluon angle) must run over the full
  // [-1, 1]: the Maris-Tandy IR peak sits near y = ±1, and stopping at the
  // outermost Legendre nodes costs an O(1/y_steps²) WTI violation.
  static const LegendrePolynomial<y_steps> lp_y_nodes;
  constexpr double y_a = -1.;
  constexpr double y_b = 1.;
  const double y_dx = 0.5 * (y_b - y_a);
  const vec_double y_weights = lp_y_nodes.weights();
  const vec_double y_nodes = linearMapTo(lp_y_nodes.zeroes(), -1., 1., y_a, y_b);

  // Hoist the k² = exp(k_grid), √k² and √(1-z²) tables out of the inner
  // loops — they only depend on the (fixed) k/z grids, so the K-kernel
  // precalculation goes from ~10⁹ exp/sqrt evaluations per q-point to a few
  // hundred table lookups.
  vec_double k_sq_table(k_steps);
  vec_double k_table(k_steps);
  for (unsigned i = 0; i < k_steps; ++i) {
    k_sq_table[i] = std::exp(k_grid[i]);
    k_table[i] = std::sqrt(k_sq_table[i]);
  }
  vec_double s_z_table(z_steps);
  for (unsigned i = 0; i < z_steps; ++i)
    s_z_table[i] = std::sqrt(1. - z_grid[i] * z_grid[i]);

  // Tables for the inner k'-quadrature nodes: K_prime is sampled DIRECTLY
  // at the Legendre quadrature points used by a_iteration_step (not at
  // k_grid). This eliminates the linear-in-log(k') interpolation of
  // K_prime that was the remaining BSE-side discretization residual after
  // splining b. The nodes/weights are the same as the qIntegral2d's k slot.
  static const LegendrePolynomial<k_steps> lp_k_quad;
  const auto& k_unit_zeros = lp_k_quad.zeroes();
  const double k_a_log = k_grid[0];
  const double k_b_log = k_grid[k_steps - 1];
  const double k_dk = 0.5 * (k_b_log - k_a_log);
  const double k_mid_log = 0.5 * (k_a_log + k_b_log);

  vec_double k_quad_log(k_steps);
  vec_double k_sq_quad_table(k_steps);
  vec_double k_quad_table(k_steps);
  for (unsigned i = 0; i < k_steps; ++i) {
    k_quad_log[i] = k_mid_log + k_dk * k_unit_zeros[i];
    k_sq_quad_table[i] = std::exp(k_quad_log[i]);
    k_quad_table[i] = std::sqrt(k_sq_quad_table[i]);
  }

  // K_prime is sparse over (i,j): only 26 of 144 (i,j) entries are non-zero
  // (8+4 block decoupling + within-block sparsity from K::isZeroIndex), so
  // SparseTensor6 allocates just the 26 slots (~5.5x memory reduction) —
  // still 1.03 GiB at the production grid.
  //
  // Allocated ONCE for the whole run, not per q-point: precalculate_K_kernel
  // overwrites every stored element for each q (all 26 non-zero (i,j) slots
  // over all k, z, k', z'), so the buffer is safe to reuse and its zero-fill
  // is pure waste. Constructing it inside the q-loop cost 24 allocate +
  // zero-fill + free cycles of 1.03 GiB — ~24.7 GiB of pointless memset.
  ijtens2_double K_prime(k_steps, z_steps, k_steps, z_steps,
                         [](unsigned i, unsigned j) { return K::isZeroIndex(i, j); });

  // loop over q
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const auto iter_start_time = steady_clock::now();

    const double &q_sq = q_grid[q_iter];
    std::cout << "\n_____________________________________\n\n"
      << "Calculation for q^2 = " << q_sq << "\n";

    // define new a, b
    tens_cmplx a(n_structs, k_steps, z_steps);
    tens_cmplx b(n_structs, k_steps, z_steps);

    // Refill the K kernel for this q (K_prime is allocated once, above).
    std::cout << " Calculating K'_ij..." << std::flush;
    const auto kprecalc_t0 = steady_clock::now();
    precalculate_K_kernel(y_nodes, y_weights, y_dx, q_sq, z_grid, k_grid,
                          k_sq_table, k_table,
                          k_sq_quad_table, k_quad_table,
                          s_z_table,
                          K_prime, quark, use_PauliVillars);
    const double kprecalc_ms =
        duration_cast<microseconds>(steady_clock::now() - kprecalc_t0).count() / 1000.;
    std::cout << " done (" << kprecalc_ms << " ms)\n";

    // Initialize a with bare vertex
    a_initialize(a, quark);

    // Main self-consistent iteration: a → b → a until convergence.
    // Implements the coupled system in Eq. (43) of the project description.
    std::cout << "  Starting iteration...\n";
    double current_acc = 1.0;
    unsigned current_step = 0;
    while (max_steps > current_step++ && current_acc > target_acc)
    {
      debug_out("\n    Started a step...\n", debug);
      
      // copy for checking the convergence
      const auto a_old = average_array_full(a);

      debug_out("    Calculating b_i...", debug);
      b_iteration_step(a, q_sq, z_grid, k_grid, b, quark);
      debug_out(" done\n", debug);

      // Refresh the cubic splines of b[j](log k'²) at each z'-grid point.
      // Used inside a_iteration_step to interpolate b at the Legendre
      // quadrature nodes in k' — replaces the linear log-k' interpolation
      // that limited WTI2 accuracy at the IR k corner.
      const auto spline_b = build_b_splines(b, k_grid);

      debug_out("    Calculating a_i...", debug);
      a_iteration_step(b, K_prime, z_grid, k_grid, a, quark, spline_b,
                       k_quad_log, k_sq_quad_table);
      debug_out(" done\n", debug);

      // check the convergence
      current_acc = update_accuracy(z_0, a, a_old);
      debug_out("    current_step = " + std::to_string(current_step) + "\n    current_acc = " + std::to_string(current_acc) + "\n", debug);
    }
    if (current_acc < target_acc)
      std::cout << "  Converged after " << current_step << " steps.\n";
    else
      std::cout << "  ! Did not converge !\n";

    std::cout << " Saving results..." << std::flush;
    // resync b with the converged a, since the last loop step recomputes a
    // from the previous-iteration b. With this call b = G_kernel · a holds for
    // the converged a, which is what the HVP loop expects (b ≡ S·Γ·S).
    b_iteration_step(a, q_sq, z_grid, k_grid, b, quark);
    // stash b into the [q,12,k,z] buffer
    for (unsigned i = 0; i < n_structs; ++i)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
        for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
          b_buf[idx_kz(q_iter, i, n_structs, k_idx, z_idx)] = b[i][k_idx][z_idx];
    // transform to g,f (almost in place!)
    transform_a_to_fg(a, q_sq, k_grid, z_grid);
    const auto& fg = a;
    // stash fg into the [q,12,k,z] buffer
    for (unsigned i = 0; i < n_structs; ++i)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
        for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
          fg_buf[idx_kz(q_iter, i, n_structs, k_idx, z_idx)] = fg[i][k_idx][z_idx];
    // z-average and stash into the [q,12,k] buffer
    const auto fg_z0 = average_array_z0(fg, z_0);
    for (unsigned i = 0; i < n_structs; ++i)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
        fg_z0_buf[idx_k(q_iter, i, n_structs, k_idx)] = fg_z0[i][k_idx];
    std::cout << "  done\n";

    const auto iter_end_time = steady_clock::now();
    std::cout << "Calculation finished after " << duration_cast<milliseconds>(iter_end_time - iter_start_time).count()/1000.<< "s\n";
  }

  std::cout << "\n_____________________________________\n\n"
    << "\nCalculating the WTIs...\n";

  // check WTI
  for (unsigned q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const double &q_sq = q_grid[q_iter];
    tens_cmplx w(3, k_steps, z_steps);
    for (unsigned k_idx = 0; k_idx < parameters::numerical::k_steps; ++k_idx)
    {
      for (unsigned z_idx = 0; z_idx < parameters::numerical::z_steps; ++z_idx)
      {
        const double q_sq = q_grid[q_iter];
        const double Q = std::sqrt(q_sq);
        const double k_sq = std::exp(k_grid[k_idx]);
        const double k = std::sqrt(k_sq);
        const double& z = z_grid[z_idx];

        double kplus2 = k_sq+ q_sq/4. + k*Q*z;
        double kminus2 = k_sq + q_sq/4. - k*Q*z;

        w[0][k_idx][z_idx] = Sigma_A<Quark>(kminus2,kplus2,quark);
        w[1][k_idx][z_idx] = Delta_A<Quark>(kminus2,kplus2,quark);
        w[2][k_idx][z_idx] = Delta_B<Quark>(kminus2,kplus2,quark);
      }
    }
    // stash w into the [q,3,k,z] buffer
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
        for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
          w_buf[idx_kz(q_iter, i, 3, k_idx, z_idx)] = w[i][k_idx][z_idx];
    const auto w_z0 = average_array_z0(w, z_0);
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
        w_z0_buf[idx_k(q_iter, i, 3, k_idx)] = w_z0[i][k_idx];
  }

  // Write everything to a single HDF5 file.
  std::cout << "\nWriting results to " << h5_path << " ..." << std::flush;
  qpv_hdf5::RunMeta meta{
      fp.name, fp.m_c, fp.mu, fp.quark_a0, fp.eta_mt, fp.lambda_mt,
      parameters::physical::lambda_pv, quark.z2(), min_q_sq, max_q_sq,
      n_structs, q_steps, k_steps, z_steps, y_steps,
      use_quark_DSE, use_PauliVillars};
  qpv_hdf5::write_run(h5_path, meta, q_grid, k_sq_table, k_grid, z_grid, y_grid,
                      fg_buf, fg_z0_buf, b_buf, w_buf, w_z0_buf);

  // When the quark propagator was solved from its DSE, also store the raw
  // dressing functions A(p²), B(p²), M(p²) on the DSE grid (group /quark_dse).
  // The analytic quark_model has no such solution, so this is DSE-only.
  if constexpr (std::is_same_v<Quark, quark_DSE>)
    qpv_hdf5::append_quark_dse(h5_path, quark.dse_log_p_sq(),
                               quark.dse_A(), quark.dse_B());
  std::cout << " done\n";

  // Build the b1/b7/b10 slices (struct indices 0, 6, 9) that the HVP loop
  // needs, indexed [q_iter][k_idx][z_idx], straight from the b buffer.
  HvpBInput hvp_in{tens_cmplx(q_steps, k_steps, z_steps),
                   tens_cmplx(q_steps, k_steps, z_steps),
                   tens_cmplx(q_steps, k_steps, z_steps),
                   quark.z2()};
  for (unsigned q_it = 0; q_it < q_steps; ++q_it)
    for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx) {
        hvp_in.b1[q_it][k_idx][z_idx]  = b_buf[idx_kz(q_it, 0, n_structs, k_idx, z_idx)];
        hvp_in.b7[q_it][k_idx][z_idx]  = b_buf[idx_kz(q_it, 6, n_structs, k_idx, z_idx)];
        hvp_in.b10[q_it][k_idx][z_idx] = b_buf[idx_kz(q_it, 9, n_structs, k_idx, z_idx)];
      }

  auto end_time = steady_clock::now();
  std::cout << "\nProgram finished after " << duration_cast<milliseconds>(end_time - start_time).count()/1000.<< "s\n";
  return hvp_in;
}

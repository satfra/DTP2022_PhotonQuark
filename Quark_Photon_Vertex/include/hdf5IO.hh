#pragma once

/*
 * HDF5 output for the quark-photon-vertex run. Replaces the per-index .dat
 * files written by fileIO.hh with a single self-describing .h5 file per run.
 *
 * Layout (one file per run, named output_<flavour>[_dse][_pv].h5):
 *   /grids/q_grid   [q]         external Q^2 grid (linear)
 *   /grids/k_sq     [k]         physical k^2 = exp(k_grid)
 *   /grids/k_grid   [k]         raw log(k^2) grid (for reproducibility)
 *   /grids/z_grid   [z]         z = cos angle grid
 *   /grids/y_grid   [y]         y angular grid
 *   /fg     [q,12,k,z]  complex   f1..f8, g1..g4 dressing functions
 *   /fg_z0  [q,12,k]    complex   z-averaged f/g
 *   /b      [q,12,k,z]  complex   b_i = (S.Gamma.S) products
 *   /w      [q,3,k,z]   complex   WTI quantities Sigma_A, Delta_A, Delta_B
 *   /w_z0   [q,3,k]     complex   z-averaged WTI quantities
 *   /hvp/p_sq [q], /hvp/Pi [q]    (written later by append_hvp, only with -h)
 *
 * Root attributes carry all run metadata (flavour, masses, MT params, grid
 * sizes, flags, and the idx->name maps). Complex data uses a compound
 * {r,i} datatype so h5py reads it as numpy complex128 (see hdf5_complex.hh).
 *
 * The whole run is accumulated in memory and written once (~21 MB at the
 * default grid) — see iterate_a_and_b in iteration.hh.
 */

#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include "hdf5_complex.hh"   // <hdf5lib/hdf5.hh> + std::complex TypeTrait
#include "types.hh"

namespace qpv_hdf5
{
  using cvec = std::vector<std::complex<double>>;

  // All scalar metadata for one run, written as root attributes.
  struct RunMeta {
    std::string flavor;
    double m_c, mu, quark_a0, eta_mt, lambda_mt, lambda_pv;
    double z2;                 // renormalisation constant Z_2 actually used
    double min_q_sq, max_q_sq;
    unsigned n_structs, q_steps, k_steps, z_steps, y_steps;
    bool use_dse, use_pauli_villars;
  };

  namespace detail
  {
    // Create a dataset shaped `dims` and write a flat complex buffer into it.
    inline void write_complex(hdf5::Group& g, const std::string& name,
                              const hdf5::Dims& dims, const cvec& buf)
    {
      auto ds = g.create_dataset(name, hdf5::type_of<std::complex<double>>(),
                                 hdf5::Dataspace::simple(dims));
      ds.write(buf.data(), buf.size());
    }

    inline void write_double(hdf5::Group& g, const std::string& name,
                             const vec_double& buf)
    {
      auto ds = g.create_dataset(name, hdf5::type_of<double>(),
                                 hdf5::Dataspace::simple({buf.size()}));
      ds.write(buf);
    }
  } // namespace detail

  // Write the full run to `path`, truncating any existing file. The complex
  // buffers are flat and row-major with q as the leading dimension.
  inline void write_run(const std::string& path, const RunMeta& m,
                        const vec_double& q_grid, const vec_double& k_sq,
                        const vec_double& k_grid, const vec_double& z_grid,
                        const vec_double& y_grid,
                        const cvec& fg, const cvec& fg_z0, const cvec& b,
                        const cvec& w, const cvec& w_z0)
  {
    using hdf5::Dims;
    auto file = hdf5::File::open(path, hdf5::Access::Truncate);
    auto root = file.root();

    // grids
    auto grids = root.create_group("grids");
    detail::write_double(grids, "q_grid", q_grid);
    detail::write_double(grids, "k_sq",   k_sq);
    detail::write_double(grids, "k_grid", k_grid);
    detail::write_double(grids, "z_grid", z_grid);
    detail::write_double(grids, "y_grid", y_grid);

    // result tensors (q as leading dim)
    const hsize_t q = m.q_steps, ns = m.n_structs, nw = 3;
    const hsize_t k = m.k_steps, z = m.z_steps;
    detail::write_complex(root, "fg",    Dims{q, ns, k, z}, fg);
    detail::write_complex(root, "fg_z0", Dims{q, ns, k},    fg_z0);
    detail::write_complex(root, "b",     Dims{q, ns, k, z}, b);
    detail::write_complex(root, "w",     Dims{q, nw, k, z}, w);
    detail::write_complex(root, "w_z0",  Dims{q, nw, k},    w_z0);

    // metadata attributes
    root.write_attribute("flavor", m.flavor);
    root.write_attribute("m_c", m.m_c);
    root.write_attribute("mu", m.mu);
    root.write_attribute("quark_a0", m.quark_a0);
    root.write_attribute("eta_mt", m.eta_mt);
    root.write_attribute("lambda_mt", m.lambda_mt);
    root.write_attribute("lambda_pv", m.lambda_pv);
    root.write_attribute("z2", m.z2);
    root.write_attribute("min_q_sq", m.min_q_sq);
    root.write_attribute("max_q_sq", m.max_q_sq);
    root.write_attribute("n_structs", m.n_structs);
    root.write_attribute("q_steps", m.q_steps);
    root.write_attribute("k_steps", m.k_steps);
    root.write_attribute("z_steps", m.z_steps);
    root.write_attribute("y_steps", m.y_steps);
    root.write_attribute("use_dse", int(m.use_dse));
    root.write_attribute("use_pauli_villars", int(m.use_pauli_villars));

    const std::vector<std::string> fg_names =
        {"f1","f2","f3","f4","f5","f6","f7","f8","g1","g2","g3","g4"};
    const std::vector<std::string> w_names = {"Sigma_A","Delta_A","Delta_B"};
    {
      auto ds = root.create_dataset("fg_index_names", hdf5::type_of<std::string>(),
                                    hdf5::Dataspace::simple({fg_names.size()}));
      ds.write(fg_names);
    }
    {
      auto ds = root.create_dataset("w_index_names", hdf5::type_of<std::string>(),
                                    hdf5::Dataspace::simple({w_names.size()}));
      ds.write(w_names);
    }
  }

  // Append the HVP result to an existing run file (group /hvp).
  inline void append_hvp(const std::string& path, const vec_double& p_sq,
                         const vec_double& pi)
  {
    auto file = hdf5::File::open(path, hdf5::Access::ReadWrite);
    auto root = file.root();
    auto g = root.create_group("hvp");
    detail::write_double(g, "p_sq", p_sq);
    detail::write_double(g, "Pi", pi);
  }

  // Append the quark-DSE solution to an existing run file (group /quark_dse):
  // the propagator dressing functions A(p^2), B(p^2) and M = B/A on the DSE's
  // own momentum grid. `log_p_sq` is log(p^2); we also store p_sq = exp(...).
  inline void append_quark_dse(const std::string& path, const vec_double& log_p_sq,
                               const vec_double& A, const vec_double& B)
  {
    vec_double p_sq(log_p_sq.size()), M(A.size());
    for (std::size_t i = 0; i < log_p_sq.size(); ++i) p_sq[i] = std::exp(log_p_sq[i]);
    for (std::size_t i = 0; i < A.size(); ++i) M[i] = (A[i] != 0.0) ? B[i] / A[i] : 0.0;

    auto file = hdf5::File::open(path, hdf5::Access::ReadWrite);
    auto root = file.root();
    auto g = root.create_group("quark_dse");
    detail::write_double(g, "p_sq", p_sq);
    detail::write_double(g, "log_p_sq", log_p_sq);
    detail::write_double(g, "A", A);
    detail::write_double(g, "B", B);
    detail::write_double(g, "M", M);
  }
} // namespace qpv_hdf5

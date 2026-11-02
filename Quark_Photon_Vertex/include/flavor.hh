#pragma once

#include <iostream>
#include <string>

/*
 * Per-flavour physical parameters for the quark-photon-vertex calculation.
 *
 * These used to be compile-time constants in parameters.hh (m_c, mu, quark_a0,
 * eta_mt, lambda_mt). They are now a runtime value object so a single binary
 * can compute any flavour, selected by a CLI flag. The constexpr globals in
 * parameters.hh are left in place (dead) so parameters_bench.hh is untouched.
 *
 * Parameter provenance — all masses at the renormalisation scale mu = 19 GeV,
 * in the Fischer/Giessen (Lambda, eta) convention (omega = Lambda/eta,
 * D = eta*Lambda^2). The whole table is sourced consistently from one group:
 *
 *   flavour   m_c      eta     Lambda  fit target   source
 *   light     0.0037   1.8     0.72    pi_0         arXiv:1406.4370 Sec. 5.1
 *   strange   0.085    1.8     0.72    kaon         arXiv:1406.4370 Sec. 5.2
 *   charm     0.870    1.157   0.72    J/psi,...    arXiv:1409.5076 Sec. 3.1.1
 *   bottom    3.790    1.357   0.72    Upsilon(1S)  arXiv:1409.5076 Sec. 3.2
 *
 * Light/strange keep the light-sector eta = 1.8; charm/bottom use the
 * Fischer-retuned eta. Lambda = 0.72 GeV is fixed for every flavour. Edit the
 * values in flavor_from_flag() below to retune.
 */
struct FlavorParams {
  std::string name;    // "light", "strange", "charm", "bottom"
  double m_c;          // current-quark mass at mu (GeV)
  double mu;           // renormalisation scale (GeV)
  double quark_a0;     // initial A dressing value
  double eta_mt;       // Maris-Tandy IR shape parameter (eta)
  double lambda_mt;    // Maris-Tandy IR scale Lambda (GeV)
};

// Select a flavour from a single character (l/s/c/b). Unknown -> light + warn.
inline FlavorParams flavor_from_flag(char c)
{
  switch (c) {
    case 'l': return {"light",   0.0037, 19.0, 1.0, 1.8,   0.72};
    case 's': return {"strange", 0.085,  19.0, 1.0, 1.8,   0.72};
    case 'c': return {"charm",   0.870,  19.0, 1.0, 1.157, 0.72};
    case 'b': return {"bottom",  3.790,  19.0, 1.0, 1.357, 0.72};
    default:
      std::cerr << "flavor_from_flag: unknown flavour '" << c
                << "', defaulting to light.\n";
      return {"light", 0.0037, 19.0, 1.0, 1.8, 0.72};
  }
}

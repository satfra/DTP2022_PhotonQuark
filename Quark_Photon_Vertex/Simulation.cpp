#include <iostream>
#include <numeric>

#include "Utils.hh"
#include "LegendrePolynomials.hh"
#include "quark_model_functions.hh"
#include "iteration.hh"
#include "parameters.hh"

using namespace std;

#define LOG(x) cout << x << endl

int main(int argc, char*argv[])
{
    std::cout << "\n\n#############     CALCULATION OF QUARK-PHOTON-VERTEX     #############\n\n\n";

    // ++++++++++++++++++++++++++++ get flags from shell ++++++++++++++++++++++++++++

    std::string flags = argc > 1 ? argv[1] : "";

    const bool debug = flags.find('v') < flags.length() ? true : false;
    if(debug) std::cout << "Showing debug output.\n";

    const bool use_quark_DSE = flags.find('d') < flags.length() ? true : false;
    if(use_quark_DSE) std::cout << "Using the quark DSE.\n";

    const bool use_PauliVillars = flags.find('p') < flags.length() ? true : false;
    if(use_PauliVillars) std::cout << "Using Pauli-Villars regularisation.\n";


    // ++++++++++++++++++++++++++++ create grids ++++++++++++++++++++++++++++

    // avoid z == 0 in a grid, which would lead to division by zero.
    static_assert(parameters::numerical::z_steps % 2 == 0);
    
    // create k^2 grid - use legendre nodes on logarithmic scale
    LegendrePolynomial<parameters::numerical::k_steps> lp_k;
    const vector<double> gauleg_nodes_k = lp_k.zeroes();
    std::vector<double> log_k_sq_grid_gauleg = linearMapTo(gauleg_nodes_k, gauleg_nodes_k.front(), gauleg_nodes_k.back(), 
        std::log(parameters::physical::lambda_IR), std::log(parameters::physical::lambda_UV));

    // create logarithmically spaced q grid
    const vector<double> q_sq_grid = logGrid(parameters::numerical::q_steps,parameters::numerical::min_q_sq, parameters::numerical::max_q_sq);
    
    // We use the zeroes of LegendrePolynomials for the z and y grids
    LegendrePolynomial<parameters::numerical::z_steps> lp_z;
    LegendrePolynomial<parameters::numerical::y_steps> lp_y;
    const vector<double> z_grid = lp_z.zeroes();
    const vector<double> y_grid = lp_y.zeroes();


    // ++++++++++++++++++++++++++++ solve quark photon vertex ++++++++++++++++++++++++++++

    if(use_quark_DSE)
        solve_BSE<quark_DSE>(q_sq_grid, z_grid, log_k_sq_grid_gauleg, y_grid, use_PauliVillars, debug);
    else
        solve_BSE<quark_model>(q_sq_grid, z_grid, log_k_sq_grid_gauleg, y_grid, use_PauliVillars, debug);


    return 0;

}

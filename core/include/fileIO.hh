#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include "Utils.hh"
#include "parameters.hh"

void saveToFile(const std::vector<double>& data, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(file);
  for(unsigned i = 0; i < data.size(); ++i)
    fileStream << data[i] << "\n";
  fileStream.close();
}

void saveToFile(const qtens_cmplx& data, const std::string& file, const std::string& header)
{
  // TODO: Make me print the values for k, z and q instead of only their indices!
  std::ofstream fileStream;
  fileStream.open(file);
  fileStream << header << "\n";
  for(unsigned i = 0; i < data.size(); ++i)
    for(unsigned j = 0; j < data[i].size(); ++j)
      for(unsigned k = 0; k < data[i][j].size(); ++k)
        for(unsigned l = 0; l < data[i][j][k].size(); ++l)
        fileStream << i << " " << j << " " << k << " " << l << " " << data[i][j][k][l] << "\n";
  fileStream.close();
}

void saveToFile_withGrids(const qtens_cmplx& data, const std::string& file, const std::string& header,
    const vec_double &q_grid, const vec_double &k_grid, const vec_double &z_grid)
{
  std::ofstream fileStream;
  fileStream.open(file);

  fileStream << header << "\n";

  using namespace parameters::numerical;
  for (unsigned int q_iter = 0; q_iter < q_steps; q_iter++)
  {
    const double &q_sq = q_grid[q_iter];
    for(unsigned i = 0; i < n_structs; ++i)
    {
      for (unsigned k_idx = 0; k_idx < k_steps; ++k_idx)
      {
        const double k_sq = std::exp(k_grid[k_idx]);
        for (unsigned z_idx = 0; z_idx < z_steps; ++z_idx)
        {
          const double& z = z_grid[z_idx];
          fileStream << q_sq << " " << i << " " << k_sq << " " << z << " " << data[q_iter][i][k_idx][z_idx].real() << " " << data[q_iter][i][k_idx][z_idx].imag() << "\n";
        }
      }
    }
  }
  fileStream.close();
}

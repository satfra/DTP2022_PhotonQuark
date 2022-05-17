#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <complex>

#include "Utils.hh"
#include "types.hh"

void saveToFile(const std::vector<double>& data, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(file);
  for(unsigned i = 0; i < data.size(); ++i)
    fileStream << i << " " << data[i] << "\n";
  fileStream.close();
}

void emptyFile(const std::string& file, const std::string& header)
{
  std::ofstream fileStream;
  fileStream.open(file, std::ofstream::out | std::ofstream::trunc);
  fileStream << header << "\n";
  fileStream.close();
}

template<unsigned N>
void saveToFile_withGrids(const tens_cmplx& data, const std::string& file,
    const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  std::ofstream fileStream;
  fileStream.open(file, std::ios_base::app);

  for(unsigned i = 0; i < N; ++i)
  {
    for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      for (unsigned z_idx = 0; z_idx < z_grid.size(); ++z_idx)
      {
        const double& z = z_grid[z_idx];
        fileStream << q_sq << " " << i << " " << k_sq << " " << z << " " << data[i][k_idx][z_idx].real() << " " << data[i][k_idx][z_idx].imag() << "\n";
      }
    }
    fileStream << "\n";
  }
  fileStream.close();
}

template<unsigned N>
void saveToFile_withGrids(const mat_cmplx& data, const std::string& file,
    const double& q_sq, const vec_double &k_grid)
{
  std::ofstream fileStream;
  fileStream.open(file, std::ios_base::app);

  for(unsigned i = 0; i < N; ++i)
  {
    for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      fileStream << q_sq << " " << i << " " << k_sq << " " << data[i][k_idx].real() << " " << data[i][k_idx].imag() << "\n";
    }
  }
  fileStream.close();
}

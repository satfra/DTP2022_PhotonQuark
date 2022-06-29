#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <complex>

#include "Utils.hh"
#include "types.hh"

// constant relative path where the files are saved
const std::string pathRel = "./data/";

void saveToFile(const std::vector<double>& data, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(pathRel+file);
  for(unsigned i = 0; i < data.size(); ++i)
    fileStream << i << " " << data[i] << "\n";
  fileStream.close();
}

void emptyFile(const std::string& file, const std::string& header)
{
  std::ofstream fileStream;
  fileStream.open(pathRel + file + ".dat", std::ofstream::out | std::ofstream::trunc);
  fileStream << header << "\n";
  fileStream.close();
}

template<unsigned N>
void emptyIdxFile(const std::string& file, const std::string& header)
{
  for(unsigned i = 0; i < N; ++i)
  {
    std::ofstream fileStream;
    fileStream.open(pathRel+file + "_idx_" + std::to_string(i) + ".dat", std::ofstream::out | std::ofstream::trunc);
    fileStream << header << "\n";
    fileStream.close();
  }
}

void saveToFile_M_Z(const std::vector<double>& M, const std::vector<double>& Z, const vec_double &k_grid, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(pathRel + file + ".dat", std::ios_base::app);
  for(unsigned i = 0; i < M.size(); ++i)
  {
    const double k_sq = std::exp(k_grid[i]);
    fileStream << k_sq << " " << M[i] << " " << Z[i] << "\n";
  }
  fileStream.close();
}

template<unsigned N>
void saveToFile_withGrids(const tens_cmplx& data, const std::string& file,
    const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  for(unsigned i = 0; i < N; ++i)
  {
    std::ofstream fileStream;
    fileStream.open(pathRel + file + "_idx_" + std::to_string(i) + ".dat", std::ios_base::app);
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
    fileStream.close();
  }
}

template<unsigned N>
void saveToFile_withGrids(const dressing* data, const std::string& file,
    const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  for(unsigned i = 0; i < N; ++i)
  {
    std::ofstream fileStream;
    fileStream.open(pathRel + file + "_idx_" + std::to_string(i) + ".dat", std::ios_base::app);
    for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      for (unsigned z_idx = 0; z_idx < z_grid.size(); ++z_idx)
      {
        const double& z = z_grid[z_idx];
        fileStream << q_sq << " " << i << " " << k_sq << " " << z << " " << (*data)[i][k_idx][z_idx].real() << " " << (*data)[i][k_idx][z_idx].imag() << "\n";
      }
    }
    fileStream << "\n";
    fileStream.close();
  }
}

void saveToFile_withGrids_double(const mat_double& data, const std::string& file,
    const double& q_sq, const vec_double &k_grid, const vec_double &z_grid)
{
  
    std::ofstream fileStream;
    fileStream.open(pathRel + file + ".dat", std::ios_base::app);
    for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      for (unsigned z_idx = 0; z_idx < z_grid.size(); ++z_idx)
      {
        const double& z = z_grid[z_idx];
        fileStream << q_sq << " " << k_sq << " " << z << " " << data[k_idx][z_idx] << "\n";
      }
    }
    fileStream << "\n";
    fileStream.close();
  
}

template<unsigned N>
void saveToFile_withGrids(const mat_cmplx& data, const std::string& file,
    const double& q_sq, const vec_double &k_grid)
{
  for(unsigned i = 0; i < N; ++i)
  {
    std::ofstream fileStream;
    fileStream.open(pathRel + file + "_idx_" + std::to_string(i) + ".dat", std::ios_base::app);
    for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
    {
      const double k_sq = std::exp(k_grid[k_idx]);
      fileStream << q_sq << " " << i << " " << k_sq << " " << data[i][k_idx].real() << " " << data[i][k_idx].imag() << "\n";
    }
    fileStream << "\n";
    fileStream.close();
  }
}

template<unsigned N>
void saveToFile_withGrids(const jtens2_cmplx& data, const std::string& file,
    const double& q_sq, const vec_double &k_prime_grid, const vec_double &z_prime_grid, const vec_double &k_grid, const vec_double &z_grid)
{
  for(unsigned i = 0; i < N; ++i)
  {
    std::ofstream fileStream;
    fileStream.open(pathRel + file + "_idx_" + std::to_string(i) + ".dat", std::ios_base::app);
    for (unsigned k_idx1 = 0; k_idx1 < k_prime_grid.size(); ++k_idx1)
    {
      const double k_prime_sq = std::exp(k_prime_grid[k_idx1]);
      for (unsigned z_idx1 = 0; z_idx1 < z_prime_grid.size(); ++z_idx1)
      {  
        const double& z_prime = z_prime_grid[z_idx1];
        for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
        {
          const double k_sq = std::exp(k_grid[k_idx]);
          for (unsigned z_idx = 0; z_idx < z_grid.size(); ++z_idx)
          {
            const double& z = z_grid[z_idx];
            fileStream << q_sq << " " << i << " " << k_prime_sq << " " << z_prime << " " << k_sq << " " << z << " " << data[k_idx1][z_idx1][8+i][k_idx][z_idx].real() << " " << data[k_idx1][z_idx1][8+i][k_idx][z_idx].imag() << "\n";
          }
        }
      }
    }
  fileStream << "\n";
  fileStream.close();
  }
}

void saveToFile_withGrids_mod(const tens_cmplx& data, const std::string& file,
    const vec_double& q_sq_grid, const vec_double &k_grid, const vec_double &z_grid)
{
    std::ofstream fileStream;
    fileStream.open(pathRel + file + ".dat", std::ios_base::app);
    
      for (unsigned q_idx = 0; q_idx < q_sq_grid.size(); ++q_idx)
      {  
        const double& q_sq = q_sq_grid[q_idx];
        for (unsigned k_idx = 0; k_idx < k_grid.size(); ++k_idx)
        {
          const double k_sq = std::exp(k_grid[k_idx]);
          for (unsigned z_idx = 0; z_idx < z_grid.size(); ++z_idx)
          {
            const double& z = z_grid[z_idx];
            fileStream << q_sq << " " << k_sq << " " << z << " " << data[q_idx][k_idx][z_idx].real() << " " << data[q_idx][k_idx][z_idx].imag() << "\n";
          }
        }
      }
  
  fileStream << "\n";
  fileStream.close();
  
}
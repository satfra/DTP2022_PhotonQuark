#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include "Utils.hh"

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

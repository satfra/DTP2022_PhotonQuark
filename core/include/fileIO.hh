#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include "iteration.hh"

void saveToFile(const std::vector<double>& data, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(file, std::ios_base::app);
  for(unsigned i = 0; i < data.size(); ++i)
    fileStream << data[i] << "\n";
  fileStream.close();
}

void saveToFile(const tens_cmplx& data, const std::string& file)
{
  std::ofstream fileStream;
  fileStream.open(file, std::ios_base::app);
  for(unsigned i = 0; i < data.size(); ++i)
    for(unsigned j = 0; j < data[i].size(); ++j)
      for(unsigned k = 0; k < data[j].size(); ++k)
        fileStream << i << " " << j << " " << k << " " << data[i] << "\n";
  fileStream.close();
}

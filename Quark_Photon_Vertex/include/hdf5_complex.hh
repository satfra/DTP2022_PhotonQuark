#pragma once

/*
 * std::complex<double> support for the vendored hdf5lib wrapper.
 *
 * hdf5lib ships TypeTraits for the builtin scalars and std::string but not for
 * complex numbers (that specialisation lived outside the copied directory).
 * We add it here as a compound datatype {r: double, i: double}.
 *
 * std::complex<double> is guaranteed by the standard to have the same layout
 * and alignment as double[2] with [0] = real, [1] = imag, so a raw write of a
 * flat std::complex buffer is valid against this compound type. h5py detects a
 * compound of two floats named exactly "r"/"i" as numpy complex128, so the
 * Python readers get native complex arrays with no post-processing.
 */

#include <hdf5lib/hdf5.hh>

#include <complex>

namespace hdf5
{
  template <> struct TypeTrait<std::complex<double>> {
    static Datatype get()
    {
      Datatype t = Datatype::compound(sizeof(std::complex<double>)); // 16 bytes
      t.insert("r", 0,             type_of<double>());
      t.insert("i", sizeof(double), type_of<double>());
      return t;
    }
  };
} // namespace hdf5

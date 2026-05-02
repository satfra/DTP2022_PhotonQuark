#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

using vec_cmplx = std::vector<std::complex<double>>;
using mat_cmplx = std::vector<vec_cmplx>;

using vec_double = std::vector<double>;
using mat_double = std::vector<vec_double>;

// Nested-vector aliases that remain off the hot path (e.g. the small DSE
// angular matrix in quark_dse.hh). Hot tensors use Tensor3/Tensor6 below.
using tens_double = std::vector<mat_double>;

namespace tensor_detail {

// Lightweight 1-D row view over flat tensor storage. Behaves like a slice of
// a vector for read/write: operator[], size(), data(), begin/end. Trivially
// copyable so passing by value is cheap.
template<typename T>
class Row1D {
  T* p_;
  std::size_t n_;
public:
  Row1D(T* p, std::size_t n) : p_(p), n_(n) {}
  T& operator[](std::size_t k) const { return p_[k]; }
  std::size_t size() const { return n_; }
  T* data() const { return p_; }
  T* begin() const { return p_; }
  T* end() const { return p_ + n_; }
};

// Lightweight 2-D row view; chained operator[] yields a Row1D.
template<typename T>
class Row2D {
  T* p_;
  std::size_t d1_, d2_;
public:
  Row2D(T* p, std::size_t d1, std::size_t d2) : p_(p), d1_(d1), d2_(d2) {}
  Row1D<T> operator[](std::size_t j) const {
    return Row1D<T>(p_ + j * d2_, d2_);
  }
  std::size_t size() const { return d1_; }
  T* data() const { return p_; }
};

} // namespace tensor_detail

// Flat 3-D tensor with row-major (i, j, k) layout. Element access via
// operator() or via chained operator[] returning views (a[i][j][k] still
// works exactly like the old nested vector). Single contiguous allocation,
// so the inner k-loop walks adjacent doubles instead of pointer-chasing.
template<typename T>
class Tensor3 {
  std::vector<T> data_;
  std::size_t d0_, d1_, d2_;
public:
  Tensor3() : d0_(0), d1_(0), d2_(0) {}
  Tensor3(std::size_t d0, std::size_t d1, std::size_t d2, T fill = T{})
    : data_(d0 * d1 * d2, fill), d0_(d0), d1_(d1), d2_(d2) {}

  T& operator()(std::size_t i, std::size_t j, std::size_t k) {
    return data_[(i * d1_ + j) * d2_ + k];
  }
  const T& operator()(std::size_t i, std::size_t j, std::size_t k) const {
    return data_[(i * d1_ + j) * d2_ + k];
  }

  tensor_detail::Row2D<T> operator[](std::size_t i) {
    return tensor_detail::Row2D<T>(data_.data() + i * d1_ * d2_, d1_, d2_);
  }
  tensor_detail::Row2D<const T> operator[](std::size_t i) const {
    return tensor_detail::Row2D<const T>(data_.data() + i * d1_ * d2_, d1_, d2_);
  }

  std::size_t size() const { return d0_; }
  std::size_t d0() const { return d0_; }
  std::size_t d1() const { return d1_; }
  std::size_t d2() const { return d2_; }
};

// Flat 6-D tensor for K' with layout (i, k, z, j, k', z'). The slice2d
// helper hands out a 2-D (k', z') view that lInterpolator2d can interpolate
// over directly without any pointer-chase through nested heap allocations.
template<typename T>
class Tensor6 {
  std::vector<T> data_;
  std::size_t d0_, d1_, d2_, d3_, d4_, d5_;

  std::size_t flat(std::size_t i, std::size_t k, std::size_t z,
                   std::size_t j, std::size_t kp, std::size_t zp) const {
    return ((((i * d1_ + k) * d2_ + z) * d3_ + j) * d4_ + kp) * d5_ + zp;
  }

public:
  Tensor6() : d0_(0), d1_(0), d2_(0), d3_(0), d4_(0), d5_(0) {}
  Tensor6(std::size_t d0, std::size_t d1, std::size_t d2,
          std::size_t d3, std::size_t d4, std::size_t d5, T fill = T{})
    : data_(d0 * d1 * d2 * d3 * d4 * d5, fill),
      d0_(d0), d1_(d1), d2_(d2), d3_(d3), d4_(d4), d5_(d5) {}

  T& operator()(std::size_t i, std::size_t k, std::size_t z,
                std::size_t j, std::size_t kp, std::size_t zp) {
    return data_[flat(i, k, z, j, kp, zp)];
  }
  const T& operator()(std::size_t i, std::size_t k, std::size_t z,
                      std::size_t j, std::size_t kp, std::size_t zp) const {
    return data_[flat(i, k, z, j, kp, zp)];
  }

  // 2-D view of K'_{ij}(k, z; ·, ·) for (k', z') interpolation.
  tensor_detail::Row2D<T> slice2d(std::size_t i, std::size_t k,
                                  std::size_t z, std::size_t j) {
    const std::size_t off = flat(i, k, z, j, 0, 0);
    return tensor_detail::Row2D<T>(data_.data() + off, d4_, d5_);
  }
  tensor_detail::Row2D<const T> slice2d(std::size_t i, std::size_t k,
                                        std::size_t z, std::size_t j) const {
    const std::size_t off = flat(i, k, z, j, 0, 0);
    return tensor_detail::Row2D<const T>(data_.data() + off, d4_, d5_);
  }
};

// Sparse 6-D tensor — same logical shape as Tensor6 but only allocates
// storage for (i,j) pairs flagged non-zero by the predicate passed at
// construction. The kernel matrix K_ij has only 26 non-zero (i,j) pairs
// out of 144 (8+4 block decoupling + within-block sparsity, see
// Kernels_K.hh::isZeroIndex), so storage drops by ~5.5×.
//
// Behaviour:
//  - operator()(i,k,z,j,kp,zp): defined only for non-zero (i,j); returns a
//    reference into the compressed buffer. Calling for a zero (i,j) is UB
//    in production builds (asserted in debug). The BSE iteration always
//    short-circuits via `if (K::isZeroIndex(i,j)) continue;` before access.
//  - slice2d(i,k,z,j): same restriction, returns a Row2D view.
template<typename T>
class SparseTensor6 {
  std::vector<T> data_;
  std::array<int, 144> ij_to_slot_{};   // -1 means zero
  std::size_t n_slots_ = 0;
  std::size_t d1_ = 0, d2_ = 0, d4_ = 0, d5_ = 0;

  std::size_t flat(std::size_t i, std::size_t k, std::size_t z,
                   std::size_t j, std::size_t kp, std::size_t zp) const {
    const int s = ij_to_slot_[i * 12 + j];
    return ((((static_cast<std::size_t>(s) * d1_) + k) * d2_ + z) * d4_ + kp) * d5_ + zp;
  }

public:
  SparseTensor6() { ij_to_slot_.fill(-1); }

  // Construct with shape (k, z, kp, zp); is_zero(i,j) returns true for
  // entries that should not be stored.
  template<typename ZeroPredicate>
  SparseTensor6(std::size_t k, std::size_t z, std::size_t kp, std::size_t zp,
                ZeroPredicate is_zero, T fill = T{})
    : d1_(k), d2_(z), d4_(kp), d5_(zp)
  {
    ij_to_slot_.fill(-1);
    int slot = 0;
    for (unsigned i = 0; i < 12; ++i)
      for (unsigned j = 0; j < 12; ++j)
        if (!is_zero(i, j))
          ij_to_slot_[i * 12 + j] = slot++;
    n_slots_ = static_cast<std::size_t>(slot);
    data_.assign(n_slots_ * d1_ * d2_ * d4_ * d5_, fill);
  }

  T& operator()(std::size_t i, std::size_t k, std::size_t z,
                std::size_t j, std::size_t kp, std::size_t zp) {
    return data_[flat(i, k, z, j, kp, zp)];
  }
  const T& operator()(std::size_t i, std::size_t k, std::size_t z,
                      std::size_t j, std::size_t kp, std::size_t zp) const {
    return data_[flat(i, k, z, j, kp, zp)];
  }

  tensor_detail::Row2D<T> slice2d(std::size_t i, std::size_t k,
                                  std::size_t z, std::size_t j) {
    const std::size_t off = flat(i, k, z, j, 0, 0);
    return tensor_detail::Row2D<T>(data_.data() + off, d4_, d5_);
  }
  tensor_detail::Row2D<const T> slice2d(std::size_t i, std::size_t k,
                                        std::size_t z, std::size_t j) const {
    const std::size_t off = flat(i, k, z, j, 0, 0);
    return tensor_detail::Row2D<const T>(data_.data() + off, d4_, d5_);
  }

  std::size_t bytes() const { return data_.size() * sizeof(T); }
  std::size_t n_slots() const { return n_slots_; }
};

using tens_cmplx = Tensor3<std::complex<double>>;
using ijtens2_double = SparseTensor6<double>;

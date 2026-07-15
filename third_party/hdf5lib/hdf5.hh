#pragma once

/// hdf5lib — header-only C++ wrapper around the HDF5 C API.
///
/// A narrow, purpose-built wrapper (vendored from the DiFfRG project and
/// re-namespaced). Everything lives under `hdf5`. See the individual headers
/// for details. The umbrella include below pulls in all public types.

#include <hdf5lib/dataset.hh>
#include <hdf5lib/dataspace.hh>
#include <hdf5lib/datatype.hh>
#include <hdf5lib/file.hh>
#include <hdf5lib/group.hh>
#include <hdf5lib/handle.hh>

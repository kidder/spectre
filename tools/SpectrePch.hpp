// Distributed under the MIT License.
// See LICENSE.txt for details.

#ifndef SPECTRE_PCH_HPP
#define SPECTRE_PCH_HPP

// Include STL headers that are included 100+ times or show up as expensive
// headers with ClangBuildAnalyzer
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <optional>
#include <ostream>
#include <random>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// For Blaze error handling (see SetupBlaze.cmake)
#include <Utilities/BlazeExceptions.hpp>

#include <Utilities/ErrorHandling/Assert.hpp>
#include <blaze/math/CustomVector.h>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/DynamicMatrix.h>
#include <blaze/math/DenseVector.h>
#include <blaze/math/GroupTag.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Version.h>
#include <blaze/util/typetraits/RemoveConst.h>

// Include Brigand related headers
#include <Utilities/TMPL.hpp>

#include <charm++.h>
#include <pup_stl.h>

#endif  // SPECTRE_PCH_HPP

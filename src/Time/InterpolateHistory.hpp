// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "NumericalAlgorithms/Interpolation/LagrangePolynomial.hpp"
#include "Time/History.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ErrorHandling/Assert.hpp"
#include "Utilities/Gsl.hpp"

namespace TimeSteppers {

/// \brief Interpolate the History `source_history` to the set of TimeStepId%s
/// `target_time_step_ids`
///
/// \details Uses Lagrange polynomials at the full set of source times to
/// interpolate to the target times.
///
/// \note Assumes that History has values at all times.  If this is not true,
/// a std::exception will be thrown.  Also all the target TimeStepId%s should be
/// at the start of a step
template <typename Vars>
History<Vars> interpolate_history(
    const History<Vars>& source_history,
    const std::vector<TimeStepId>& target_time_step_ids) {
  History<Vars> result{source_history.integration_order()};
  using DerivVars = db::prefix_variables<::Tags::dt, Vars>;
  std::vector<double> source_times(source_history.size());
  alg::transform(source_history, source_times.begin(), [](const auto& record) {
    return record.time_step_id.substep_time();
  });
  for (const auto& target_time_step_id : target_time_step_ids) {
    ASSERT(target_time_step_id.substep() == 0, "Bad target substep");
    const double target_time = target_time_step_id.substep_time();
    std::vector<double> coefficients(source_times.size());
    for (size_t i = 0; i < source_times.size(); ++i) {
      coefficients[i] = lagrange_polynomial(
          i, target_time, source_times.begin(), source_times.end());
    }
    Vars value = coefficients[0] * source_history[0].value.value();
    DerivVars derivative = coefficients[0] * source_history[0].derivative;
    for (size_t i = 1; i < source_times.size(); ++i) {
      value += coefficients[i] * source_history[i].value.value();
      derivative += coefficients[i] * source_history[i].derivative;
    }
    result.insert(target_time_step_id, std::move(value), std::move(derivative));
  }
  return result;
}
}  // namespace TimeSteppers

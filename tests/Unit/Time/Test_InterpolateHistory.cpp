// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Helpers//DataStructures/TestTags.hpp"
#include "Time/History.hpp"
#include "Time/InterpolateHistory.hpp"
#include "Time/Slab.hpp"
#include "Time/Time.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/TMPL.hpp"

namespace {
void test_interpolation() {
  const auto f = [](const double t) { return 0.5 + 1.5 * t + 2.5 * square(t); };
  const auto df_dt = [](const double t) { return 1.5 + 5.0 * t; };
  DataVector x{1.5, 2.75, 3.25};
  const auto v = [&f, &x](const double t) {
    return f(t) * (x + 2.0 * square(x) - 4.0);
  };
  const auto dv_dt = [&df_dt, &x](const double t) {
    return df_dt(t) * (x + 2.0 * square(x) - 4.0);
  };

  using VariablesType =
      Variables<tmpl::list<TestHelpers::Tags::Scalar<DataVector>>>;

  using DtVariablesType =
      Variables<tmpl::list<::Tags::dt<TestHelpers::Tags::Scalar<DataVector>>>>;

  const auto add_to_history = [&f, &df_dt](
                                  TimeSteppers::History<double>& the_history,
                                  const Time& time) {
    the_history.insert(TimeStepId(true, 0, time), f(time.value()),
                       df_dt(time.value()));
  };

  const auto add_to_history_vars =
      [&v, &dv_dt, &x](TimeSteppers::History<VariablesType>& the_history,
                       const Time& time) {
        VariablesType vars{x.size()};
        DtVariablesType dt_vars{x.size()};
        get(get<TestHelpers::Tags::Scalar<DataVector>>(vars)) = v(time.value());
        get(get<::Tags::dt<TestHelpers::Tags::Scalar<DataVector>>>(dt_vars)) =
            dv_dt(time.value());
        the_history.insert(TimeStepId(true, 0, time), std::move(vars),
                           std::move(dt_vars));
      };

  const Slab slab(1., 3.);
  TimeSteppers::History<double> initial_history{4};
  TimeSteppers::History<double> expected_history{4};
  TimeSteppers::History<VariablesType> initial_history_vars{4};
  TimeSteppers::History<VariablesType> expected_history_vars{4};
  for (int i = 0; i < 4; ++i) {
    add_to_history(initial_history, slab.start() + slab.duration() * i / 4);
    add_to_history(expected_history,
                   slab.start() + slab.duration() * (4 + i) / 8);
    add_to_history_vars(initial_history_vars,
                        slab.start() + slab.duration() * i / 4);
    add_to_history_vars(expected_history_vars,
                        slab.start() + slab.duration() * (4 + i) / 8);
  }

  std::vector<TimeStepId> target_time_step_ids(4);
  alg::transform(expected_history, target_time_step_ids.begin(),
                 [](const auto& record) { return record.time_step_id; });

  const auto interpolated_history =
      TimeSteppers::interpolate_history(initial_history, target_time_step_ids);
  CHECK(interpolated_history == expected_history);

  const auto interpolated_history_vars = TimeSteppers::interpolate_history(
      initial_history_vars, target_time_step_ids);
  CHECK(interpolated_history_vars == expected_history_vars);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Time.InterpolateHistory", "[Unit][Time]") {
  test_interpolation();
}

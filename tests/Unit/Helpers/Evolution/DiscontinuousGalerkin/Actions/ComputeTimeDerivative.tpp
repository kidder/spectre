// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <memory>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/BoundaryConditions/BoundaryCondition.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/Initialization/ConservativeSystem.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/BoundaryCondition.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/BoundaryCorrection.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Component.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ConservativeSystem.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/DomainCreator.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/SystemType.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "Time/AdvanceTime.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/Tags/StepperErrorEstimatesEnabled.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/TimeSteppers/AdamsBashforth.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

#include "Parallel/Printf/Printf.hpp"

namespace TestHelpers::evolution::dg::Actions {
namespace {
template <typename Component, typename Runner, size_t Dim>
void print_databox_tags(const Runner& runner, const ElementId<Dim>& id) {
  const auto& box = get_databox<Component>(runner, id);
  Parallel::printf("%s\n", box.print_tags());
}

template <typename Component, typename Runner, size_t Dim>
void print_databox_items(const Runner& runner, const ElementId<Dim>& id) {
  const auto& box = get_databox<Component>(runner, id);
  Parallel::printf("%s\n", box.print_items());
}

template <typename Metavariables>
ActionTesting::MockRuntimeSystem<Metavariables> make_runner(
    const DomainCreator<Metavariables::volume_dim>& domain_creator) {
  static constexpr size_t volume_dim = Metavariables::volume_dim;
  std::unique_ptr<TimeStepper> time_stepper =
      std::make_unique<TimeSteppers::AdamsBashforth>(5);
  return ActionTesting::MockRuntimeSystem<Metavariables>{
      {std::move(time_stepper), std::move(domain_creator.create_domain()),
       ::dg::Formulation::StrongInertial,
       std::make_unique<BoundaryCorrection<volume_dim>>(),
       domain_creator.external_boundary_conditions()}};
}

template <size_t Dim, typename System>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;
  using system = System;
  static constexpr bool local_time_stepping = false;

  using initialization_actions = tmpl::list<
      ActionTesting::InitializeDataBox<
          tmpl::list<domain::Tags::InitialExtents<volume_dim>,
                     domain::Tags::InitialRefinementLevels<volume_dim>,
                     ::evolution::dg::Tags::Quadrature, ::Tags::Time,
                     Initialization::Tags::InitialTimeDelta,
                     Initialization::Tags::InitialSlabSize<local_time_stepping>,
                     ::Tags::StepperErrorEstimatesEnabled, Var3>,
          tmpl::list<>>,
      Initialization::Actions::InitializeItems<
          Initialization::TimeStepping<Metavariables, TimeStepper>,
          ::evolution::dg::Initialization::Domain<Metavariables>,
          Initialization::TimeStepperHistory<Metavariables>>,
      Initialization::Actions::ConservativeSystem<System>,
      InitializeVars<system>,
      ::evolution::dg::Initialization::Mortars<volume_dim, System>,
      ::Actions::MutateApply<AdvanceTime<>>>;
  using test_actions =
      tmpl::list<::evolution::dg::Actions::ComputeTimeDerivative<
          Dim, System, AllStepChoosers, local_time_stepping, false>>;
  using component_list = tmpl::list<
      Component<Metavariables, initialization_actions, test_actions>>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<domain::BoundaryConditions::BoundaryCondition,
                             tmpl::list<BoundaryCondition<Dim>>>,
                  tmpl::pair<::evolution::BoundaryCorrection,
                             tmpl::list<BoundaryCorrection<Dim>>>>;
  };
};

template <SystemType system_type, bool HasPrims, size_t Dim>
void test_impl() {
  CAPTURE(system_type);
  CAPTURE(HasPrims);
  CAPTURE(Dim);

  domain::creators::register_derived_with_charm();
  const auto creator = domain_creator<Dim>();
  const auto domain = creator->create_domain();
  const size_t num_blocks = domain.blocks().size();

  using system = ConservativeSystem<Dim>;

  using metavariables = Metavariables<Dim, system>;
  register_classes_with_charm<TimeSteppers::AdamsBashforth>();
  register_factory_classes_with_charm<metavariables>();

  auto runner = make_runner<metavariables>(*creator);

  using component =
      Component<metavariables, typename metavariables::initialization_actions,
                typename metavariables::test_actions>;

  const auto initial_refinement = creator->initial_refinement_levels();
  const auto initial_extents = creator->initial_extents();
  const double initial_time = 1.5;
  const double initial_dt = 0.5;
  const double initial_slab_size = initial_dt;

  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::emplace_component_and_initialize<component>(
        &runner, element_id,
        {initial_extents, initial_refinement,
         Spectral::Quadrature::GaussLobatto, initial_time, initial_dt,
         initial_slab_size, false, Scalar<DataVector>{}});
    for (size_t i = 0; i < 5; ++i) {
      ActionTesting::next_action<component>(make_not_null(&runner), element_id);
    }
  }

  BoundaryCondition<Dim>::number_of_times_called = 0;
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  using variables_tag = typename system::variables_tag;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;

  const auto& first_box = get_databox<component>(runner, ElementId<Dim>{0});
  print_databox_tags<component>(runner, ElementId<Dim>{0});
  using tags_for_uncomparable_items =
      tmpl::list<domain::Tags::ElementMap<Dim, Frame::Grid>,
                 domain::CoordinateMaps::Tags::CoordinateMap<Dim, Frame::Grid,
                                                             Frame::Inertial>,
                 Parallel::Tags::MetavariablesImpl<metavariables>>;
  using tags_for_unchanged_items =
      tmpl::list<Parallel::Tags::ArrayIndex<ElementId<Dim>>,
                 Parallel::Tags::GlobalCache<metavariables>,
                 domain::Tags::InitialExtents<Dim>,
                 domain::Tags::InitialRefinementLevels<Dim>,
                 ::evolution::dg::Tags::Quadrature, Tags::Time,
                 Initialization::Tags::InitialTimeDelta,
                 Initialization::Tags::InitialSlabSize<false>,
                 Tags::StepperErrorEstimatesEnabled, Var3, Tags::TimeStepId,
                 Tags::StepNumberWithinSlab, Tags::AdaptiveSteppingDiagnostics,
                 Tags::Next<Tags::TimeStepId>, Tags::TimeStep,
                 Tags::ChangeSlabSize::SlabSizeGoal, domain::Tags::Mesh<Dim>,
                 domain::Tags::Element<Dim>, domain::Tags::NeighborMesh<Dim>,
                 variables_tag, Tags::HistoryEvolvedVariables<variables_tag>,
                 ::evolution::dg::Tags::MortarData<Dim>,
                 ::evolution::dg::Tags::MortarMesh<Dim>,
                 ::evolution::dg::Tags::MortarInfo<Dim>,
                 ::evolution::dg::Tags::MortarNextTemporalId<Dim>>;
  using tags_to_check =
      tmpl::list<dt_variables_tag,
                 ::evolution::dg::Tags::NormalCovectorAndMagnitude<Dim>,
                 ::evolution::dg::Tags::MortarDataHistory<
                     Dim, typename dt_variables_tag::type>>;
  using tested_tags = tmpl::append<tags_for_uncomparable_items,
                                   tags_for_unchanged_items, tags_to_check>;
  using box_type = std::decay_t<decltype(first_box)>;
  using mutable_tags = box_type::mutable_item_creation_tags;
  using untested_tags = tmpl::list_difference<mutable_tags, tested_tags>;
  static_assert(std::is_same_v<untested_tags, tmpl::list<>>,
                "Please put each tag in untested_tags in the appropriate list "
                "for tested_tags");
  const auto items_before = copy_items<tags_for_unchanged_items>(first_box);
  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::next_action<component>(make_not_null(&runner), element_id);
  }
  const auto items_after = copy_items<tags_for_unchanged_items>(first_box);
  CHECK(items_before == items_after);
  CHECK(BoundaryCondition<Dim>::number_of_times_called == num_blocks);
  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    //    print_databox_items<component>(runner, element_id);
    const auto& box = get_databox<component>(runner, element_id);
    const auto& x =
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
    const DataVector r_squared = get(dot_product(x, x));
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    Variables<tmpl::list<::Tags::dt<Var1>, ::Tags::dt<Var2<Dim>>>>
        expected_dt_evolved_vars{mesh.number_of_grid_points()};
    Var1::dt_value(
        make_not_null(&get<Tags::dt<Var1>>(expected_dt_evolved_vars)),
        r_squared, x);
    Var2<Dim>::dt_value(
        make_not_null(&get<Tags::dt<Var2<Dim>>>(expected_dt_evolved_vars)),
        r_squared, x);
    const auto& dt_evolved_vars =
        ActionTesting::get_databox_tag<component, dt_variables_tag>(runner,
                                                                    element_id);
    const auto approx = Approx::custom().scale(100.0).epsilon(1.e-13);
    CHECK_VARIABLES_CUSTOM_APPROX(dt_evolved_vars, expected_dt_evolved_vars,
                                  approx);
  }
}
}  // namespace

template <SystemType system_type, bool UsePrims, size_t Dim>
void test() {
  test_impl<system_type, UsePrims, Dim>();
}
}  // namespace TestHelpers::evolution::dg::Actions

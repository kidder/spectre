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
#include "Domain/Creators/TimeDependence/CubicScale.hpp"
#include "Domain/Creators/TimeDependence/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/TimeDependence/UniformTranslation.hpp"
#include "Domain/FunctionsOfTime/RegisterDerivedWithCharm.hpp"
#include "Domain/InterfaceLogicalCoordinates.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ApplyBoundaryCorrections.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarData.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarDataHolder.hpp"
#include "Evolution/DiscontinuousGalerkin/MortarTags.hpp"
#include "Evolution/DiscontinuousGalerkin/NormalVectorTags.hpp"
#include "Evolution/Initialization/ConservativeSystem.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/BoundaryCondition.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Component.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/DomainCreator.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "Time/AdvanceTime.hpp"
#include "Time/StepChoosers/Constant.hpp"
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

template <bool UseMovingMesh, size_t Dim>
void check_dt_evolved_vars(
    const Mesh<Dim>& mesh, const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
    const double t,
    const Variables<tmpl::list<::Tags::dt<Var1>, ::Tags::dt<Var2<Dim>>>>&
        dt_evolved_vars) {
  Variables<tmpl::list<::Tags::dt<Var1>, ::Tags::dt<Var2<Dim>>>>
      expected_dt_evolved_vars{mesh.number_of_grid_points(), 0.0};
  Var1::dt_value<UseMovingMesh>(
      make_not_null(&get<Tags::dt<Var1>>(expected_dt_evolved_vars)), x, t);
  Var2<Dim>::template dt_value<UseMovingMesh>(
      make_not_null(&get<Tags::dt<Var2<Dim>>>(expected_dt_evolved_vars)), x, t);
  const auto approx = Approx::custom().scale(1.0).epsilon(1.e-11);
  CHECK_VARIABLES_CUSTOM_APPROX(dt_evolved_vars, expected_dt_evolved_vars,
                                approx);
}

template <size_t Dim>
tnsr::I<DataVector, Dim, Frame::Inertial> inertial_interface_coordinates(
    const Mesh<Dim - 1>& mesh, const Direction<Dim>& direction,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  return coordinate_map(
      element_map(interface_logical_coordinates(mesh, direction)), t,
      functions_of_time);
}

template <typename System, bool LocalTimeStepping, bool UseMovingMesh,
          typename HistoryMap, typename NormalsMap, size_t Dim>
void check_mortar_data_and_inboxes(
    const DirectionalIdMap<Dim, ::evolution::dg::MortarDataHolder<Dim>>&
        mortar_data,
    const HistoryMap& mortar_data_history,
    const DirectionalIdMap<Dim, ::evolution::dg::BoundaryData<Dim>>&
        boundary_data,
    const NormalsMap& normal_covector_and_magnitude,
    const Mesh<Dim>& volume_mesh,
    const DirectionalIdMap<Dim, Mesh<Dim - 1>>& mortar_meshes,
    const DirectionalIdMap<Dim, Mesh<Dim>>& neighbor_meshes,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const Scalar<DataVector>& det_inv_jacobian, const TimeStepId& time_step_id,
    const TimeStepId& next_time_step_id, const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time,
    const ::dg::Formulation dg_formulation) {
  if constexpr (not LocalTimeStepping) {
    REQUIRE(mortar_data_history.size() == 0);
  }
  const double sign =
      dg_formulation == ::dg::Formulation::StrongInertial ? 1.0 : -1.0;
  for (const auto& [mortar_id, mortar_data_holder] : mortar_data) {
    CAPTURE(mortar_id);
    const auto& mortar_mesh = mortar_meshes.at(mortar_id);
    const auto x = inertial_interface_coordinates(
        mortar_mesh, mortar_id.direction(), element_map, coordinate_map, t,
        functions_of_time);
    const auto& normal_covector =
        get<::evolution::dg::Tags::NormalCovector<Dim>>(
            normal_covector_and_magnitude.at(mortar_id.direction()).value());
    auto expected_mortar_data_vars =
        System::boundary_correction::template expected_mortar_data<
            UseMovingMesh>(normal_covector, x, t, sign);
    DataVector expected_mortar_data(expected_mortar_data_vars.data(),
                                    expected_mortar_data_vars.size());

    const auto& local_mortar_data = mortar_data_holder.local();
    CAPTURE(local_mortar_data);
    if constexpr (LocalTimeStepping) {
      const size_t sliced_dim = mortar_id.direction().dimension();
      const Spectral::Quadrature sliced_quadrature =
          volume_mesh.quadrature(sliced_dim);
      REQUIRE(not local_mortar_data.mortar_data.has_value());
      REQUIRE(not local_mortar_data.mortar_mesh.has_value());
      REQUIRE(not local_mortar_data.face_mesh.has_value());
      const auto& boundary_history = mortar_data_history.at(mortar_id);
      const auto& lts_mortar_data = boundary_history.local().data(time_step_id);
      CAPTURE(lts_mortar_data);
      REQUIRE(lts_mortar_data.mortar_data.has_value());
      CHECK_ITERABLE_APPROX(lts_mortar_data.mortar_data.value(),
                            expected_mortar_data);
      REQUIRE(lts_mortar_data.mortar_mesh.has_value());
      CHECK(lts_mortar_data.mortar_mesh.value() == mortar_mesh);
      REQUIRE(lts_mortar_data.face_mesh.has_value());
      CHECK(lts_mortar_data.face_mesh.value() ==
            volume_mesh.slice_away(mortar_id.direction().dimension()));
      REQUIRE(lts_mortar_data.face_normal_magnitude.has_value());
      const auto& expected_magnitude =
          get<::evolution::dg::Tags::MagnitudeOfNormal>(
              normal_covector_and_magnitude.at(mortar_id.direction()).value());
      CHECK(lts_mortar_data.face_normal_magnitude == expected_magnitude);
      if (sliced_quadrature == Spectral::Quadrature::Gauss) {
        REQUIRE(lts_mortar_data.face_det_jacobian.has_value());
        REQUIRE(lts_mortar_data.volume_det_inv_jacobian.has_value());
        CHECK_ITERABLE_APPROX(lts_mortar_data.volume_det_inv_jacobian.value(),
                              det_inv_jacobian);
        REQUIRE(lts_mortar_data.volume_mesh.has_value());
        CHECK(lts_mortar_data.volume_mesh.value() == volume_mesh);
      } else {
        REQUIRE(not lts_mortar_data.face_det_jacobian.has_value());
        REQUIRE(not lts_mortar_data.volume_det_inv_jacobian.has_value());
        REQUIRE(not lts_mortar_data.volume_mesh.has_value());
      }
    } else {
      REQUIRE(local_mortar_data.mortar_data.has_value());
      CHECK_ITERABLE_APPROX(local_mortar_data.mortar_data.value(),
                            expected_mortar_data);
      REQUIRE(local_mortar_data.mortar_mesh.has_value());
      CHECK(local_mortar_data.mortar_mesh.value() == mortar_mesh);
      REQUIRE(local_mortar_data.face_mesh.has_value());
      CHECK(local_mortar_data.face_mesh.value() ==
            volume_mesh.slice_away(mortar_id.direction().dimension()));
    }
    REQUIRE(not local_mortar_data.face_normal_magnitude.has_value());
    REQUIRE(not local_mortar_data.face_det_jacobian.has_value());
    REQUIRE(not local_mortar_data.volume_det_inv_jacobian.has_value());
    REQUIRE(not local_mortar_data.volume_mesh.has_value());

    const ::evolution::dg::BoundaryData<Dim>& received_data =
        boundary_data.at(mortar_id);
    auto expected_received_data_vars =
        System::boundary_correction::template expected_mortar_data<
            UseMovingMesh>(normal_covector, x, t, -sign);
    DataVector expected_received_data(expected_received_data_vars.data(),
                                      expected_received_data_vars.size());
    CAPTURE(received_data);
    CHECK(received_data.volume_mesh == neighbor_meshes.at(mortar_id));
    REQUIRE(not received_data.volume_mesh_ghost_cell_data.has_value());
    REQUIRE(received_data.boundary_correction_mesh.has_value());
    CHECK(received_data.boundary_correction_mesh.value() == mortar_mesh);
    REQUIRE(not received_data.ghost_cell_data.has_value());
    REQUIRE(received_data.boundary_correction_data.has_value());
    CHECK_ITERABLE_APPROX(received_data.boundary_correction_data.value(),
                          expected_received_data);
    CHECK(received_data.validity_range == next_time_step_id);
    CHECK(received_data.tci_status == 0);
    CHECK(received_data.integration_order == 1);
    REQUIRE(not received_data.interpolated_boundary_data.has_value());
  }
}

template <typename Metavariables>
ActionTesting::MockRuntimeSystem<Metavariables> make_runner(
    const DomainCreator<Metavariables::volume_dim>& domain_creator,
    const double initial_dt, const ::dg::Formulation dg_formulation,
    const domain::creators::time_dependence::TimeDependence<
        Metavariables::volume_dim>& time_dependence) {
  std::unique_ptr<typename Metavariables::TimeStepperBase> time_stepper =
      std::make_unique<TimeSteppers::AdamsBashforth>(5);
  using boundary_correction =
      typename Metavariables::system::boundary_correction;
  auto functions_of_time = time_dependence.functions_of_time();
  if constexpr (Metavariables::local_time_stepping) {
    std::vector<std::unique_ptr<StepChooser<StepChooserUse::LtsStep>>>
        step_choosers;
    step_choosers.emplace_back(
        std::make_unique<StepChoosers::Constant>(initial_dt));
    const double minimum_time_step = 1.0e-8;
    return ActionTesting::MockRuntimeSystem<Metavariables>{
        {std::move(time_stepper), std::move(domain_creator.create_domain()),
         dg_formulation, std::make_unique<boundary_correction>(),
         domain_creator.external_boundary_conditions(), minimum_time_step,
         std::move(step_choosers)},
        {std::move(functions_of_time)}};
  } else {
    return ActionTesting::MockRuntimeSystem<Metavariables>{
        {std::move(time_stepper), std::move(domain_creator.create_domain()),
         dg_formulation, std::make_unique<boundary_correction>(),
         domain_creator.external_boundary_conditions()},
        {std::move(functions_of_time)}};
  }
}

template <size_t Dim, typename System, bool LocalTimeStepping,
          bool UseMovingMesh, bool UseNodegroupDgElements>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;
  static constexpr bool local_time_stepping = LocalTimeStepping;
  using system = System;
  using boundary_correction = typename system::boundary_correction;

  using TimeStepperBase =
      tmpl::conditional_t<LocalTimeStepping, LtsTimeStepper, TimeStepper>;

  using initialization_actions = tmpl::list<
      ActionTesting::InitializeDataBox<
          tmpl::list<domain::Tags::InitialExtents<Dim>,
                     domain::Tags::InitialRefinementLevels<Dim>,
                     ::evolution::dg::Tags::Quadrature, ::Tags::Time,
                     Initialization::Tags::InitialTimeDelta,
                     Initialization::Tags::InitialSlabSize<LocalTimeStepping>,
                     ::Tags::StepperErrorEstimatesEnabled, Source1,
                     Source2<Dim>>,
          tmpl::list<>>,
      Initialization::Actions::InitializeItems<
          Initialization::TimeStepping<Metavariables, TimeStepperBase>,
          ::evolution::dg::Initialization::Domain<Metavariables>,
          Initialization::TimeStepperHistory<Metavariables>>,
      Initialization::Actions::ConservativeSystem<System>,
      InitializeVars<System>,
      ::evolution::dg::Initialization::Mortars<Dim, System>,
      ::Actions::MutateApply<AdvanceTime<>>>;

  using test_actions = tmpl::list<
      ::evolution::dg::Actions::ComputeTimeDerivative<
          Dim, System, AllStepChoosers, LocalTimeStepping,
          UseNodegroupDgElements>,
      ::evolution::dg::Actions::ApplyBoundaryCorrectionsToTimeDerivative<
          Dim, UseNodegroupDgElements>>;

  using component_list = tmpl::list<Component<Metavariables>>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<domain::BoundaryConditions::BoundaryCondition,
                             tmpl::list<BoundaryCondition<Dim>>>,
                  tmpl::pair<::evolution::BoundaryCorrection,
                             tmpl::list<boundary_correction>>,
                  tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                             tmpl::list<StepChoosers::Constant>>>;
  };
};

template <size_t Dim, typename System, bool LocalTimeStepping,
          bool UseMovingMesh, bool UseNodegroupDgElements>
void test_impl(const ::dg::Formulation dg_formulation,
               const Spectral::Quadrature quadrature) {
  CAPTURE(Dim);
  CAPTURE(LocalTimeStepping);
  CAPTURE(UseMovingMesh);
  CAPTURE(dg_formulation);
  CAPTURE(quadrature);

  domain::creators::register_derived_with_charm();
  domain::creators::time_dependence::register_derived_with_charm();
  domain::FunctionsOfTime::register_derived_with_charm();

  std::unique_ptr<domain::creators::time_dependence::TimeDependence<Dim>>
      time_dependence{nullptr};
  if constexpr (UseMovingMesh) {
    time_dependence =
        std::make_unique<domain::creators::time_dependence::CubicScale<Dim>>(
            0.0, 1000.0, true, std::array{grid::a_0, grid::a_0},
            std::array{grid::a_dot, grid::a_dot}, std::array{0.0, 0.0});
  } else {
    time_dependence =
        std::make_unique<domain::creators::time_dependence::None<Dim>>();
  }

  const auto creator = domain_creator(*time_dependence);
  const auto domain = creator->create_domain();
  const size_t num_blocks = domain.blocks().size();

  using metavariables = Metavariables<Dim, System, LocalTimeStepping,
                                      UseMovingMesh, UseNodegroupDgElements>;
  register_classes_with_charm<TimeSteppers::AdamsBashforth>();
  register_factory_classes_with_charm<metavariables>();

  const double initial_dt = 0.5;

  auto runner = make_runner<metavariables>(*creator, initial_dt, dg_formulation,
                                           *time_dependence);

  using component = Component<metavariables>;

  const auto initial_refinement = creator->initial_refinement_levels();
  const auto initial_extents = creator->initial_extents();
  const double initial_time = 1.5;
  const double initial_slab_size = initial_dt;

  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::emplace_component_and_initialize<component>(
        &runner, element_id,
        {initial_extents, initial_refinement, quadrature, initial_time,
         initial_dt, initial_slab_size, false, Scalar<DataVector>{},
         tnsr::I<DataVector, Dim, Frame::Inertial>{}});
    for (size_t i = 0; i < 5; ++i) {
      ActionTesting::next_action<component>(make_not_null(&runner), element_id);
    }
  }

  BoundaryCondition<Dim>::number_of_times_called = 0;
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  using variables_tag = typename System::variables_tag;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;

  const auto& first_box = get_databox<component>(runner, ElementId<Dim>{0});
  print_databox_tags<component>(runner, ElementId<Dim>{0});
  using tags_for_uncomparable_items =
      tmpl::list<domain::Tags::ElementMap<Dim, Frame::Grid>,
                 domain::CoordinateMaps::Tags::CoordinateMap<Dim, Frame::Grid,
                                                             Frame::Inertial>,
                 Parallel::Tags::MetavariablesImpl<metavariables>>;
  using tags_for_unchanged_items = tmpl::append<
      tmpl::list<Parallel::Tags::ArrayIndex<ElementId<Dim>>,
                 Parallel::Tags::GlobalCache<metavariables>,
                 domain::Tags::InitialExtents<Dim>,
                 domain::Tags::InitialRefinementLevels<Dim>,
                 ::evolution::dg::Tags::Quadrature, Tags::Time,
                 Initialization::Tags::InitialTimeDelta,
                 Initialization::Tags::InitialSlabSize<LocalTimeStepping>,
                 Tags::StepperErrorEstimatesEnabled, Source1, Source2<Dim>,
                 Tags::TimeStepId, Tags::StepNumberWithinSlab,
                 Tags::AdaptiveSteppingDiagnostics,
                 Tags::Next<Tags::TimeStepId>, Tags::TimeStep,
                 Tags::ChangeSlabSize::SlabSizeGoal, domain::Tags::Mesh<Dim>,
                 domain::Tags::Element<Dim>, domain::Tags::NeighborMesh<Dim>,
                 variables_tag, Tags::HistoryEvolvedVariables<variables_tag>,
                 ::evolution::dg::Tags::MortarMesh<Dim>,
                 ::evolution::dg::Tags::MortarInfo<Dim>,
                 ::evolution::dg::Tags::MortarNextTemporalId<Dim>>,
      tmpl::conditional_t<System::has_primitive_and_conservative_vars,
                          tmpl::list<Tags::Variables<tmpl::list<PrimVar>>>,
                          tmpl::list<>>>;
  using tags_to_check =
      tmpl::list<dt_variables_tag,
                 ::evolution::dg::Tags::NormalCovectorAndMagnitude<Dim>,
                 ::evolution::dg::Tags::MortarData<Dim>,
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
  size_t number_of_element_external_boundaries = 0;
  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::next_action<component>(make_not_null(&runner), element_id);
    const auto& box = get_databox<component>(runner, element_id);
    const Element<Dim>& element = db::get<domain::Tags::Element<Dim>>(box);
    number_of_element_external_boundaries +=
        element.external_boundaries().size();
  }
  const auto items_after = copy_items<tags_for_unchanged_items>(first_box);
  CHECK(items_before == items_after);
  CHECK(BoundaryCondition<Dim>::number_of_times_called ==
        number_of_element_external_boundaries);

  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    //    print_databox_items<component>(runner, element_id);
    const auto& box = get_databox<component>(runner, element_id);
    const auto& x =
        db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
    const double t = db::get<Tags::Time>(box);
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    if (dg_formulation == ::dg::Formulation::StrongInertial) {
      const auto& dt_evolved_vars =
          ActionTesting::get_databox_tag<component, dt_variables_tag>(
              runner, element_id);
      check_dt_evolved_vars<UseMovingMesh>(mesh, x, t, dt_evolved_vars);
    }
    const auto& mortar_data =
        db::get<::evolution::dg::Tags::MortarData<Dim>>(box);
    const auto& mortar_data_history =
        db::get<::evolution::dg::Tags::MortarDataHistory<
            Dim, typename dt_variables_tag::type>>(box);
    const auto& time_step_id = db::get<Tags::TimeStepId>(box);
    const auto& boundary_data =
        ActionTesting::get_inbox_tag<
            component,
            ::evolution::dg::Tags::BoundaryCorrectionAndGhostCellsInbox<
                Dim, UseNodegroupDgElements>>(runner, element_id)
            .messages.at(time_step_id);
    const auto& normal_covector_and_magnitude =
        db::get<::evolution::dg::Tags::NormalCovectorAndMagnitude<Dim>>(box);
    const auto& mortar_meshes =
        db::get<::evolution::dg::Tags::MortarMesh<Dim>>(box);
    const auto& neighbor_meshes = db::get<domain::Tags::NeighborMesh<Dim>>(box);
    const auto& element_map =
        db::get<domain::Tags::ElementMap<Dim, Frame::Grid>>(box);
    const auto& coordinate_map =
        db::get<domain::CoordinateMaps::Tags::CoordinateMap<Dim, Frame::Grid,
                                                            Frame::Inertial>>(
            box);
    const auto& det_inv_jacobian = db::get<
        domain::Tags::DetInvJacobian<Frame::ElementLogical, Frame::Inertial>>(
        box);
    const auto& next_time_step_id = db::get<Tags::Next<Tags::TimeStepId>>(box);
    const auto& functions_of_time = db::get<domain::Tags::FunctionsOfTime>(box);
    check_mortar_data_and_inboxes<System, LocalTimeStepping, UseMovingMesh>(
        mortar_data, mortar_data_history, boundary_data,
        normal_covector_and_magnitude, mesh, mortar_meshes, neighbor_meshes,
        element_map, coordinate_map, det_inv_jacobian, time_step_id,
        next_time_step_id, t, functions_of_time, dg_formulation);
  }

  if (dg_formulation == ::dg::Formulation::WeakInertial and
      not LocalTimeStepping) {
    for (size_t b = 0; b < num_blocks; ++b) {
      ElementId<Dim> element_id{b};
      ActionTesting::next_action<component>(make_not_null(&runner), element_id);
    }
    for (size_t b = 0; b < num_blocks; ++b) {
      ElementId<Dim> element_id{b};
      const auto& box = get_databox<component>(runner, element_id);
      const auto& x =
          db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
      const double t = db::get<Tags::Time>(box);
      const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
      const auto& dt_evolved_vars =
          ActionTesting::get_databox_tag<component, dt_variables_tag>(
              runner, element_id);
      check_dt_evolved_vars<UseMovingMesh>(mesh, x, t, dt_evolved_vars);
    }
  }
}

template <size_t Dim, typename System, bool LocalTimeStepping,
          bool UseMovingMesh, bool UseNodegroupDgElements>
void test_dg_formulation(const ::dg::Formulation dg_formulation) {
  test_impl<Dim, System, LocalTimeStepping, UseMovingMesh,
            UseNodegroupDgElements>(dg_formulation,
                                    Spectral::Quadrature::GaussLobatto);
  // test_impl<Dim, System, LocalTimeStepping, UseMovingMesh,
  //           UseNodegroupDgElements>(dg_formulation,
  //                                   Spectral::Quadrature::Gauss);
}

template <size_t Dim, typename System, bool LocalTimeStepping,
          bool UseMovingMesh, bool UseNodegroupDgElements>
void test_use_nodegroup_dg_elements() {
  // test_dg_formulation<Dim, System, LocalTimeStepping, UseMovingMesh,
  //                     UseNodegroupDgElements>(
  //     ::dg::Formulation::StrongInertial);
  test_dg_formulation<Dim, System, LocalTimeStepping, UseMovingMesh,
                      UseNodegroupDgElements>(::dg::Formulation::WeakInertial);
}

template <size_t Dim, typename System, bool LocalTimeStepping,
          bool UseMovingMesh>
void test_moving_mesh() {
  test_use_nodegroup_dg_elements<Dim, System, LocalTimeStepping, UseMovingMesh,
                                 false>();
  // Action testing framework does not support node groups, so we cannot
  // currently test this branch...
  // test_use_nodegroup_dg_elements<Dim, System, LocalTimeStepping,
  // UseMovingMesh, true>();
}

template <size_t Dim, typename System, bool LocalTimeStepping>
void test_lts() {
  //  test_moving_mesh<Dim, System, LocalTimeStepping, false>();
  test_moving_mesh<Dim, System, LocalTimeStepping, true>();
}
}  // namespace

template <size_t Dim, typename System>
void test() {
  test_lts<Dim, System, false>();
  //  test_lts<Dim, System, true>();
}
}  // namespace TestHelpers::evolution::dg::Actions

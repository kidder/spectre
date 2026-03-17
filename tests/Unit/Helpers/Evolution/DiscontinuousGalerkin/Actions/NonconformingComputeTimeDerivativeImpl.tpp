// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/MathWrapper.hpp"
#include "DataStructures/TaggedContainers.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/CreateInitialElement.hpp"
#include "Domain/Creators/NonconformingSphericalShells.hpp"
#include "Domain/Creators/RotatedBricks.hpp"
#include "Domain/Creators/RotatedIntervals.hpp"
#include "Domain/Creators/RotatedRectangles.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Creators/Tags/Domain.hpp"
#include "Domain/Domain.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Evolution/BoundaryCorrection.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivative.hpp"
#include "Evolution/DiscontinuousGalerkin/Actions/VolumeTermsImpl.tpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/Mortars.hpp"
#include "Evolution/DiscontinuousGalerkin/Initialization/QuadratureTag.hpp"
#include "Evolution/Initialization/ConservativeSystem.hpp"
#include "Evolution/Initialization/DgDomain.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/PassVariables.hpp"
#include "Framework/ActionTesting.hpp"
#include "Helpers/DataStructures/MathWrapper.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/ComputeTimeDerivativeImpl.hpp"
#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/SystemType.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Formulation.hpp"
#include "NumericalAlgorithms/Spectral/Basis.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Quadrature.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "ParallelAlgorithms/Actions/InitializeItems.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "Time/AdvanceTime.hpp"
#include "Time/BoundaryHistory.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/StepChoosers/Constant.hpp"
#include "Time/StepChoosers/StepChooser.hpp"
#include "Time/Tags/StepperErrorEstimatesEnabled.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsBashforth.hpp"
#include "Utilities/CloneUniquePtrs.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

#include "Parallel/Printf/Printf.hpp"

namespace TestHelpers::evolution::dg::Actions {
// v1 = r^2
// dx_v1 = 2 x
// dy_v1 = 2 y
// dz_v1 = 2 z
//
// v2_x = 3 x^2
// dx_v2_x = 6 x
// dy_v2_x = 0
// dz_v2_x = 0
//
// v2_y = 2 y^2
// dx_v2_y = 0
// dy_v2_y = 4 y
// dz_v2_z = 0
//
// v2_z = z^2
// dx_v2_z = 0
// dy_v2_z = 0
// dz_v2_z = 2 z
//
// v3 = 3 r^2 - 1
//
// f1_i = v1^2 v2^i
// f2_ij = v1 v2^i v2^j + v1^3 delta_ij
//
// dt_v1 = - d_i f1_i + v3^3
//   = - v1^2 (dx_v2_x + dy_v2_y + dz_v2_z)
//     - 2 v1 (v2_x dx_v1 + v2_y dy_v1 + v2_z dz_v1) + v3^3
//   = -r^4 (6 x + 4 y + 2 z) - 2 r^2 (6 x^3 + 4 y^3 + 2 z^3) + (3 r^2 - 1)^3
//
// dt_v2_x = - d_i f2_xi
//   = - dx_v1 v2_x v2_x - 2 v1 dx_v2_x v2_x
//     - dy_v1 v2_x v2_y - v1 dy_v2_x v2_y - v1 v2_x dy_v2_y
//     - dz_v1 v2_x v2_z - v1 dz_v2_x v2_z - v1 v2_x dz_v2_z - 3 v1^2 dx_v1
//   = - 18 x^5 - 36 r^2 x^3 - 12 x^2 y^3 - 12 r^2 x^2 y - 6 x^3 z^2
//     - 6 r^2 x^2 z - 6 r^4 x
//
// dt_v2_y = - d_i f2_yi + v3
//   = - dx_v1 v2_y v2_x - v1 dx_v2_y v2_x - v1 v2_y dy_v2_x
//     - dy_v1 v2_y v2_y - 2 v1 dy_v2_y v2_y
//     - dz_v1 v2_y v2_z - v1 dz_v2_y v2_z - v1 v2_x dz_v2_y - 3 v1^2 dy_v1 + v3
//   = -
struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;

  static void value(const gsl::not_null<type*> var1,
                    const DataVector& r_squared) {
    var1->get() = r_squared;
  }
  // dt var1 = -d_i (var1**2 * var2^i) + var3**2
  //   = -r^4 (6 x + 4 y + 2 z) - 2 r^2 (6 x^3 + 4 y^3 + 2 z^3) + (3 r^2 - 1)^3
  template <size_t Dim>
  static void dt_value(const gsl::not_null<type*> dt_var1,
                    const DataVector& r_squared,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x) {
    dt_var1->get() = -6.0*get<0>(x)*square(r_squared)
      - 12.0*cube(get<0>(x))*r_squared + cube(3.0*r_squared-1.0);
    if constexpr (Dim > 1) {
      dt_var1->get() -= 4.0*get<1>(x)*square(r_squared)
        + 8.0*cube(get<1>(x))*r_squared;
    }
    if constexpr (Dim > 2) {
      dt_var1->get() -= 2.0*get<2>(x)*square(r_squared)
        + 4.0*cube(get<2>(x))*r_squared;
    }
  }
};

template <size_t Dim>
struct Var2 : db::SimpleTag {
  using type = tnsr::I<DataVector, Dim, Frame::Inertial>;

  static void value(const gsl::not_null<type*> var2,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x) {
    get<0>(*var2) = 3.0*square(get<0>(x));
    if constexpr(Dim > 1) {
      get<1>(*var2) = 2.0*square(get<1>(x));
    }
    if constexpr(Dim > 2) {
      get<2>(*var2) = square(get<2>(x));
    }
  }
};

struct Var3 : db::SimpleTag {
  using type = Scalar<DataVector>;

  static void value(const gsl::not_null<type*> var3,
                    const DataVector& r_squared) {
    var3->get() = 3.0*r_squared - 1.0;
  }
};

struct PrimVar1 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template <size_t Dim>
struct PrimVar2 : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim, Frame::Inertial>;
};

struct Var3Squared : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template <size_t Dim, SystemType system_type, bool HasPrimitiveVars>
struct TimeDerivativeTerms {
  using temporary_tags = tmpl::list<Var3Squared>;
  using common_argument_tags = tmpl::list<Var1, Var2<Dim>, Var3>;
  using argument_tags =
      tmpl::conditional_t<HasPrimitiveVars,
                          tmpl::push_back<common_argument_tags, PrimVar1>,
                          common_argument_tags>;

  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Scalar<DataVector>*> dt_var1,
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> dt_var2,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux_var1,
      const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*>
          flux_var2,
      const gsl::not_null<Scalar<DataVector>*> square_var3,
      const Scalar<DataVector>& var1,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
      const Scalar<DataVector>& var3) {
        get(*square_var3) = square(get(var3));

        // Set source terms
        get(*dt_var1) = get(*square_var3);
        for (size_t d = 0; d < Dim; ++d) {
          dt_var2->get(d) = get(var3) * d;
        }

        // Set fluxes
        for (size_t i = 0; i < Dim; ++i) {
          flux_var1->get(i) = square(get(var1)) * var2.get(i);
          for (size_t j = 0; j < Dim; ++j) {
            flux_var2->get(i, j) = var2.get(i) * var2.get(j) * get(var1);
            if (i == j) {
              flux_var2->get(i, j) += cube(get(var1));
            }
          }
        }
        return {true};
      }
};

template <size_t Dim, SystemType system_type, bool HasPrimitiveVars>
struct TimeDerivativeTermsWithVariables
    : public ::evolution::PassVariables,
      private TimeDerivativeTerms<Dim, system_type, HasPrimitiveVars> {
  using base = TimeDerivativeTerms<Dim, system_type, HasPrimitiveVars>;

  using temporary_tags = typename base::temporary_tags;
  using argument_tags = typename base::argument_tags;

  using dt_var1 = ::Tags::dt<Var1>;
  using dt_var2 = ::Tags::dt<Var2<Dim>>;
  using flux_var1 = ::Tags::Flux<Var1, tmpl::size_t<Dim>, Frame::Inertial>;
  using flux_var2 = ::Tags::Flux<Var2<Dim>, tmpl::size_t<Dim>, Frame::Inertial>;

  static ::evolution::dg::TimeDerivativeDecisions<Dim> apply(
      // Time derivatives returned by reference. All the tags in the
      // variables_tag in the system struct.
      const gsl::not_null<Variables<tmpl::list<dt_var1, dt_var2>>*> dt_vars,

      // Fluxes returned by reference. Listed in the system struct as
      // flux_variables.
      const gsl::not_null<Variables<tmpl::list<flux_var1, flux_var2>>*>
          flux_vars,

      // Temporaries returned by reference. Listed in temporary_tags above.
      const gsl::not_null<Variables<tmpl::list<>>*> temporaries) {
    // just forward to other implementation to reduce code duplication
    base::apply(get<dt_var1>(dt_vars), get<dt_var2>(dt_vars),
                get<flux_var1>(flux_vars), get<flux_var2>(flux_vars));
    return {true};
  }
};

template <size_t Dim, bool HasPrims>
struct BoundaryTerms final : public ::evolution::BoundaryCorrection {
  struct MaxAbsCharSpeed : db::SimpleTag {
    using type = Scalar<DataVector>;
  };

  /// \cond
  explicit BoundaryTerms(CkMigrateMessage* /*unused*/) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(BoundaryTerms);  // NOLINT
  /// \endcond
  BoundaryTerms() = default;
  BoundaryTerms(const BoundaryTerms&) = default;
  BoundaryTerms& operator=(const BoundaryTerms&) = default;
  BoundaryTerms(BoundaryTerms&&) = default;
  BoundaryTerms& operator=(BoundaryTerms&&) = default;
  ~BoundaryTerms() override = default;

  using variables_tags = tmpl::list<Var1, Var2<Dim>>;
  using variables_tag = Tags::Variables<variables_tags>;

  std::unique_ptr<BoundaryCorrection> get_clone() const override {
    return std::make_unique<BoundaryTerms>(*this);
  }

  void pup(PUP::er& p) override {  // NOLINT
    BoundaryCorrection::pup(p);
  }

  using dg_package_field_tags = tmpl::push_back<
       tmpl::append<db::wrap_tags_in<::Tags::NormalDotFlux, variables_tags>,
                    variables_tags>,
       MaxAbsCharSpeed>;
  using dg_package_data_temporary_tags = tmpl::list<>;
  using dg_package_data_primitive_tags =
      tmpl::conditional_t<HasPrims, tmpl::list<>, tmpl::list<>>;
  using dg_package_data_volume_tags =
      tmpl::conditional_t<HasPrims, tmpl::list<>, tmpl::list<>>;
  using dg_boundary_terms_volume_tags =
      tmpl::conditional_t<HasPrims, tmpl::list<>, tmpl::list<>>;

  double dg_package_data(
    const gsl::not_null<Scalar<DataVector>*> out_normal_dot_flux_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        out_normal_dot_flux_var2,
    const gsl::not_null<Scalar<DataVector>*> out_var1,
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> out_var2,
    const gsl::not_null<Scalar<DataVector>*> max_abs_char_speed,

    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,

    const tnsr::I<DataVector, Dim, Frame::Inertial>& flux_var1,
    const tnsr::IJ<DataVector, Dim, Frame::Inertial>& flux_var2,

    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
        mesh_velocity,
    const std::optional<Scalar<DataVector>>& normal_dot_mesh_velocity) const {
        return -999.0;
    }
};

template <size_t Dim>
class BoundaryCondition : public domain::BoundaryConditions::BoundaryCondition {
 public:
  BoundaryCondition() = default;
  BoundaryCondition(BoundaryCondition&&) = default;
  BoundaryCondition& operator=(BoundaryCondition&&) = default;
  BoundaryCondition(const BoundaryCondition&) = default;
  BoundaryCondition& operator=(const BoundaryCondition&) = default;
  ~BoundaryCondition() override = default;
  explicit BoundaryCondition(CkMigrateMessage* msg)
      : domain::BoundaryConditions::BoundaryCondition(msg) {}

  void pup(PUP::er& p) override {
    domain::BoundaryConditions::BoundaryCondition::pup(p);
  }
};

template <size_t Dim>
class CountWhenCalled : public BoundaryCondition<Dim> {
 public:
  CountWhenCalled() = default;
  CountWhenCalled(CountWhenCalled&&) = default;
  CountWhenCalled& operator=(CountWhenCalled&&) = default;
  CountWhenCalled(const CountWhenCalled&) = default;
  CountWhenCalled& operator=(const CountWhenCalled&) =
      default;
  ~CountWhenCalled() override = default;

  explicit CountWhenCalled(CkMigrateMessage* msg)
      : BoundaryCondition<Dim>(msg) {}

  WRAPPED_PUPable_decl_base_template(
      domain::BoundaryConditions::BoundaryCondition, CountWhenCalled);

  auto get_clone() const -> std::unique_ptr<
      domain::BoundaryConditions::BoundaryCondition> override {
    return std::make_unique<CountWhenCalled<Dim>>(*this);
  }

  // This determines the function called in BoundaryConditionsImpl.hpp
  // namely dg_demand_outgoing_char_speeds
  static constexpr ::evolution::BoundaryConditions::Type bc_type =
      ::evolution::BoundaryConditions::Type::DemandOutgoingCharSpeeds;

  void pup(PUP::er& p) override { BoundaryCondition<Dim>::pup(p); }

  using dg_interior_evolved_variables_tags = tmpl::list<>;
  using dg_interior_primitive_variables_tags = tmpl::list<>;
  using dg_interior_temporary_tags = tmpl::list<>;
  using dg_gridless_tags = tmpl::list<>;

  // Just count the number of times the boundary condition is called
  static std::optional<std::string> dg_demand_outgoing_char_speeds(
      const std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>>&
      /*face_mesh_velocity*/,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
      /*outward_directed_normal_covector*/) {
    CountWhenCalled<Dim>::number_of_times_called += 1;
    return std::nullopt;
  }

  static size_t number_of_times_called;
};

template <size_t Dim>
PUP::able::PUP_ID CountWhenCalled<Dim>::my_PUP_ID = 0;

template <size_t Dim>
size_t CountWhenCalled<Dim>::number_of_times_called = 0;

template <size_t Dim, SystemType system_type, bool HasPrimitiveVariables,
          bool PassVariables>
struct System {
  static constexpr bool has_primitive_and_conservative_vars =
      HasPrimitiveVariables;
  static constexpr size_t volume_dim = Dim;

  using boundary_conditions_base = BoundaryCondition<Dim>;

  using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim>>>;
  using flux_variables = tmpl::conditional_t<
      system_type == SystemType::Conservative, tmpl::list<Var1, Var2<Dim>>,
      tmpl::conditional_t<system_type == SystemType::Nonconservative,
                          tmpl::list<>, tmpl::list<Var2<Dim>>>>;
  using gradient_variables = tmpl::conditional_t<
      system_type == SystemType::Conservative, tmpl::list<>,
      tmpl::conditional_t<system_type == SystemType::Nonconservative,
                          tmpl::list<Var1, Var2<Dim>>, tmpl::list<Var1>>>;
  using compute_volume_time_derivative_terms = tmpl::conditional_t<
      PassVariables,
      TimeDerivativeTermsWithVariables<Dim, system_type,
                                       has_primitive_and_conservative_vars>,
      TimeDerivativeTerms<Dim, system_type,
                          has_primitive_and_conservative_vars>>;

};

template <typename System>
struct InitializeVars {
   template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static Parallel::iterable_action_return_t apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) {
    static constexpr size_t volume_dim = System::volume_dim;
    using variables_tag = System::variables_tag;
    const auto& x =
      db::get<domain::Tags::Coordinates<volume_dim, Frame::Inertial>>(box);
    db::mutate<variables_tag, Var3> ([&x](
      const gsl::not_null<typename variables_tag::type*> evolved_vars,
      const gsl::not_null<Scalar<DataVector>*> var3) {
        const DataVector r_squared = get(dot_product(x, x));
        Var1::value(make_not_null(&get<Var1>(*evolved_vars)), r_squared);
        Var2<volume_dim>::value(
          make_not_null(&get<Var2<volume_dim>>(*evolved_vars)), x);
        Var3::value(var3, r_squared);
      }, make_not_null(&box));
    return {Parallel::AlgorithmExecution::Continue, std::nullopt};
  }
};

template <typename Metavariables>
struct component {
  using metavariables = Metavariables;
  static constexpr size_t volume_dim = metavariables::volume_dim;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = ElementId<volume_dim>;
  using time_stepper_base = metavariables::TimeStepperBase;
  using system = metavariables::system;

  using variables_tag = typename Metavariables::system::variables_tag;
  using simple_tags = tmpl::list<
    domain::Tags::InitialExtents<volume_dim>,
    domain::Tags::InitialRefinementLevels<volume_dim>,
    ::evolution::dg::Tags::Quadrature,
    ::Tags::Time,
    Initialization::Tags::InitialTimeDelta,
    Initialization::Tags::InitialSlabSize<metavariables::local_time_stepping>,
    ::Tags::StepperErrorEstimatesEnabled,
    Var3
    >;
  using compute_tags = tmpl::list<>;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          Parallel::Phase::Initialization,
          tmpl::flatten<tmpl::list<
             ActionTesting::InitializeDataBox<simple_tags, compute_tags>,
             Initialization::Actions::InitializeItems<
             Initialization::TimeStepping<metavariables, time_stepper_base>,
             ::evolution::dg::Initialization::Domain<metavariables>,
             Initialization::TimeStepperHistory<metavariables>>,
             Initialization::Actions::ConservativeSystem<system>,
             InitializeVars<system>,
             ::evolution::dg::Initialization::Mortars<volume_dim, system>,
             ::Actions::MutateApply<AdvanceTime>>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<
              ::evolution::dg::Actions::ComputeTimeDerivative<
              volume_dim, typename Metavariables::system,
              AllStepChoosers, Metavariables::local_time_stepping,
              Metavariables::use_nodegroup_dg_elements>>>>;
};

template <size_t Dim, SystemType SystemTypeIn, bool LocalTimeStepping,
          bool UseMovingMesh, bool HasPrimitiveVariables, bool PassVariables,
          bool UseNodegroupDgElements>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;
  static constexpr SystemType system_type = SystemTypeIn;
  static constexpr bool use_moving_mesh = UseMovingMesh;
  static constexpr bool local_time_stepping = LocalTimeStepping;
  static constexpr bool pass_variables = PassVariables;
  static constexpr bool use_nodegroup_dg_elements = UseNodegroupDgElements;
  using system =
      System<Dim, system_type, HasPrimitiveVariables, pass_variables>;
  using TimeStepperBase =
      tmpl::conditional_t<LocalTimeStepping, LtsTimeStepper, TimeStepper>;

  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes = tmpl::map<
        tmpl::pair<BoundaryCondition<Dim>,
                   tmpl::list<CountWhenCalled<Dim>>>,
        tmpl::pair<::evolution::BoundaryCorrection,
                   tmpl::list<BoundaryTerms<Dim, HasPrimitiveVariables>>>,
        tmpl::pair<StepChooser<StepChooserUse::LtsStep>,
                   tmpl::list<StepChoosers::Constant>>>;
  };
  using component_list = tmpl::list<component<Metavariables>>;
  using const_global_cache_tags =
      tmpl::list<domain::Tags::Domain<Dim>,
                 ::Tags::ConcreteTimeStepper<TimeStepperBase>>;
};

template <size_t Dim, bool UseMovingMesh>
std::unique_ptr<DomainCreator<Dim>> make_domain_creator(const bool conforming)
{
  if (conforming) {
    if constexpr (Dim == 1) {
      return std::make_unique<domain::creators::RotatedIntervals>(
        std::array{-3.0}, std::array{1.5}, std::array{2.0}, std::array{0_st},
        std::array<std::array<size_t, 2>, 1>{std::array{6_st, 6_st}},
        std::make_unique<CountWhenCalled<Dim>>(),
        std::make_unique<CountWhenCalled<Dim>>(),
        UseMovingMesh ? nullptr : nullptr);
    }
    if constexpr (Dim == 2) {
      return std::make_unique<domain::creators::RotatedRectangles>(
        std::array{-3.0, -2.0}, std::array{1.5, 4.5}, std::array{2.0, 7.0},
        std::array{0_st, 0_st},
        std::array<std::array<size_t, 2>, 2>{std::array{6_st, 6_st},
        std::array{6_st, 6_st}},
        std::make_unique<CountWhenCalled<Dim>>());
    }
    if constexpr (Dim == 3) {
      return std::make_unique<domain::creators::RotatedBricks>(
        std::array{-3.0, -2.0, -1.0}, std::array{1.5, 4.5, 2.5},
        std::array{2.0, 7.0, 3.0}, std::array{0_st, 0_st, 0_st},
        std::array<std::array<size_t, 2>, 3>{std::array{6_st, 6_st},
        std::array{6_st, 6_st}, std::array{6_st, 6_st}},
        std::make_unique<CountWhenCalled<Dim>>());
    }
  } else {
    if constexpr (Dim == 3) {
        return std::make_unique<domain::creators::NonconformingSphericalShells>(
          2.0, 3.0, 4.0, 0, 0, 5, 5, 9,
          std::make_unique<CountWhenCalled<Dim>>(),
          std::make_unique<CountWhenCalled<Dim>>());
    }
  }
}

template <bool LocalTimeStepping, bool UseMovingMesh, size_t Dim,
          SystemType system_type, bool HasPrims, bool PassVariables,
          bool UseNodegroupDgElements>
void test_impl(const Spectral::Quadrature quadrature,
               const ::dg::Formulation dg_formulation,
               const bool domain_is_conforming) {
  CAPTURE(LocalTimeStepping);
  CAPTURE(UseMovingMesh);
  CAPTURE(Dim);
  CAPTURE(system_type);
  CAPTURE(HasPrims);
  CAPTURE(PassVariables);
  CAPTURE(quadrature);
  CAPTURE(dg_formulation);
  using metavars =
      Metavariables<Dim, system_type, LocalTimeStepping, UseMovingMesh,
                    HasPrims, PassVariables, UseNodegroupDgElements>;
  register_classes_with_charm<TimeSteppers::AdamsBashforth>();
  register_factory_classes_with_charm<metavars>();
  domain::creators::register_derived_with_charm();

  using system = typename metavars::system;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavars>;
  using variables_tag = typename system::variables_tag;
  using flux_variables = typename system::flux_variables;
  using flux_variables_tag = ::Tags::Variables<flux_variables>;
  using fluxes_tag = db::add_tag_prefix<::Tags::Flux, flux_variables_tag,
                                        tmpl::size_t<Dim>, Frame::Inertial>;
  using dt_variables_tag = db::add_tag_prefix<::Tags::dt, variables_tag>;

  const auto creator =
    make_domain_creator<Dim, UseMovingMesh>(domain_is_conforming);
  auto domain = creator->create_domain();
  const size_t num_blocks = domain.blocks().size();
  const auto initial_refinement = creator->initial_refinement_levels();
  const auto initial_extents = creator->initial_extents();
  const double initial_time = 1.5;
  const double initial_dt = 0.5;
  const double initial_slab_size = LocalTimeStepping ? 4.5 : initial_dt;
  std::unique_ptr<typename metavars::TimeStepperBase> time_stepper =
      std::make_unique<TimeSteppers::AdamsBashforth>(5);
  std::unordered_map<std::string,
                     std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>
      functions_of_time{};
  tuples::TaggedTuple<domain::Tags::FunctionsOfTimeInitialize>
     mutable_global_cache_items(std::move(functions_of_time));

  ActionTesting::MockRuntimeSystem<metavars> runner =
    [&initial_dt, &domain, &dg_formulation, &time_stepper, &creator,
     &mutable_global_cache_items]() {
    if constexpr(LocalTimeStepping) {
      std::vector<std::unique_ptr<StepChooser<StepChooserUse::LtsStep>>>
        step_choosers;
      step_choosers.emplace_back(std::make_unique<StepChoosers::Constant>(
        initial_dt));
      tuples::TaggedTuple<
        domain::Tags::Domain<Dim>,
        ::Tags::ConcreteTimeStepper<typename metavars::TimeStepperBase>,
        ::dg::Tags::Formulation,
        ::evolution::Tags::BoundaryCorrection,
        domain::Tags::ExternalBoundaryConditions<Dim>,
        Tags::StepChoosers,
        Tags::MinimumTimeStep>
          const_global_cache_items(std::move(domain), std::move(time_stepper),
            dg_formulation, std::make_unique<BoundaryTerms<Dim, HasPrims>>(),
            creator->external_boundary_conditions(), std::move(step_choosers),
            1.0e-8);
      return MockRuntimeSystem{std::move(const_global_cache_items),
                               std::move(mutable_global_cache_items)};
    } else {
      tuples::TaggedTuple<
        domain::Tags::Domain<Dim>,
        ::Tags::ConcreteTimeStepper<typename metavars::TimeStepperBase>,
        ::dg::Tags::Formulation,
        ::evolution::Tags::BoundaryCorrection,
        domain::Tags::ExternalBoundaryConditions<Dim>>
          const_global_cache_items(std::move(domain), std::move(time_stepper),
            dg_formulation, std::make_unique<BoundaryTerms<Dim, HasPrims>>(),
            creator->external_boundary_conditions());
      return MockRuntimeSystem{std::move(const_global_cache_items),
                               std::move(mutable_global_cache_items)};
    }
  }();

  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::emplace_component_and_initialize<component<metavars>>(
      &runner,
      element_id,
      {initial_extents,
       initial_refinement,
       quadrature,
       initial_time,
       initial_dt,
       initial_slab_size,
       false,
       Scalar<DataVector>{}
       });
  }
  for (size_t i = 0; i < 5; ++i) {
    for (size_t b = 0; b < num_blocks; ++b) {
      ElementId<Dim> element_id{b};
        ActionTesting::next_action<component<metavars>>(
           make_not_null(&runner), element_id);
    }
  }
  // Start testing the actual dg::ComputeTimeDerivative action
  const ElementId<Dim> id{0};
  const auto& first_box = get_databox<component<metavars>>(runner, id);
  Parallel::printf("%s\n", first_box.print_tags());
  CountWhenCalled<Dim>::number_of_times_called = 0;
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);
  using tags_for_unchanged_items = tmpl::list<
    Tags::Time,
    Tags::TimeStepId,
    Tags::StepNumberWithinSlab,
    Tags::AdaptiveSteppingDiagnostics,
    Tags::Next<Tags::TimeStepId>,
    Tags::TimeStep,
    Tags::ChangeSlabSize::SlabSizeGoal,
    domain::Tags::Mesh<Dim>,
    domain::Tags::Element<Dim>,
    domain::Tags::NeighborMesh<Dim>,
    variables_tag,
    Tags::HistoryEvolvedVariables<variables_tag>,
    ::evolution::dg::Tags::MortarMesh<Dim>,
    ::evolution::dg::Tags::MortarInfo<Dim>,
    ::evolution::dg::Tags::MortarNextTemporalId<Dim>
    >;
  const auto items_before = copy_items<tags_for_unchanged_items>(first_box);
  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    ActionTesting::next_action<component<metavars>>(make_not_null(&runner),
                                                    element_id);
  }
  const auto items_after = copy_items<tags_for_unchanged_items>(first_box);
  CHECK(items_before == items_after);
  CHECK(CountWhenCalled<Dim>::number_of_times_called == num_blocks);
  // Check dt_evolved_vars
  for (size_t b = 0; b < num_blocks; ++b) {
    ElementId<Dim> element_id{b};
    const auto& box = get_databox<component<metavars>>(runner, element_id);
    const auto& x =
      db::get<domain::Tags::Coordinates<Dim, Frame::Inertial>>(box);
    const DataVector r_squared = get(dot_product(x, x));
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    Variables<tmpl::list<::Tags::dt<Var1>, ::Tags::dt<Var2<Dim>>>>
      expected_dt_evolved_vars{mesh.number_of_grid_points()};
    Var1::dt_value(make_not_null(
                     &get<Tags::dt<Var1>>(expected_dt_evolved_vars)),
                     r_squared, x);
    const auto& dt_evolved_vars = ActionTesting::get_databox_tag<
                  component<metavars>, dt_variables_tag>(
                  runner, element_id);
    CHECK_ITERABLE_APPROX(get<Tags::dt<Var1>>(dt_evolved_vars),
                          get<Tags::dt<Var1>>(expected_dt_evolved_vars));
    // CHECK_VARIABLES_APPROX(
    //   SINGLE_ARG(ActionTesting::get_databox_tag<
    //              component<metavars>, dt_variables_tag>(
    //                runner, element_id)),
    //   expected_dt_evolved_vars);
  }
//      dt_variables_tag,
//      ::evolution::dg::Tags::MortarData<Dim>,
//      ::evolution::dg::Tags::NormalCovectorAndMagnitude<Dim>,
//      ::evolution::dg::Tags::MortarDataHistory<Dim,
//                                           typename dt_variables_tag::type>
    // Check NormalCovectorAndMagnitude
    // Check MortarData
    // Check TimeStep (if LTS)
    // Check MortarDataHistory (if LTS)
    // Check Inboxes
}

template <SystemType system_type, bool UsePrims, size_t Dim>
void test_nonconforming() {
  constexpr bool use_nodegroup_dg_elements = false;

  const auto invoke_tests_with_quadrature_and_formulation =
      [](const Spectral::Quadrature quadrature,
         const ::dg::Formulation local_dg_formulation) {
        const auto moving_mesh_helper = [&local_dg_formulation,
                                         &quadrature](auto moving_mesh) {
          // PassVariables == false
          test_impl<false, std::decay_t<decltype(moving_mesh)>::value, Dim,
                    system_type, UsePrims, false, use_nodegroup_dg_elements>(
              quadrature, local_dg_formulation, true);
          // test_impl<true, std::decay_t<decltype(moving_mesh)>::value, Dim,
          //           system_type, UsePrims, false, use_nodegroup_dg_elements>(
          //     quadrature, local_dg_formulation);

          // PassVariables == true
          // test_impl<false, std::decay_t<decltype(moving_mesh)>::value, Dim,
          //           system_type, UsePrims, true, use_nodegroup_dg_elements>(
          //     quadrature, local_dg_formulation);
          // test_impl<true, std::decay_t<decltype(moving_mesh)>::value, Dim,
          //           system_type, UsePrims, true, use_nodegroup_dg_elements>(
          //     quadrature, local_dg_formulation);
        };
        moving_mesh_helper(std::integral_constant<bool, false>{});
        // moving_mesh_helper(std::integral_constant<bool, true>{});
      };
  invoke_tests_with_quadrature_and_formulation(
        Spectral::Quadrature::GaussLobatto, ::dg::Formulation::StrongInertial);

}
}  // namespace TestHelpers::evolution::dg::Actions

// Distributed under the MIT License.
// See LICENSE.txt for details.
/// \cond
#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "Domain/Amr/Actions/AdjustDomain.hpp"
#include "Domain/Amr/Actions/EvaluateRefinementCriteria.hpp"
#include "Domain/Amr/Actions/InitializeAmr.hpp"
#include "Domain/Amr/Actions/InitializeDomain.hpp"
#include "Domain/Amr/Actions/SendGlobalAmrDiagnostics.hpp"
#include "Domain/Amr/Components/AmrMonitor.hpp"
#include "Domain/Amr/Components/Element.hpp"
#include "Domain/Creators/DomainCreator.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Options/Options.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct DummySystem {};
}  // namespace

/// \page ExploreAmrExecutablePage ExploreAmr Executable
/// The ExploreAmr executable is being used to develop the mechanics of
/// adaptive mesh refinement.
///
/// See ExploreAmrMetavars for a description of the metavariables of this
/// executable.

/// \brief The metavariables for the ExploreAmr executable
template <size_t Dim>
struct ExploreAmrMetavars {
  static constexpr size_t volume_dim = Dim;
  using system = DummySystem;

  static constexpr Options::String help{
      "Test anisotropic refinement by randomly refining a grid.\n"};

  /// The phases of the executable
  enum class Phase {
    Initialization,      /**< Initial phase used to initialize components */
    CheckDomain,         /**< Verify that the Domain is valid */
    EvaluateAmrCriteria, /**< Decide where to do adaptive mesh refinment */
    AdjustDomain,        /**< Adjust the Domain */
    Exit                 /**< Final phase of a successful execution */
  };

  /// \brief Items in the Parallel::ConstGlobalCache
  ///
  /// The Parallel::ConstGlobalCache contains the items in this list plus the
  /// items in the `const_global_cache_tag_list` of each component in
  /// ExploreAmrMetavars::component_list.  The items in the
  /// Parallel::ConstGlobalCache are determined from the input file.
  using const_global_cache_tags = tmpl::list<amr::Tags::NumberOfIterations>;

  /// The parallel components used in the executable
  using component_list = tmpl::list<
      Parallel::Components::AmrMonitor<ExploreAmrMetavars>,
      Parallel::Components::Element<
          ExploreAmrMetavars,
          tmpl::list<
              Parallel::PhaseActions<
                  Phase, Phase::Initialization,
                  tmpl::list<Actions::SetupDataBox,
                             ::amr::Actions::InitializeDomain<volume_dim>,
                             ::amr::Actions::InitializeAmr<volume_dim>,
                             ::Initialization::Actions::
                                 RemoveOptionsAndTerminatePhase>>,
              Parallel::PhaseActions<
                  Phase, Phase::CheckDomain,
                  tmpl::list<::amr::Actions::SendGlobalAmrDiagnostics,
                             Parallel::Actions::TerminatePhase>>,
              Parallel::PhaseActions<
                  Phase, Phase::EvaluateAmrCriteria,
                  tmpl::list<::amr::Actions::EvaluateRefinementCriteria,
                             Parallel::Actions::TerminatePhase>>,
              Parallel::PhaseActions<
                  Phase, Phase::AdjustDomain,
                  tmpl::list<::amr::Actions::AdjustDomain,
                             Parallel::Actions::TerminatePhase>>>>>;

  /// \brief Determine the next phase of the executable
  ///
  /// The phases are executed in the following order:
  /// - Phase::Initialization
  /// - Phase::Exit
  static Phase determine_next_phase(
      const Phase& current_phase,
      const Parallel::CProxy_GlobalCache<
          ExploreAmrMetavars>& /*cache_proxy*/) noexcept;
};

template <size_t Dim>
std::string name(
    const typename ExploreAmrMetavars<Dim>::Phase& phase) noexcept {
  switch (phase) {
    case ExploreAmrMetavars<Dim>::Phase::Initialization:
      return "Initialization";
    case ExploreAmrMetavars<Dim>::Phase::CheckDomain:
      return "CheckDomain";
    case ExploreAmrMetavars<Dim>::Phase::EvaluateAmrCriteria:
      return "EvaluateAmrCriteria";
    case ExploreAmrMetavars<Dim>::Phase::AdjustDomain:
      return "AdjustDomain";
    case ExploreAmrMetavars<Dim>::Phase::Exit:
      return "Exit";
  }
  ERROR("Switch failed. Did you static_cast<Phase> an integral value?");
}

template <size_t Dim>
typename ExploreAmrMetavars<Dim>::Phase
ExploreAmrMetavars<Dim>::determine_next_phase(
    const ExploreAmrMetavars::Phase& current_phase,
    const Parallel::CProxy_GlobalCache<ExploreAmrMetavars>&
        cache_proxy) noexcept {
  auto& local_cache = *(cache_proxy.ckLocalBranch());

  static size_t iterations = 0;
  const auto& number_of_iterations =
      get<amr::Tags::NumberOfIterations>(local_cache);
  switch (current_phase) {
    case ExploreAmrMetavars::Phase::Initialization:
      Parallel::printf("\nEntering phase: CheckDomain\n\n");
      return Phase::CheckDomain;
    case ExploreAmrMetavars::Phase::CheckDomain:
      if (iterations++ < number_of_iterations) {
        Parallel::printf("\nEntering phase: EvaluateAmrCriteria\n\n");
        return Phase::EvaluateAmrCriteria;
      }
      Parallel::printf("\nEntering phase:  Exit\n\n");
      return Phase::Exit;
    case ExploreAmrMetavars::Phase::EvaluateAmrCriteria:
      Parallel::printf("\nEntering phase: AdjustDomain\n\n");
      return Phase::AdjustDomain;
    case ExploreAmrMetavars::Phase::AdjustDomain:
      Parallel::printf(
          "\nEntering phase:  CheckDomain for iteration "
          "%zu\n\n--------------------------------------\n\n",
          iterations);
      return Phase::CheckDomain;
    case ExploreAmrMetavars::Phase::Exit:
      ERROR(
          "Should never call determine_next_phase with the current phase "
          "being 'Exit'");
  }
  ERROR("Switch failed. Did you static_cast<Phase> an integral value?");
}

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling, &domain::creators::register_derived_with_charm};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};

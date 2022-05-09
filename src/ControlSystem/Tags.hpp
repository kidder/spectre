// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <utility>

#include "ControlSystem/Averager.hpp"
#include "ControlSystem/Controller.hpp"
#include "ControlSystem/Protocols/ControlSystem.hpp"
#include "ControlSystem/TimescaleTuner.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Options/Options.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace ah {
enum class ObjectLabel;
}  // namespace ah
namespace control_system {
template <typename ControlSystem>
struct OptionHolder;
}  // namespace control_system
namespace OptionTags {
struct InitialTime;
}  // namespace OptionTags
/// \endcond

namespace control_system {
/// \ingroup ControlSystemGroup
/// All tags that will be used in the LinkedMessageQueue's within control
/// systems.
///
/// These tags will be used to retreive the results of the measurements that
/// were sent to the control system which have been placed inside a
/// LinkedMessageQueue.
namespace QueueTags {
/// \ingroup ControlSystemGroup
/// Holds the centers of each horizon from measurements as DataVectors
template <::ah::ObjectLabel Horizon>
struct Center {
  using type = DataVector;
};
}  // namespace QueueTags

/// \ingroup ControlSystemGroup
/// All option tags related to the control system
namespace OptionTags {
/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Options group for all control system options
struct ControlSystemGroup {
  static std::string name() { return "ControlSystems"; }
  static constexpr Options::String help = {
      "Options for all control systems used in a simulation."};
};

/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Option tag for each individual control system. The name of this option is
/// the name of the \p ControlSystem struct it is templated on. This way all
/// control systems will have a unique name.
template <typename ControlSystem>
struct ControlSystemInputs {
  using type = control_system::OptionHolder<ControlSystem>;
  static constexpr Options::String help{"Options for a control system."};
  static std::string name() { return ControlSystem::name(); }
  using group = ControlSystemGroup;
};

/// \ingroup OptionTagsGroup
/// \ingroup ControlSystemGroup
/// Option tag on whether to write data to disk.
struct WriteDataToDisk {
  using type = bool;
  static constexpr Options::String help = {
      "Whether control system data should be saved during an evolution."};
  using group = ControlSystemGroup;
};
}  // namespace OptionTags

/// \ingroup ControlSystemGroup
/// Alias to get all the option holders from a list of control systems. This is
/// useful in the `option_tags` alias of simple tags for getting all the options
/// from control systems.
template <typename ControlSystems>
using inputs =
    tmpl::transform<ControlSystems,
                    tmpl::bind<OptionTags::ControlSystemInputs, tmpl::_1>>;

/// \ingroup ControlSystemGroup
/// All DataBox tags related to the control system
namespace Tags {
/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for writing control system data to disk
struct WriteDataToDisk : db::SimpleTag {
  using type = bool;
  using option_tags = tmpl::list<OptionTags::WriteDataToDisk>;

  static constexpr bool pass_metavariables = false;
  static type create_from_options(const type& option) { return option; }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the averager
///
/// To compute the `deriv_order`th derivative of a control error, the max
/// derivative we need from the averager is the `deriv_order - 1`st derivative.
template <typename ControlSystem>
struct Averager : db::SimpleTag {
  using type = ::Averager<ControlSystem::deriv_order - 1>;

  using option_tags =
      tmpl::list<OptionTags::ControlSystemInputs<ControlSystem>>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(
      const control_system::OptionHolder<ControlSystem>& option_holder) {
    return option_holder.averager;
  }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the timescale tuner
template <typename ControlSystem>
struct TimescaleTuner : db::SimpleTag {
  using type = ::TimescaleTuner;

  using option_tags =
      tmpl::list<OptionTags::ControlSystemInputs<ControlSystem>>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(
      const control_system::OptionHolder<ControlSystem>& option_holder) {
    return option_holder.tuner;
  }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the controller
template <typename ControlSystem>
struct Controller : db::SimpleTag {
  using type = ::Controller<ControlSystem::deriv_order>;

  using option_tags =
      tmpl::list<::OptionTags::InitialTime,
                 OptionTags::ControlSystemInputs<ControlSystem>>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(
      const double initial_time,
      const control_system::OptionHolder<ControlSystem>& option_holder) {
    type controller = option_holder.controller;
    const ::TimescaleTuner tuner = option_holder.tuner;

    controller.set_initial_update_time(initial_time);
    controller.assign_time_between_updates(min(tuner.current_timescale()));

    return controller;
  }
};

/// \ingroup DataBoxTagsGroup
/// \ingroup ControlSystemGroup
/// DataBox tag for the control error
template <typename ControlSystem>
struct ControlError : db::SimpleTag {
  using type = typename ControlSystem::control_error;

  using option_tags =
      tmpl::list<OptionTags::ControlSystemInputs<ControlSystem>>;
  static constexpr bool pass_metavariables = false;
  static type create_from_options(
      const control_system::OptionHolder<ControlSystem>& option_holder) {
    return option_holder.control_error;
  }
};
}  // namespace Tags

/// \ingroup ControlSystemGroup
/// Holds all options for a single control system
///
/// This struct collects all the options for a given control system during
/// option parsing. Then during initialization, the options can be retrieved via
/// their public member names and assigned to their corresponding DataBox tags.
template <typename ControlSystem>
struct OptionHolder {
  static_assert(tt::assert_conforms_to<
                ControlSystem, control_system::protocols::ControlSystem>);
  using control_system = ControlSystem;
  static constexpr size_t deriv_order = control_system::deriv_order;
  struct Averager {
    using type = ::Averager<deriv_order - 1>;
    static constexpr Options::String help = {
        "Averages the derivatives of the control error and possibly the "
        "control error itself."};
  };

  struct Controller {
    using type = ::Controller<deriv_order>;
    static constexpr Options::String help = {
        "Computes the control signal which will be used to reset the functions "
        "of time."};
  };

  struct TimescaleTuner {
    using type = ::TimescaleTuner;
    static constexpr Options::String help = {
        "Keeps track of the damping timescales for the control system upon "
        "which other timescales are based of off."};
  };

  struct ControlError {
    using type = typename ControlSystem::control_error;
    static constexpr Options::String help = {
        "Computes the control error for the control system based on quantities "
        "in the simulation."};
  };

  using options =
      tmpl::list<Averager, Controller, TimescaleTuner, ControlError>;
  static constexpr Options::String help = {"Options for a control system."};

  OptionHolder(::Averager<deriv_order - 1> input_averager,
               ::Controller<deriv_order> input_controller,
               ::TimescaleTuner input_tuner,
               typename ControlSystem::control_error input_control_error)
      : averager(std::move(input_averager)),
        controller(std::move(input_controller)),
        tuner(std::move(input_tuner)),
        control_error(std::move(input_control_error)) {}

  OptionHolder() = default;
  OptionHolder(const OptionHolder& /*rhs*/) = default;
  OptionHolder& operator=(const OptionHolder& /*rhs*/) = default;
  OptionHolder(OptionHolder&& /*rhs*/) = default;
  OptionHolder& operator=(OptionHolder&& /*rhs*/) = default;
  ~OptionHolder() = default;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) {
    p | averager;
    p | controller;
    p | tuner;
    p | control_error;
  };

  // These members are specifically made public for easy access during
  // initialization
  ::Averager<deriv_order - 1> averager{};
  ::Controller<deriv_order> controller{};
  ::TimescaleTuner tuner{};
  typename ControlSystem::control_error control_error{};
};
}  // namespace control_system

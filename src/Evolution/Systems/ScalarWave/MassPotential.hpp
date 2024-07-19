// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
namespace gsl {
template <typename T>
class not_null;
}  // namespace gsl
/// \endcond

namespace ScalarWave {
/// @{
/*!
 * \brief Computes the mass potential of the scalar wave system.
 *
 * Below is the function used to calculate the mass potential.
 *
 * \f{align*}
 * U = m^2 * \psi
 * \f}
 */
void mass_potential(gsl::not_null<Scalar<DataVector>*> result,
                    const Scalar<DataVector>& psi, const double& mass);

Scalar<DataVector> mass_potential(const Scalar<DataVector>& psi,
                                  const double& mass);
/// @}

namespace Tags {
/// \brief Computes the mass potential using ScalarWave::mass_potential()
struct MassPotentialCompute : MassPotential, db::ComputeTag {
  using argument_tags = tmpl::list<Psi, MassValue>;

  using return_type = Scalar<DataVector>;

  static constexpr auto function =
      static_cast<void (*)(gsl::not_null<Scalar<DataVector>*> result,
                           const Scalar<DataVector>&, const double&)>(
          &mass_potential);

  using base = MassPotential;
};
}  // namespace Tags
}  // namespace ScalarWave

// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/ScalarWave/MassPotential.hpp"

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

namespace ScalarWave {
void mass_potential(gsl::not_null<Scalar<DataVector>*> result,
                    const Scalar<DataVector>& psi, const double& mass) {
  get(*result) = square(mass) * get(psi);
}

Scalar<DataVector> mass_potential(const Scalar<DataVector>& psi,
                                  const double& mass) {
  Scalar<DataVector> result{get(psi).size()};
  mass_potential(make_not_null(&result), psi, mass);
  return result;
}
.
}  // namespace ScalarWave


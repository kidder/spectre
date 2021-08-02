// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Evolution/Systems/GrMhd/ValenciaDivClean/Subcell/TciOnFdGrid.hpp"

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Evolution/DgSubcell/PerssonTci.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"

namespace grmhd::ValenciaDivClean::subcell {
bool TciOnFdGrid::apply(const Scalar<DataVector>& tilde_d,
                        const Scalar<DataVector>& tilde_tau,
                        const bool vars_needed_fixing, const Mesh<3>& dg_mesh,
                        const TciOptions& tci_options,
                        const double persson_exponent) noexcept {
  return vars_needed_fixing or
         min(get(tilde_d)) <
             tci_options.minimum_rest_mass_density_times_lorentz_factor or
         min(get(tilde_tau)) < tci_options.minimum_tilde_tau or
         evolution::dg::subcell::persson_tci(tilde_d, dg_mesh, persson_exponent,
                                             1.0e-18) or
         evolution::dg::subcell::persson_tci(tilde_tau, dg_mesh,
                                             persson_exponent, 1.0e-18);
}
}  // namespace grmhd::ValenciaDivClean::subcell

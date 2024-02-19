// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Index.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/Gsl.hpp"

namespace evolution::dg::subcell {
namespace detail {
template <size_t Dim>
bool persson_tci_impl(gsl::not_null<DataVector*> filtered_component,
                      const DataVector& component, const Mesh<Dim>& dg_mesh,
                      double alpha, size_t num_highest_modes);
}  // namespace detail

/*!
 * \brief Troubled cell indicator using the idea of spectral falloff by
 *    \cite Persson2006sub
 *
 * Consider a discontinuity sensing quantity \f$U\f$, which is typically a
 * scalar but could be a tensor of any rank. Let \f$U\f$ have the 1d spectral
 * decomposition (generalization to higher-dimensional tensor product bases is
 * done dimension-by-dimension):
 *
 * \f{align*}{
 *   U(x)=\sum_{i=0}^{N}c_i P_i(x),
 * \f}
 *
 * where \f$P_i(x)\f$ are the basis functions, in our case the Legendre
 * polynomials, and \f$c_i\f$ are the spectral coefficients. We then define a
 * filtered solution \f$\hat{U}\f$ as
 *
 * \f{align*}{
 *   \hat{U}(x)=\sum_{i=N+1-M}^{N} c_i P_i(x).
 * \f}
 *
 * where $M$ is the number of highest modes to include in the filtered solution.
 *
 * Note that when an exponential filter is being used to deal with aliasing,
 * lower modes can be included in \f$\hat{U}\f$. The main goal of \f$\hat{U}\f$
 * is to measure how much power is in the highest modes, which are the modes
 * responsible for Gibbs phenomena.
 *
 * A cell is troubled if
 *
 * \f{align*}{
 *    \frac{(\hat{U}, \hat{U})}{(U, U)} > (N + 1 - M)^{-\alpha}
 * \f}
 *
 * where \f$(\cdot,\cdot)\f$ is an inner product, which we take to be the
 * Euclidean \f$L_2\f$ norm (i.e. we do not divide by the number of grid points
 * since that cancels out anyway) computed as
 *
 * \f{align*}{
 *    (U, U) \equiv \sqrt{\sum_{i=0}^N U_i^2}
 * \f}
 *
 * where $U_i$ are nodal values of the quantity $U$ at grid points.
 *
 * Typically, \f$\alpha=4.0\f$ and $M=1$ is a good choice.
 *
 */
template <size_t Dim, typename SymmList, typename IndexList>
bool persson_tci(const Tensor<DataVector, SymmList, IndexList>& tensor,
                 const Mesh<Dim>& dg_mesh, const double alpha,
                 const size_t num_highest_modes) {
  DataVector filtered_component(dg_mesh.number_of_grid_points());
  for (size_t component_index = 0; component_index < tensor.size();
       ++component_index) {
    if (detail::persson_tci_impl(make_not_null(&filtered_component),
                                 tensor[component_index], dg_mesh, alpha,
                                 num_highest_modes)) {
      return true;
    }
  }
  return false;
}
}  // namespace evolution::dg::subcell

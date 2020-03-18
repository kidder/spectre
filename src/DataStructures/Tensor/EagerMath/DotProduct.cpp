// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"


#include <cstddef>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/MakeWithValue.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Gsl.hpp"

/// \cond
template <typename DataType, typename Index>
void dot_product(
    const gsl::not_null<Scalar<DataType>*> dot_product,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_b) noexcept {
  get(*dot_product) = get<0>(vector_a) * get<0>(vector_b);
  for (size_t d = 1; d < Index::dim; ++d) {
    get(*dot_product) += vector_a.get(d) * vector_b.get(d);
  }
}

template <typename DataType, typename Index>
Scalar<DataType> dot_product(
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_b) noexcept {
  Scalar<DataType> dot_product(vector_a.size());
  ::dot_product(make_not_null(&dot_product), vector_a, vector_b);
  return dot_product;
}

template <typename DataType, typename Index>
void dot_product(
    const gsl::not_null<Scalar<DataType>*> dot_product,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<change_index_up_lo<Index>>>&
        vector_b) noexcept {
  get(*dot_product) = get<0>(vector_a) * get<0>(vector_b);
  for (size_t d = 1; d < Index::dim; ++d) {
    get(*dot_product) += vector_a.get(d) * vector_b.get(d);
  }
}

template <typename DataType, typename Index>
Scalar<DataType> dot_product(
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<change_index_up_lo<Index>>>&
        vector_b) noexcept {
  Scalar<DataType> dot_product(vector_a.size());
  ::dot_product(make_not_null(&dot_product), vector_a, vector_b);
  return dot_product;
}
template <typename DataType, typename Index>
void dot_product(
    const gsl::not_null<Scalar<DataType>*> dot_product,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_b,
    const Tensor<DataType, Symmetry<1, 1>,
                 index_list<change_index_up_lo<Index>,
                            change_index_up_lo<Index>>>& metric) noexcept {
  get(*dot_product) = get<0>(vector_a) * get<0>(vector_b) * get<0, 0>(metric);
  for (size_t b = 1; b < Index::dim; ++b) {
    get(*dot_product) += get<0>(vector_a) * vector_b.get(b) * metric.get(0, b);
  }

  for (size_t a = 1; a < Index::dim; ++a) {
    for (size_t b = 0; b < Index::dim; ++b) {
      get(*dot_product) += vector_a.get(a) * vector_b.get(b) * metric.get(a, b);
    }
  }
}

template <typename DataType, typename Index>
Scalar<DataType> dot_product(
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_a,
    const Tensor<DataType, Symmetry<1>, index_list<Index>>& vector_b,
    const Tensor<DataType, Symmetry<1, 1>,
                 index_list<change_index_up_lo<Index>,
                            change_index_up_lo<Index>>>& metric) noexcept {
  Scalar<DataType> dot_product(vector_a.size());
  ::dot_product(make_not_null(&dot_product), vector_a, vector_b, metric);
  return dot_product;
}

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)
#define DTYPE(data) BOOST_PP_TUPLE_ELEM(1, data)
#define FRAME(data) BOOST_PP_TUPLE_ELEM(2, data)
#define INDEXTYPE(data) BOOST_PP_TUPLE_ELEM(3, data)
#define UPORLO(data) BOOST_PP_TUPLE_ELEM(4, data)
#define INDEX(data) INDEXTYPE(data)<DIM(data), UPORLO(data), FRAME(data)>

#define INSTANTIATE(_, data)                                           \
  template Scalar<DTYPE(data)> dot_product(                            \
      const Tensor<DTYPE(data), Symmetry<1>, index_list<INDEX(data)>>& \
          vector_a,                                                    \
      const Tensor<DTYPE(data), Symmetry<1>, index_list<INDEX(data)>>& \
          vector_b) noexcept;                                          \
  template Scalar<DTYPE(data)> dot_product(                            \
      const Tensor<DTYPE(data), Symmetry<1>, index_list<INDEX(data)>>& \
          vector_a,                                                    \
      const Tensor<DTYPE(data), Symmetry<1>,                           \
                   index_list<change_index_up_lo<INDEX(data)>>>&       \
          vector_b) noexcept;                                          \
  template Scalar<DTYPE(data)> dot_product(                            \
      const Tensor<DTYPE(data), Symmetry<1>, index_list<INDEX(data)>>& \
          vector_a,                                                    \
      const Tensor<DTYPE(data), Symmetry<1>, index_list<INDEX(data)>>& \
          vector_b,                                                    \
      const Tensor<DTYPE(data), Symmetry<1, 1>,                        \
                   index_list<change_index_up_lo<INDEX(data)>,         \
                              change_index_up_lo<INDEX(data)>>>&       \
          metric) noexcept;

GENERATE_INSTANTIATIONS(INSTANTIATE, (1, 2, 3), (double, DataVector),
                        (Frame::Grid, Frame::Inertial),
                        (SpatialIndex, SpacetimeIndex), (UpLo::Lo, UpLo::Up))

#undef INSTANTIATE
#undef INDEX
#undef UPORLO
#undef INDEXTYPE
#undef FRAME
#undef DTYPE
#undef DIM
/// \endcond

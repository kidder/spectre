// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"

#include <array>
#include <cstddef>
#include <optional>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/StdArrayHelpers.hpp"

namespace TestHelpers::evolution::dg::Actions {
namespace wave {
template <size_t Dim>
std::array<double, Dim> k() {
  if constexpr (Dim == 1) {
    return std::array{1.0};
  } else if constexpr (Dim == 2) {
    return std::array{0.6, 0.8};
  } else {
    return std::array{3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0};
  }
}

template std::array<double, 1> k<1>();
template std::array<double, 2> k<2>();
template std::array<double, 3> k<3>();

template <size_t Dim>
std::array<double, Dim> comoving_v() {
  return c_s * k<Dim>();
}

template std::array<double, 1> comoving_v<1>();
template std::array<double, 2> comoving_v<2>();
template std::array<double, 3> comoving_v<3>();

template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  if constexpr (Dim == 1) {
    return DataVector{k<1>()[0] * get<0>(x) - c_s * t};
  } else if constexpr (Dim == 2) {
    return DataVector{k<2>()[0] * get<0>(x) + k<2>()[1] * get<1>(x) - c_s * t};
  } else {
    return DataVector{k<3>()[0] * get<0>(x) + k<3>()[1] * get<1>(x) +
                      k<3>()[2] * get<2>(x) - c_s * t};
  }
}

template <size_t Dim>
DataVector f(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  return DataVector{0.01 * square(u(x, t)) + 1.0};
}

template <size_t Dim>
DataVector df(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
              const double t) {
  return DataVector{0.02 * u(x, t)};
}

template DataVector f(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                      const double t);
template DataVector f(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                      const double t);
template DataVector f(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                      const double t);

template DataVector df(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                       const double t);
template DataVector df(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                       const double t);
template DataVector df(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                       const double t);
}  // namespace wave

namespace grid {
template <Is grid_is, size_t Dim>
std::optional<tnsr::I<DataVector, Dim, Frame::Inertial>> v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  if constexpr (grid_is == Is::Stationary) {
    return std::nullopt;
  } else if constexpr (grid_is == Is::Comoving) {
    auto result = make_with_value<tnsr::I<DataVector, Dim, Frame::Inertial>>(
        x, wave::c_s);
    for (size_t d = 0; d < Dim; ++d) {
      result.get(d) *= gsl::at(wave::k<Dim>(), d);
    }
    return result;
  } else {
    tnsr::I<DataVector, Dim, Frame::Inertial> result = x;
    const double a = a_0 + a_dot * t;
    for (size_t d = 0; d < Dim; ++d) {
      result.get(d) *= a_dot / a;
    }
    return result;
  }
}

template <Is grid_is, size_t Dim>
std::optional<Scalar<DataVector>> div_v(
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  if constexpr (grid_is == Is::Stationary) {
    return std::nullopt;
  } else if constexpr (grid_is == Is::Comoving) {
    return make_with_value<Scalar<DataVector>>(x, 0.0);
  } else {
    return make_with_value<Scalar<DataVector>>(
        x, a_dot * static_cast<double>(Dim) / (a_0 + a_dot * t));
  }
}

template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>
v<Is::Stationary>(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                  const double t);
template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>> v<Is::Comoving>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<tnsr::I<DataVector, 1, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                 const double t);
template std::optional<tnsr::I<DataVector, 2, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                 const double t);
template std::optional<tnsr::I<DataVector, 3, Frame::Inertial>>
v<Is::Expanding>(const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                 const double t);

template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Stationary>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Comoving>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template std::optional<Scalar<DataVector>> div_v<Is::Expanding>(
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
}  // namespace grid

namespace {
template <size_t Dim>
DataVector n_dot_k(
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector) {
  DataVector n_dot_k = get<0>(normal_covector) * wave::k<Dim>()[0];
  for (size_t i = 1; i < Dim; ++i) {
    n_dot_k += normal_covector.get(i) * gsl::at(wave::k<Dim>(), i);
  }
  return n_dot_k;
}

template <grid::Is grid_is, size_t Dim>
DataVector n_dot_v(
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  DataVector n_dot_v =
      get<0>(normal_covector) * get<0>(grid::v<grid_is>(x, t).value());
  for (size_t i = 1; i < Dim; ++i) {
    n_dot_v += normal_covector.get(i) * grid::v<grid_is>(x, t).value().get(i);
  }
  return n_dot_v;
}
}  // namespace

template <size_t Dim>
void Var1::value(const gsl::not_null<type*> var1,
                 const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                 const double t) {
  var1->get() = wave::f(x, t);
}

template <grid::Is grid_is, size_t Dim>
void Var1::normal_dot_flux(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
    const double sign) {
  get(*normal_dot_flux) =
      sign * n_dot_k(normal_covector) * square(wave::f(x, t));
  if constexpr (grid_is != grid::Is::Stationary) {
    get(*normal_dot_flux) -=
        sign * n_dot_v<grid_is>(normal_covector, x, t) * wave::f(x, t);
  }
}

template <size_t Dim>
void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2) {
  *flux1 = var2;
}

template <grid::Is grid_is, size_t Dim>
void Var1::dt_value(const gsl::not_null<type*> dt_var1,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t) {
  dt_var1->get() = -wave::c_s * wave::df(x, t);
  if constexpr (grid_is != grid::Is::Stationary) {
    for (size_t i = 0; i < Dim; ++i) {
      dt_var1->get() += grid::v<grid_is>(x, t).value().get(i) *
                        gsl::at(wave::k<Dim>(), i) * wave::df(x, t);
    }
  }
}

template <size_t Dim>
void Var2<Dim>::value(const gsl::not_null<type*> var2,
                      const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                      const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    var2->get(d) = gsl::at(wave::k<Dim>(), d) * square(wave::f(x, t));
  }
}

template <size_t Dim>
void Var2<Dim>::flux(
    const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux2,
    const Scalar<DataVector>& var1,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& var2,
    const Scalar<DataVector>& prim_var) {
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      flux2->get(i, j) = var2.get(i) * var2.get(j) / get(var1);
      if (i == j) {
        flux2->get(i, j) += get(prim_var);
      }
    }
  }
}

template <size_t Dim>
template <grid::Is grid_is>
void Var2<Dim>::normal_dot_flux(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t,
    const double sign) {
  for (size_t i = 0; i < Dim; ++i) {
    normal_dot_flux->get(i) = sign * gsl::at(wave::k<Dim>(), i) *
                              n_dot_k(normal_covector) * cube(wave::f(x, t));
    normal_dot_flux->get(i) +=
        0.5 * sign * square(wave::c_s) * wave::f(x, t) * normal_covector.get(i);
    if constexpr (grid_is != grid::Is::Stationary) {
      normal_dot_flux->get(i) -=
          sign * n_dot_v<grid_is>(normal_covector, x, t) *
          gsl::at(wave::k<Dim>(), i) * square(wave::f(x, t));
    }
  }
}

template <size_t Dim>
template <grid::Is grid_is>
void Var2<Dim>::dt_value(const gsl::not_null<type*> dt_var2,
                         const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                         const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    dt_var2->get(d) = -2.0 * wave::c_s * gsl::at(wave::k<Dim>(), d) *
                      wave::f(x, t) * wave::df(x, t);
  }
  if constexpr (grid_is != grid::Is::Stationary) {
    for (size_t d = 0; d < Dim; ++d) {
      for (size_t i = 0; i < Dim; ++i) {
        dt_var2->get(d) += 2.0 * gsl::at(wave::k<Dim>(), d) * wave::f(x, t) *
                           grid::v<grid_is>(x, t).value().get(i) *
                           gsl::at(wave::k<Dim>(), i) * wave::df(x, t);
      }
    }
  }
}

void TempVar::value(const gsl::not_null<type*> temp_var,
                    const Scalar<DataVector>& var1) {
  get(*temp_var) = 0.5 * square(wave::c_s) * get(var1);
}

template <size_t Dim>
void TempVar::value(const gsl::not_null<type*> temp_var,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t) {
  get(*temp_var) = 0.5 * square(wave::c_s) * wave::f(x, t);
}

void PrimVar::value(const gsl::not_null<type*> prim_var,
                    const Scalar<DataVector>& var1) {
  get(*prim_var) = 0.5 * square(wave::c_s) * get(var1);
}

template <size_t Dim>
void PrimVar::value(const gsl::not_null<type*> prim_var,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t) {
  get(*prim_var) = 0.5 * square(wave::c_s) * wave::f(x, t);
}

template <size_t Dim>
void Source1::value(const gsl::not_null<type*> source1,
                    const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                    const double t) {
  source1->get() = wave::df(x, t) * (2.0 * wave::f(x, t) - wave::c_s);
}

template <size_t Dim>
void Source2<Dim>::value(const gsl::not_null<type*> source2,
                         const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                         const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    source2->get(d) =
        gsl::at(wave::k<Dim>(), d) * wave::df(x, t) *
        (3.0 * square(wave::f(x, t)) - 2.0 * wave::c_s * wave::f(x, t) +
         0.5 * square(wave::c_s));
  }
}

template struct Var2<1>;
template struct Var2<2>;
template struct Var2<3>;

template void Var1::value(const gsl::not_null<type*> var1,
                          const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                          const double t);
template void Var1::value(const gsl::not_null<type*> var1,
                          const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                          const double t);
template void Var1::value(const gsl::not_null<type*> var1,
                          const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                          const double t);

template void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, 1, Frame::Inertial>*> flux1,
    const tnsr::I<DataVector, 1, Frame::Inertial>& var2);
template void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, 2, Frame::Inertial>*> flux1,
    const tnsr::I<DataVector, 2, Frame::Inertial>& var2);
template void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> flux1,
    const tnsr::I<DataVector, 3, Frame::Inertial>& var2);

template void Var1::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var1::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<Scalar<DataVector>*> normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);

template void Var2<1>::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<tnsr::I<DataVector, 1, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<2>::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<tnsr::I<DataVector, 2, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<3>::normal_dot_flux<grid::Is::Stationary>(
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<1>::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<tnsr::I<DataVector, 1, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<2>::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<tnsr::I<DataVector, 2, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<3>::normal_dot_flux<grid::Is::Comoving>(
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<1>::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<tnsr::I<DataVector, 1, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 1, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<2>::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<tnsr::I<DataVector, 2, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 2, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t,
    const double sign);
template void Var2<3>::normal_dot_flux<grid::Is::Expanding>(
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*>
        normal_dot_flux,
    const tnsr::i<DataVector, 3, Frame::Inertial>& normal_covector,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t,
    const double sign);

template void Var1::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var1::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var1,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);

template void Var2<1>::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var2<2>::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var2<3>::dt_value<grid::Is::Stationary>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template void Var2<1>::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var2<2>::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var2<3>::dt_value<grid::Is::Comoving>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);
template void Var2<1>::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x, const double t);
template void Var2<2>::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x, const double t);
template void Var2<3>::dt_value<grid::Is::Expanding>(
    const gsl::not_null<type*> dt_var2,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x, const double t);

template void TempVar::value(const gsl::not_null<type*> temp_var,
                             const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                             const double t);
template void TempVar::value(const gsl::not_null<type*> temp_var,
                             const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                             const double t);
template void TempVar::value(const gsl::not_null<type*> temp_var,
                             const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                             const double t);

template void PrimVar::value(const gsl::not_null<type*> prim_var,
                             const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                             const double t);
template void PrimVar::value(const gsl::not_null<type*> prim_var,
                             const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                             const double t);
template void PrimVar::value(const gsl::not_null<type*> prim_var,
                             const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                             const double t);

template void Source1::value(const gsl::not_null<type*> source1,
                             const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                             const double t);
template void Source1::value(const gsl::not_null<type*> source1,
                             const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                             const double t);
template void Source1::value(const gsl::not_null<type*> source1,
                             const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                             const double t);

template struct Source2<1>;
template struct Source2<2>;
template struct Source2<3>;

}  // namespace TestHelpers::evolution::dg::Actions

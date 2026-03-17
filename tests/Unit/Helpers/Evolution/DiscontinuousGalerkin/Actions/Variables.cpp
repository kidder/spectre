// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Evolution/DiscontinuousGalerkin/Actions/Variables.hpp"

namespace TestHelpers::evolution::dg::Actions {
namespace wave {
template <>
constexpr std::array<double, 1> k<1> = std::array{1.0};

template <>
constexpr std::array<double, 2> k<2> = std::array{0.6, 0.8};

template <>
constexpr std::array<double, 3> k<3> =
    std::array{3.0 / 13.0, 4.0 / 13.0, 12.0 / 13.0};

template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t) {
  if constexpr (Dim == 1) {
    return DataVector{k<1>[0] * get<0>(x) - c_s * t};
  } else if constexpr (Dim == 2) {
    return DataVector{k<2>[0] * get<0>(x) + k<2>[1] * get<1>(x) - c_s * t};
  } else {
    return DataVector{k<3>[0] * get<0>(x) + k<3>[1] * get<1>(x) +
                      k<3>[2] * get<2>(x) - c_s * t};
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

template <size_t Dim>
void Var1::value(const gsl::not_null<type*> var1,
                 const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                 const double t) {
  var1->get() = wave::f(x, t);
}

template <size_t Dim>
void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  for (size_t i = 0; i < Dim; ++i) {
    flux->get(i) = gsl::at(wave::k<Dim>, i) * square(wave::f(x, t));
  }
}

template <size_t Dim>
void Var1::dt_value(const gsl::not_null<type*> dt_var1,
                     const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                     const double t) {
  dt_var1->get() = -wave::c_s * wave::df(x, t);
}

template <size_t Dim>
void Var2<Dim>::value(const gsl::not_null<type*> var2,
                      const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                      const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    var2->get(d) = gsl::at(wave::k<Dim>, d) * square(wave::f(x, t));
  }
}

template <size_t Dim>
void Var2<Dim>::flux(
    const gsl::not_null<tnsr::IJ<DataVector, Dim, Frame::Inertial>*> flux,
    const tnsr::I<DataVector, Dim, Frame::Inertial>& x, const double t) {
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      flux->get(i, j) = gsl::at(wave::k<Dim>, i) * gsl::at(wave::k<Dim>, j) *
                        cube(wave::f(x, t));
      if (i == j) {
        flux->get(i, j) += 0.5 * square(wave::c_s) * wave::f(x, t);
      }
    }
  }
}

template <size_t Dim>
void Var2<Dim>::dt_value(const gsl::not_null<type*> dt_var2,
                         const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                         const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    dt_var2->get(d) = -2.0 * wave::c_s * gsl::at(wave::k<Dim>, d) *
                      wave::f(x, t) * wave::df(x, t);
  }
}

template <size_t Dim>
void Var3<Dim>::value(const gsl::not_null<type*> var3,
                      const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                      const double t) {
  for (size_t d = 0; d < Dim; ++d) {
    var3->get(d) = gsl::at(wave::k<Dim>, d) * wave::f(x, t);
  }
}

void Var4::value(const gsl::not_null<type*> var4,
                 const Scalar<DataVector>& var1) {
  get(*var4) = 0.5 * square(wave::c_s) * get(var1);
}

template <size_t Dim>
void Var4::value(const gsl::not_null<type*> var4,
                 const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
                 const double t) {
  get(*var4) = 0.5 * square(wave::c_s) * wave::f(x, t);
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
        gsl::at(wave::k<Dim>, d) * wave::df(x, t) *
        (3.0 * square(wave::f(x, t)) - 2.0 * wave::c_s * wave::f(x, t) +
         0.5 * square(wave::c_s));
  }
}

template struct Var2<1>;
template struct Var2<2>;
template struct Var2<3>;
template struct Var3<1>;
template struct Var3<2>;
template struct Var3<3>;

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
    const gsl::not_null<tnsr::I<DataVector, 1, Frame::Inertial>*> flux,
    const tnsr::I<DataVector, 1, Frame::Inertial>& x,
    const double t);
template void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, 2, Frame::Inertial>*> flux,
    const tnsr::I<DataVector, 2, Frame::Inertial>& x,
    const double t);
template void Var1::flux(
    const gsl::not_null<tnsr::I<DataVector, 3, Frame::Inertial>*> flux,
    const tnsr::I<DataVector, 3, Frame::Inertial>& x,
    const double t);

template void Var1::dt_value(const gsl::not_null<type*> dt_var1,
                          const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                          const double t);
template void Var1::dt_value(const gsl::not_null<type*> dt_var1,
   const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                          const double t);
template void Var1::dt_value(const gsl::not_null<type*> dt_var1,
                          const tnsr::I<DataVector, 3, Frame::Inertial>& x,
                          const double t);

template void Var4::value(const gsl::not_null<type*> var4,
                          const tnsr::I<DataVector, 1, Frame::Inertial>& x,
                          const double t);
template void Var4::value(const gsl::not_null<type*> var4,
                          const tnsr::I<DataVector, 2, Frame::Inertial>& x,
                          const double t);
template void Var4::value(const gsl::not_null<type*> var4,
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

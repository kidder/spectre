// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>

#include "DataStructures/Tensor/TypeAliases.hpp"

/// \cond
class DataVector;
/// \endcond

/// Useful helpers for constructing a wave-like solution
namespace TestHelpers::evolution::wave {
/// The wave speed \f$c_s = \frac{1}{2}\f$
static constexpr double c_s = 0.5;

/// \brief The wave vector \f$\vec{k}\f$
///
/// \details If \f$d = 1\f$, \f$\vec{k} = (1.0)\f$.
/// If \f$d = 2\f$, \f$\vec{k} = (\frac{3}{5}, \frac{4}{5})\f$.
/// If \f$d = 3\f$, \f$\vec{k} = (\frac{3}{13}, \frac{4}{13}, \frac{12}{13})\f$.
template <size_t Dim>
std::array<double, Dim> k();

/// The velocity comoving with the wave \f$\vec{v} = c_s \vec{k}\f$
template <size_t Dim>
std::array<double, Dim> comoving_v();

/// The variables are given as functions of \f$u = k_i x^i - t\f$
template <size_t Dim>
DataVector u(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

/// The profile function \f$f = 1 + \frac{1}{100} u^2\f$
template <size_t Dim>
DataVector f(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
             const double t);

/// The derivative of the profile function \f$\frac{df}{du} = \frac{1}{50} u\f$
template <size_t Dim>
DataVector df(const tnsr::I<DataVector, Dim, Frame::Inertial>& x,
              const double t);
}  // namespace TestHelpers::evolution::wave

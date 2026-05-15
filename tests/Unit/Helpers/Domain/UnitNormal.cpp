// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Domain/UnitNormal.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/FaceNormal.hpp"

template <size_t Dim>
tnsr::i<DataVector, Dim, Frame::Inertial> unit_normal(
    const Mesh<Dim - 1>& mesh, const Direction<Dim>& direction,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  auto n = unnormalized_face_normal(mesh, element_map, coordinate_map, t,
                                    functions_of_time, direction);
  const auto magnitude_of_n = magnitude(n);
  for (size_t d = 0; d < Dim; ++d) {
    n.get(d) /= get(magnitude_of_n);
  }
  return n;
}

template tnsr::i<DataVector, 1, Frame::Inertial> unit_normal(
    const Mesh<0>& mesh, const Direction<1>& direction,
    const ElementMap<1, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 1>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template tnsr::i<DataVector, 2, Frame::Inertial> unit_normal(
    const Mesh<1>& mesh, const Direction<2>& direction,
    const ElementMap<2, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 2>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template tnsr::i<DataVector, 3, Frame::Inertial> unit_normal(
    const Mesh<2>& mesh, const Direction<3>& direction,
    const ElementMap<3, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

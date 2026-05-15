// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/Domain/NormalCovectorAndMagnitude.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/FaceNormal.hpp"
#include "Evolution/DiscontinuousGalerkin/NormalVectorTags.hpp"

template <size_t Dim>
Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                     ::evolution::dg::Tags::NormalCovector<Dim>>>
normal_covector_and_magnitude(
    const Mesh<Dim - 1>& mesh, const Direction<Dim>& direction,
    const ElementMap<Dim, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, Dim>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time) {
  Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                       ::evolution::dg::Tags::NormalCovector<Dim>>>
      result{mesh.number_of_grid_points()};
  auto& n = get<::evolution::dg::Tags::NormalCovector<Dim>>(result);
  n = unnormalized_face_normal(mesh, element_map, coordinate_map, t,
                               functions_of_time, direction);
  auto& magnitude_of_n = get<::evolution::dg::Tags::MagnitudeOfNormal>(result);
  magnitude_of_n = magnitude(n);
  for (size_t d = 0; d < Dim; ++d) {
    n.get(d) /= get(magnitude_of_n);
  }
  return result;
}

template Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                              ::evolution::dg::Tags::NormalCovector<1>>>
normal_covector_and_magnitude(
    const Mesh<0>& mesh, const Direction<1>& direction,
    const ElementMap<1, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 1>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                              ::evolution::dg::Tags::NormalCovector<2>>>
normal_covector_and_magnitude(
    const Mesh<1>& mesh, const Direction<2>& direction,
    const ElementMap<2, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 2>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

template Variables<tmpl::list<::evolution::dg::Tags::MagnitudeOfNormal,
                              ::evolution::dg::Tags::NormalCovector<3>>>
normal_covector_and_magnitude(
    const Mesh<2>& mesh, const Direction<3>& direction,
    const ElementMap<3, Frame::Grid>& element_map,
    const domain::CoordinateMapBase<Frame::Grid, Frame::Inertial, 3>&
        coordinate_map,
    const double t,
    const std::unordered_map<
        std::string, std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime>>&
        functions_of_time);

// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//

#include "custom_utilities/embedded_beam_utilities.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

#include <algorithm>

namespace Kratos
{

void EmbeddedHostInterpolation::Initialize(const Geometry<Node>&          rEmbeddedGeometry,
                                           GlobalPointersVector<Element>& rCandidateHostElements,
                                           double                         Tolerance)
{
    KRATOS_TRY

    const auto number_of_embedded_nodes = rEmbeddedGeometry.PointsNumber();
    mInterpolationNodes.clear();
    mNodeIndices.assign(number_of_embedded_nodes, {});
    mWeights.assign(number_of_embedded_nodes, {});
    mHostElements.assign(number_of_embedded_nodes, nullptr);

    for (std::size_t i = 0; i < number_of_embedded_nodes; ++i) {
        const auto& r_position = rEmbeddedGeometry[i].Coordinates();

        for (auto it = rCandidateHostElements.ptr_begin(); it != rCandidateHostElements.ptr_end(); ++it) {
            auto*       p_host          = it->get();
            const auto& r_host_geometry = p_host->GetGeometry();

            array_1d<double, 3> local_coordinates;
            if (!r_host_geometry.IsInside(r_position, local_coordinates, Tolerance)) continue;

            Vector shape_function_values;
            r_host_geometry.ShapeFunctionsValues(shape_function_values, local_coordinates);
            for (std::size_t k = 0; k < r_host_geometry.PointsNumber(); ++k) {
                mNodeIndices[i].push_back(AddInterpolationNode(r_host_geometry(k)));
                mWeights[i].push_back(shape_function_values[k]);
            }
            mHostElements[i] = p_host;
            break;
        }

        if (!mHostElements[i]) {
            mNodeIndices[i].push_back(AddInterpolationNode(rEmbeddedGeometry(i)));
            mWeights[i].push_back(1.0);
        }
    }

    KRATOS_CATCH("")
}

std::size_t EmbeddedHostInterpolation::AddInterpolationNode(const Node::Pointer& rpNode)
{
    const auto it =
        std::find_if(mInterpolationNodes.begin(), mInterpolationNodes.end(),
                     [&rpNode](const auto& rpOther) { return rpOther->Id() == rpNode->Id(); });
    if (it != mInterpolationNodes.end()) {
        return static_cast<std::size_t>(std::distance(mInterpolationNodes.begin(), it));
    }

    mInterpolationNodes.push_back(rpNode);
    return mInterpolationNodes.size() - 1;
}

bool EmbeddedHostInterpolation::IsInitialized() const { return !mInterpolationNodes.empty(); }

const std::vector<Node::Pointer>& EmbeddedHostInterpolation::GetInterpolationNodes() const
{
    return mInterpolationNodes;
}

const std::vector<std::size_t>& EmbeddedHostInterpolation::GetInterpolationNodeIndices(std::size_t EmbeddedNodeIndex) const
{
    return mNodeIndices.at(EmbeddedNodeIndex);
}

const std::vector<double>& EmbeddedHostInterpolation::GetInterpolationWeights(std::size_t EmbeddedNodeIndex) const
{
    return mWeights.at(EmbeddedNodeIndex);
}

Element* EmbeddedHostInterpolation::GetHostElement(std::size_t EmbeddedNodeIndex) const
{
    return mHostElements.at(EmbeddedNodeIndex);
}

bool EmbeddedHostInterpolation::IsEmbedded(std::size_t EmbeddedNodeIndex) const
{
    return GetHostElement(EmbeddedNodeIndex) != nullptr;
}

bool EmbeddedHostInterpolation::IsInterpolatedFromItself(const Node& rEmbeddedNode, std::size_t EmbeddedNodeIndex) const
{
    const auto& r_indices = GetInterpolationNodeIndices(EmbeddedNodeIndex);
    const auto& r_weights = GetInterpolationWeights(EmbeddedNodeIndex);
    for (std::size_t k = 0; k < r_indices.size(); ++k) {
        if (mInterpolationNodes[r_indices[k]]->Id() == rEmbeddedNode.Id() &&
            std::abs(r_weights[k] - 1.0) < 1.0e-12) {
            return true;
        }
    }
    return false;
}

array_1d<double, 3> EmbeddedHostInterpolation::Interpolate(const Variable<array_1d<double, 3>>& rVariable,
                                                           std::size_t EmbeddedNodeIndex,
                                                           IndexType   Step) const
{
    const auto& r_indices = GetInterpolationNodeIndices(EmbeddedNodeIndex);
    const auto& r_weights = GetInterpolationWeights(EmbeddedNodeIndex);

    array_1d<double, 3> result = ZeroVector(3);
    for (std::size_t k = 0; k < r_indices.size(); ++k) {
        noalias(result) +=
            r_weights[k] * mInterpolationNodes[r_indices[k]]->FastGetSolutionStepValue(rVariable, Step);
    }
    return result;
}

BoundedMatrix<double, 3, 3> EmbeddedBeamUtilities::CalculateLocalAxes(const Geometry<Node>& rBeamGeometry,
                                                                      std::size_t Dimension)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rBeamGeometry.PointsNumber() < 2)
        << "A beam geometry needs at least two points to determine its local axes" << std::endl;

    array_1d<double, 3> t = rBeamGeometry[1].GetInitialPosition().Coordinates() -
                            rBeamGeometry[0].GetInitialPosition().Coordinates();
    const auto length = norm_2(t);
    KRATOS_ERROR_IF(length <= 0.0) << "The beam geometry has a zero length" << std::endl;
    t /= length;

    array_1d<double, 3> n = ZeroVector(3);
    array_1d<double, 3> m = ZeroVector(3);
    if (Dimension == 2) {
        n[0] = -t[1];
        n[1] = t[0];
        m[2] = 1.0;
    } else {
        if (rBeamGeometry.Has(LOCAL_AXIS_2)) {
            noalias(n) = rBeamGeometry.GetValue(LOCAL_AXIS_2);
        } else {
            n[1] = 1.0;
        }

        // When the auxiliary vector is (anti-)parallel to the beam axis, use the global z axis instead
        if (std::abs(inner_prod(t, n)) >= (1.0 - 1.0e-8) * norm_2(n)) {
            n.clear();
            n[2] = 1.0;
        }

        // Gram-Schmidt orthogonalization
        noalias(n) -= inner_prod(t, n) * t;
        n /= norm_2(n);
        noalias(m) = MathUtils<double>::CrossProduct(t, n);
    }

    BoundedMatrix<double, 3, 3> result;
    for (std::size_t i = 0; i < 3; ++i) {
        result(0, i) = t[i];
        result(1, i) = n[i];
        result(2, i) = m[i];
    }
    return result;

    KRATOS_CATCH("")
}

double EmbeddedBeamUtilities::CalculateReferenceLength(const Geometry<Node>& rBeamGeometry)
{
    return norm_2(rBeamGeometry[1].GetInitialPosition().Coordinates() -
                  rBeamGeometry[0].GetInitialPosition().Coordinates());
}

} // namespace Kratos

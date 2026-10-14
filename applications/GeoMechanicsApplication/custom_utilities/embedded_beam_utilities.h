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

#pragma once

#include "containers/global_pointers_vector.h"
#include "geometries/geometry.h"
#include "includes/element.h"
#include "includes/kratos_export_api.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

#include <vector>

namespace Kratos
{

/**
 * @brief Interpolation of the translations of the nodes of an embedded line element from the
 * displacement field of the (solid) host elements in which these nodes are located.
 * @details Following Sadek & Shahrour (2004), "A three dimensional embedded beam element for
 * reinforced geomaterials", Int. J. Numer. Anal. Meth. Geomech. 28:931-946, Eq. (5), the
 * translation of embedded node i reads
 *     U_i = N^s(X_i) U^s
 * where N^s are the shape functions of the host element evaluated at the position X_i of the
 * embedded node and U^s are the nodal displacements of that host element. Every node that
 * contributes to the interpolation is stored only once, such that it can be used to set up the
 * degrees of freedom of the embedded element. An embedded node that is not located inside any of
 * the candidate host elements is 'free': its translation is given by its own displacement.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) EmbeddedHostInterpolation
{
public:
    void Initialize(const Geometry<Node>&          rEmbeddedGeometry,
                    GlobalPointersVector<Element>& rCandidateHostElements,
                    double                         Tolerance);

    [[nodiscard]] bool                              IsInitialized() const;
    [[nodiscard]] const std::vector<Node::Pointer>& GetInterpolationNodes() const;
    [[nodiscard]] const std::vector<std::size_t>& GetInterpolationNodeIndices(std::size_t EmbeddedNodeIndex) const;
    [[nodiscard]] const std::vector<double>& GetInterpolationWeights(std::size_t EmbeddedNodeIndex) const;

    /// Returns the host element of the given embedded node, or a nullptr when the node is free
    [[nodiscard]] Element* GetHostElement(std::size_t EmbeddedNodeIndex) const;
    [[nodiscard]] bool     IsEmbedded(std::size_t EmbeddedNodeIndex) const;

    /// Returns true when the translation of the embedded node is given by the node itself, i.e.
    /// when it is free or when it coincides with (is) a node of its host element
    [[nodiscard]] bool IsInterpolatedFromItself(const Node& rEmbeddedNode, std::size_t EmbeddedNodeIndex) const;

    [[nodiscard]] array_1d<double, 3> Interpolate(const Variable<array_1d<double, 3>>& rVariable,
                                                  std::size_t EmbeddedNodeIndex,
                                                  IndexType   Step = 0) const;

private:
    std::size_t AddInterpolationNode(const Node::Pointer& rpNode);

    std::vector<Node::Pointer>            mInterpolationNodes;
    std::vector<std::vector<std::size_t>> mNodeIndices;
    std::vector<std::vector<double>>      mWeights;
    std::vector<Element*>                 mHostElements;
};

class KRATOS_API(GEO_MECHANICS_APPLICATION) EmbeddedBeamUtilities
{
public:
    /**
     * @brief Calculates the local axes of a straight 2-noded beam in its reference configuration
     * @return Matrix of which the rows are the unit vectors of the local axes (x' along the beam
     * axis), i.e. the rotation matrix from global to local axes. In 2D, the local y' axis is
     * obtained by rotating x' by 90 degrees counter-clockwise and z' coincides with the global z
     * axis. In 3D, the local y' axis is (the orthogonal part of) LOCAL_AXIS_2 when given on the
     * geometry, else (the orthogonal part of) the global y axis, else the global z axis.
     */
    [[nodiscard]] static BoundedMatrix<double, 3, 3> CalculateLocalAxes(const Geometry<Node>& rBeamGeometry,
                                                                        std::size_t Dimension);

    [[nodiscard]] static double CalculateReferenceLength(const Geometry<Node>& rBeamGeometry);
};

} // namespace Kratos

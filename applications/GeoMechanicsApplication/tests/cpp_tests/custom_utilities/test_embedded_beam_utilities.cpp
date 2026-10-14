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
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "utilities/math_utils.h"

#include <cmath>

namespace
{

using namespace Kratos;

Matrix MakeMatrix(const std::vector<std::vector<double>>& rRows)
{
    Matrix result(rRows.size(), rRows.front().size());
    for (std::size_t i = 0; i < rRows.size(); ++i) {
        for (std::size_t j = 0; j < rRows[i].size(); ++j) {
            result(i, j) = rRows[i][j];
        }
    }
    return result;
}

Matrix CalculateLocalAxesOfLine3D(const std::vector<double>& rStart, const std::vector<double>& rEnd)
{
    const auto geometry = Line3D2<Node>(Kratos::make_intrusive<Node>(1, rStart[0], rStart[1], rStart[2]),
                                        Kratos::make_intrusive<Node>(2, rEnd[0], rEnd[1], rEnd[2]));
    return EmbeddedBeamUtilities::CalculateLocalAxes(geometry, 3);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(EmbeddedBeamUtilities_CalculatesLocalAxesOfInclinedBeamIn2D, KratosGeoMechanicsFastSuite)
{
    const auto geometry = Line2D2<Node>(Kratos::make_intrusive<Node>(1, 1.0, 1.0, 0.0),
                                        Kratos::make_intrusive<Node>(2, 4.0, 5.0, 0.0));

    const Matrix local_axes = EmbeddedBeamUtilities::CalculateLocalAxes(geometry, 2);

    KRATOS_EXPECT_MATRIX_NEAR(local_axes, MakeMatrix({{0.6, 0.8, 0.0}, {-0.8, 0.6, 0.0}, {0.0, 0.0, 1.0}}), 1.0e-12)
    KRATOS_EXPECT_NEAR(EmbeddedBeamUtilities::CalculateReferenceLength(geometry), 5.0, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedBeamUtilities_CalculatesLocalAxesOfBeamsIn3D, KratosGeoMechanicsFastSuite)
{
    // Along x: the local axes coincide with the global axes
    KRATOS_EXPECT_MATRIX_NEAR(CalculateLocalAxesOfLine3D({0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}),
                              MakeMatrix({{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}), 1.0e-12)

    // Vertical: y' is the global y axis
    KRATOS_EXPECT_MATRIX_NEAR(CalculateLocalAxesOfLine3D({0.0, 0.0, 0.0}, {0.0, 0.0, -3.0}),
                              MakeMatrix({{0.0, 0.0, -1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}}), 1.0e-12)

    // Along +y and -y: y' is the global z axis
    KRATOS_EXPECT_MATRIX_NEAR(CalculateLocalAxesOfLine3D({0.0, 0.0, 0.0}, {0.0, 2.0, 0.0}),
                              MakeMatrix({{0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}}), 1.0e-12)
    KRATOS_EXPECT_MATRIX_NEAR(CalculateLocalAxesOfLine3D({0.0, 2.0, 0.0}, {0.0, 0.0, 0.0}),
                              MakeMatrix({{0.0, -1.0, 0.0}, {0.0, 0.0, 1.0}, {-1.0, 0.0, 0.0}}), 1.0e-12)

    // Arbitrarily inclined: the local axes form an orthonormal, right-handed basis
    const auto local_axes = CalculateLocalAxesOfLine3D({0.3, 0.4, 0.2}, {1.5, 1.2, 1.7});
    KRATOS_EXPECT_MATRIX_NEAR(Matrix(prod(local_axes, trans(local_axes))), IdentityMatrix(3), 1.0e-12)
    KRATOS_EXPECT_NEAR(MathUtils<double>::Det(local_axes), 1.0, 1.0e-12);
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedBeamUtilities_UsesLocalAxis2OfGeometryWhenGiven, KratosGeoMechanicsFastSuite)
{
    auto                geometry = Line3D2<Node>(Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0),
                                                 Kratos::make_intrusive<Node>(2, 0.0, 0.0, 2.0));
    array_1d<double, 3> local_axis_2;
    local_axis_2[0] = 1.0;
    local_axis_2[1] = 1.0;
    local_axis_2[2] = 0.5;
    geometry.SetValue(LOCAL_AXIS_2, local_axis_2);

    const auto c = 1.0 / std::sqrt(2.0);
    KRATOS_EXPECT_MATRIX_NEAR(EmbeddedBeamUtilities::CalculateLocalAxes(geometry, 3),
                              MakeMatrix({{0.0, 0.0, 1.0}, {c, c, 0.0}, {-c, c, 0.0}}), 1.0e-12)
}

KRATOS_TEST_CASE_IN_SUITE(EmbeddedHostInterpolation_CollectsUniqueHostNodesAndHandlesFreeNodes, KratosGeoMechanicsFastSuite)
{
    // Two quadrilaterals sharing nodes 2 and 5
    std::vector<Node::Pointer> nodes;
    for (const auto& [id, x, y] : std::vector<std::tuple<std::size_t, double, double>>{
             {1, 0.0, 0.0}, {2, 1.0, 0.0}, {3, 2.0, 0.0}, {4, 0.0, 1.0}, {5, 1.0, 1.0}, {6, 2.0, 1.0}}) {
        nodes.push_back(Kratos::make_intrusive<Node>(id, x, y, 0.0));
    }
    auto p_host_1 = Kratos::make_intrusive<Element>(
        1, Kratos::make_shared<Quadrilateral2D4<Node>>(nodes[0], nodes[1], nodes[4], nodes[3]));
    auto p_host_2 = Kratos::make_intrusive<Element>(
        2, Kratos::make_shared<Quadrilateral2D4<Node>>(nodes[1], nodes[2], nodes[5], nodes[4]));
    GlobalPointersVector<Element> hosts;
    hosts.push_back(GlobalPointer<Element>(p_host_1.get()));
    hosts.push_back(GlobalPointer<Element>(p_host_2.get()));

    // The first beam node is inside the first host, the second one is outside both hosts
    const auto embedded_geometry = Line2D2<Node>(Kratos::make_intrusive<Node>(11, 0.25, 0.5, 0.0),
                                                 Kratos::make_intrusive<Node>(12, 3.0, 0.5, 0.0));
    EmbeddedHostInterpolation interpolation;
    interpolation.Initialize(embedded_geometry, hosts, 1.0e-6);

    KRATOS_EXPECT_TRUE(interpolation.IsEmbedded(0))
    KRATOS_EXPECT_EQ(interpolation.GetHostElement(0)->Id(), 1u);
    KRATOS_EXPECT_FALSE(interpolation.IsEmbedded(1))
    KRATOS_EXPECT_TRUE(interpolation.IsInterpolatedFromItself(embedded_geometry[1], 1))

    // Nodes 1, 2, 5 and 4 of the host, followed by the free beam node 12
    std::vector<std::size_t> node_ids;
    for (const auto& rp_node : interpolation.GetInterpolationNodes()) {
        node_ids.push_back(rp_node->Id());
    }
    KRATOS_EXPECT_EQ(node_ids, (std::vector<std::size_t>{1, 2, 5, 4, 12}));
    const auto&               r_weights        = interpolation.GetInterpolationWeights(0);
    const std::vector<double> expected_weights = {0.375, 0.125, 0.125, 0.375};
    ASSERT_EQ(r_weights.size(), expected_weights.size());
    for (std::size_t i = 0; i < r_weights.size(); ++i) {
        KRATOS_EXPECT_NEAR(r_weights[i], expected_weights[i], 1.0e-12);
    }
}

} // namespace Kratos::Testing

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

#include "containers/model.h"
#include "custom_processes/find_embedded_beam_host_elements_process.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/kratos_parameters.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <algorithm>

namespace
{

using namespace Kratos;

std::vector<std::size_t> GetHostIds(const Element& rElement)
{
    std::vector<std::size_t> result;
    const auto&              r_hosts = rElement.GetValue(NEIGHBOUR_ELEMENTS);
    for (auto it = r_hosts.ptr_begin(); it != r_hosts.ptr_end(); ++it) {
        result.push_back(it->get()->Id());
    }
    std::sort(result.begin(), result.end());
    return result;
}

void AddElement(ModelPart& rModelPart, std::size_t Id, const Geometry<Node>::Pointer& rpGeometry)
{
    rModelPart.AddElement(Kratos::make_intrusive<Element>(Id, rpGeometry));
}

ModelPart& CreateMainModelPart(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.CreateSubModelPart("soil");
    r_result.CreateSubModelPart("beam");
    return r_result;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(FindEmbeddedBeamHostElementsProcess_AssignsSolidElementsContainingTheBeamNodes,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_main = CreateMainModelPart(model);
    auto& r_soil = r_main.GetSubModelPart("soil");
    auto& r_beam = r_main.GetSubModelPart("beam");

    // Two quadrilaterals (Ids 1 and 2) next to each other and a line element (Id 3) on top of them
    for (const auto& [id, x, y] : std::vector<std::tuple<std::size_t, double, double>>{
             {1, 0.0, 0.0}, {2, 1.0, 0.0}, {3, 2.0, 0.0}, {4, 0.0, 1.0}, {5, 1.0, 1.0}, {6, 2.0, 1.0}}) {
        r_soil.CreateNewNode(id, x, y, 0.0);
    }
    AddElement(r_soil, 1,
               Kratos::make_shared<Quadrilateral2D4<Node>>(r_soil.pGetNode(1), r_soil.pGetNode(2),
                                                           r_soil.pGetNode(5), r_soil.pGetNode(4)));
    AddElement(r_soil, 2,
               Kratos::make_shared<Quadrilateral2D4<Node>>(r_soil.pGetNode(2), r_soil.pGetNode(3),
                                                           r_soil.pGetNode(6), r_soil.pGetNode(5)));
    AddElement(r_soil, 3, Kratos::make_shared<Line2D2<Node>>(r_soil.pGetNode(4), r_soil.pGetNode(5)));

    // Beam nodes inside the first quadrilateral (11, 14), inside the second one (12), on a shared
    // node of both (5) and outside the soil (13)
    r_beam.CreateNewNode(11, 0.5, 0.5, 0.0);
    r_beam.CreateNewNode(12, 1.5, 0.5, 0.0);
    r_beam.CreateNewNode(13, 2.5, 0.5, 0.0);
    r_beam.CreateNewNode(14, 0.7, 0.2, 0.0);
    r_beam.AddNode(r_soil.pGetNode(5));
    AddElement(r_beam, 21, Kratos::make_shared<Line2D2<Node>>(r_beam.pGetNode(11), r_beam.pGetNode(12)));
    AddElement(r_beam, 22, Kratos::make_shared<Line2D2<Node>>(r_beam.pGetNode(12), r_beam.pGetNode(13)));
    AddElement(r_beam, 23, Kratos::make_shared<Line2D2<Node>>(r_beam.pGetNode(11), r_beam.pGetNode(14)));
    AddElement(r_beam, 24, Kratos::make_shared<Line2D2<Node>>(r_beam.pGetNode(14), r_beam.pGetNode(5)));

    FindEmbeddedBeamHostElementsProcess process(
        model, Parameters{R"({"model_part_name": "Main.beam", "host_model_part_name_list": ["Main.soil"]})"});
    process.ExecuteInitialize();

    KRATOS_EXPECT_EQ(GetHostIds(r_main.GetElement(21)), (std::vector<std::size_t>{1, 2}));
    KRATOS_EXPECT_EQ(GetHostIds(r_main.GetElement(22)), (std::vector<std::size_t>{2}));
    KRATOS_EXPECT_EQ(GetHostIds(r_main.GetElement(23)), (std::vector<std::size_t>{1}));
    KRATOS_EXPECT_EQ(GetHostIds(r_main.GetElement(24)), (std::vector<std::size_t>{1}));

    // The soil elements themselves did not get any host elements assigned
    KRATOS_EXPECT_FALSE(r_main.GetElement(1).Has(NEIGHBOUR_ELEMENTS))
}

KRATOS_TEST_CASE_IN_SUITE(FindEmbeddedBeamHostElementsProcess_FindsHostsOfNodesOnElementFacesInLargerMesh,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_main = CreateMainModelPart(model);
    auto& r_soil = r_main.GetSubModelPart("soil");
    auto& r_beam = r_main.GetSubModelPart("beam");

    // A structured mesh of 4 x 4 x 4 hexahedra of 1 x 1 x 1
    constexpr std::size_t number_of_divisions = 4;
    const auto            node_id             = [](std::size_t I, std::size_t J, std::size_t K) {
        return 1 + I + (number_of_divisions + 1) * (J + (number_of_divisions + 1) * K);
    };
    for (std::size_t k = 0; k <= number_of_divisions; ++k) {
        for (std::size_t j = 0; j <= number_of_divisions; ++j) {
            for (std::size_t i = 0; i <= number_of_divisions; ++i) {
                r_soil.CreateNewNode(node_id(i, j, k), static_cast<double>(i),
                                     static_cast<double>(j), static_cast<double>(k));
            }
        }
    }
    std::size_t element_id = 1;
    for (std::size_t k = 0; k < number_of_divisions; ++k) {
        for (std::size_t j = 0; j < number_of_divisions; ++j) {
            for (std::size_t i = 0; i < number_of_divisions; ++i) {
                AddElement(
                    r_soil, element_id++,
                    Kratos::make_shared<Hexahedra3D8<Node>>(
                        r_soil.pGetNode(node_id(i, j, k)), r_soil.pGetNode(node_id(i + 1, j, k)),
                        r_soil.pGetNode(node_id(i + 1, j + 1, k)), r_soil.pGetNode(node_id(i, j + 1, k)),
                        r_soil.pGetNode(node_id(i, j, k + 1)), r_soil.pGetNode(node_id(i + 1, j, k + 1)),
                        r_soil.pGetNode(node_id(i + 1, j + 1, k + 1)),
                        r_soil.pGetNode(node_id(i, j + 1, k + 1))));
            }
        }
    }

    // An inclined beam from a point on the top surface down into the mesh
    r_beam.CreateNewNode(1001, 1.5, 2.5, 4.0);
    r_beam.CreateNewNode(1002, 2.25, 2.0, 1.7);
    AddElement(r_beam, 1001, Kratos::make_shared<Line3D2<Node>>(r_beam.pGetNode(1001), r_beam.pGetNode(1002)));

    FindEmbeddedBeamHostElementsProcess process(
        model, Parameters{R"({"model_part_name": "Main.beam", "host_model_part_name": "Main.soil"})"});
    process.ExecuteInitialize();

    const auto host_ids = GetHostIds(r_main.GetElement(1001));
    ASSERT_EQ(host_ids.size(), 2u);
    array_1d<double, 3> local_coordinates;
    KRATOS_EXPECT_TRUE(
        r_main.GetElement(host_ids[0]).GetGeometry().IsInside(r_main.GetNode(1002).Coordinates(), local_coordinates, 1.0e-6) ||
        r_main.GetElement(host_ids[1]).GetGeometry().IsInside(r_main.GetNode(1002).Coordinates(), local_coordinates, 1.0e-6))
    KRATOS_EXPECT_TRUE(
        r_main.GetElement(host_ids[0]).GetGeometry().IsInside(r_main.GetNode(1001).Coordinates(), local_coordinates, 1.0e-6) ||
        r_main.GetElement(host_ids[1]).GetGeometry().IsInside(r_main.GetNode(1001).Coordinates(), local_coordinates, 1.0e-6))
}

KRATOS_TEST_CASE_IN_SUITE(FindEmbeddedBeamHostElementsProcess_ThrowsWhenHostModelPartIsNotSpecified,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    CreateMainModelPart(model);

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        FindEmbeddedBeamHostElementsProcess(model, Parameters{R"({"model_part_name": "Main.beam"})"}),
        "Please specify 'host_model_part_name' or 'host_model_part_name_list'")
}

KRATOS_TEST_CASE_IN_SUITE(FindEmbeddedBeamHostElementsProcess_ThrowsWhenThereAreNoSolidHostElements,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_main = CreateMainModelPart(model);
    auto& r_soil = r_main.GetSubModelPart("soil");
    r_soil.CreateNewNode(1, 0.0, 0.0, 0.0);
    r_soil.CreateNewNode(2, 1.0, 0.0, 0.0);
    AddElement(r_soil, 1, Kratos::make_shared<Line2D2<Node>>(r_soil.pGetNode(1), r_soil.pGetNode(2)));

    FindEmbeddedBeamHostElementsProcess process(
        model, Parameters{R"({"model_part_name": "Main.beam", "host_model_part_name": "Main.soil"})"});
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(process.ExecuteInitialize(),
                                      "No solid elements were found in the host model part(s)")
}

} // namespace Kratos::Testing

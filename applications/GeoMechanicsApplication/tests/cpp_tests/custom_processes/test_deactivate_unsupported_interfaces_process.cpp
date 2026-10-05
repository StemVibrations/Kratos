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
#include "custom_processes/deactivate_unsupported_interfaces_process.h"
#include "includes/kratos_flags.h"
#include "test_setup_utilities/element_setup_utilities.hpp"
#include "test_setup_utilities/model_setup_utilities.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

namespace
{

using namespace Kratos;
using namespace Kratos::Testing;

Element::Pointer AddElement(ModelPart& rModelPart, Element::Pointer pElement)
{
    pElement->SetId(rModelPart.NumberOfElements() + 1);
    rModelPart.AddElement(pElement);
    return pElement;
}

// Creates a 3+3 noded line interface (nodes 1-3 at the first side, nodes 4-6 at the second side),
// of which the first side is attached to a continuum element
struct InterfaceSetup {
    Element::Pointer pInterface;
    Element::Pointer pContinuumAtFirstSide;
};

InterfaceSetup CreateInterfaceWithContinuumAtFirstSide(ModelPart& rModelPart)
{
    ModelSetupUtilities::CreateNumberOfNewNodes(rModelPart, 16);
    auto result                  = InterfaceSetup{};
    result.pContinuumAtFirstSide = AddElement(
        rModelPart, ElementSetupUtilities::Create2D6NElement(
                        ModelSetupUtilities::GetNodesFromIds(rModelPart, {1, 2, 7, 3, 8, 9}), {}));
    result.pInterface = AddElement(
        rModelPart, ElementSetupUtilities::Create2D6NInterfaceElement(
                        ModelSetupUtilities::GetNodesFromIds(rModelPart, {1, 2, 3, 4, 5, 6}), {}));
    return result;
}

std::size_t DeactivateUnsupportedInterfaces(ModelPart& rModelPart)
{
    auto process = DeactivateUnsupportedInterfacesProcess{rModelPart};
    process.Execute();
    return process.GetNumberOfDeactivatedInterfaces();
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_KeepsInterfaceBetweenActiveContinuumElements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);
    AddElement(r_model_part, ElementSetupUtilities::Create2D6NElement(
                                 ModelSetupUtilities::GetNodesFromIds(r_model_part, {4, 5, 10, 6, 11, 12}), {}));

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 0);
    KRATOS_EXPECT_TRUE(setup.pInterface->IsActive())
}

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_DeactivatesInterfaceWhenContinuumAtOneSideIsInactive,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: the continuum at the second side is excavated
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);
    auto p_excavated_element = AddElement(
        r_model_part, ElementSetupUtilities::Create2D6NElement(
                          ModelSetupUtilities::GetNodesFromIds(r_model_part, {4, 5, 10, 6, 11, 12}), {}));
    p_excavated_element->Set(ACTIVE, false);

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 1);
    KRATOS_EXPECT_FALSE(setup.pInterface->IsActive())
    KRATOS_EXPECT_TRUE(setup.pContinuumAtFirstSide->IsActive())
}

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_LineElementCanSupportASideOfALineInterface,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: a (structural) line element that contains both corner nodes of the second side
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);
    AddElement(r_model_part, ElementSetupUtilities::Create2D2NElement(
                                 ModelSetupUtilities::GetNodesFromIds(r_model_part, {4, 5}), {}));

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 0);
    KRATOS_EXPECT_TRUE(setup.pInterface->IsActive())
}

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_DeactivatesInterfaceOfWhichASideIsOnlyPartiallySupported,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: the continuum at the second side is excavated, and the only active element at that
    // side contains just one of its corner nodes
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);
    auto p_excavated_element = AddElement(
        r_model_part, ElementSetupUtilities::Create2D6NElement(
                          ModelSetupUtilities::GetNodesFromIds(r_model_part, {4, 5, 10, 6, 11, 12}), {}));
    p_excavated_element->Set(ACTIVE, false);
    AddElement(r_model_part, ElementSetupUtilities::Create2D2NElement(
                                 ModelSetupUtilities::GetNodesFromIds(r_model_part, {4, 13}), {}));

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 1);
    KRATOS_EXPECT_FALSE(setup.pInterface->IsActive())
}

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_KeepsInterfaceOfWhichASideNeverHadAnAdjacentElement,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange: the second side is not attached to any element (e.g. it is supported by boundary conditions)
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 0);
    KRATOS_EXPECT_TRUE(setup.pInterface->IsActive())
}

KRATOS_TEST_CASE_IN_SUITE(DeactivateUnsupportedInterfacesProcess_IgnoresInactiveInterfaces, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto setup   = CreateInterfaceWithContinuumAtFirstSide(r_model_part);
    setup.pInterface->Set(ACTIVE, false);

    // Act & Assert
    KRATOS_EXPECT_EQ(DeactivateUnsupportedInterfaces(r_model_part), 0);
    KRATOS_EXPECT_FALSE(setup.pInterface->IsActive())
}

} // namespace Kratos::Testing

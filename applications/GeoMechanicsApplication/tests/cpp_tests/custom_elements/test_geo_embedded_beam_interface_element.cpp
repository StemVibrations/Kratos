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

#include "custom_elements/geo_embedded_beam_interface_element.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

#include <numbers>

namespace
{

using namespace Kratos;

class StubHostElementWithStress : public Element
{
public:
    StubHostElementWithStress(IndexType NewId, GeometryType::Pointer pGeometry, const Vector& rStressVector)
        : Element(NewId, pGeometry), mStressVector(rStressVector)
    {
    }

    using Element::CalculateOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&) override
    {
        if (rVariable == CAUCHY_STRESS_VECTOR) rOutput.assign(2, mStressVector);
    }

private:
    Vector mStressVector;
};

ModelPart& CreateModelPartWithDisplacements(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.AddNodalSolutionStepVariable(VELOCITY);
    r_result.AddNodalSolutionStepVariable(ACCELERATION);
    r_result.AddNodalSolutionStepVariable(REACTION);
    return r_result;
}

void AddDisplacementDofs(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
    }
}

Properties::Pointer CreateInterfaceProperties(double ShearStiffness, double NormalStiffness)
{
    auto result = std::make_shared<Properties>(1);
    result->SetValue(INTERFACE_SHEAR_STIFFNESS, ShearStiffness);
    result->SetValue(INTERFACE_NORMAL_STIFFNESS, NormalStiffness);
    return result;
}

Geometry<Node>::Pointer CreateHexahedronGeometry(ModelPart& rModelPart)
{
    return Kratos::make_shared<Hexahedra3D8<Node>>(
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0), rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0),
        rModelPart.CreateNewNode(3, 2.0, 2.0, 0.0), rModelPart.CreateNewNode(4, 0.0, 2.0, 0.0),
        rModelPart.CreateNewNode(5, 0.0, 0.0, 2.0), rModelPart.CreateNewNode(6, 2.0, 0.0, 2.0),
        rModelPart.CreateNewNode(7, 2.0, 2.0, 2.0), rModelPart.CreateNewNode(8, 0.0, 2.0, 2.0));
}

Geometry<Node>::Pointer CreateQuadrilateralGeometry(ModelPart& rModelPart)
{
    return Kratos::make_shared<Quadrilateral2D4<Node>>(
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0), rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0),
        rModelPart.CreateNewNode(3, 2.0, 2.0, 0.0), rModelPart.CreateNewNode(4, 0.0, 2.0, 0.0));
}

void AssignHost(Element& rElement, const Element::Pointer& rpHost)
{
    GlobalPointersVector<Element> hosts;
    hosts.push_back(GlobalPointer<Element>(rpHost.get()));
    rElement.SetValue(NEIGHBOUR_ELEMENTS, hosts);
}

array_1d<double, 3> MakeArray(double X, double Y, double Z)
{
    array_1d<double, 3> result;
    result[0] = X;
    result[1] = Y;
    result[2] = Z;
    return result;
}

/// A vertical beam (nodes 9 and 10, length 1.6) inside a 2 x 2 x 2 hexahedron (nodes 1-8)
struct VerticalBeamIn3DHost {
    explicit VerticalBeamIn3DHost(const Properties::Pointer& rpProperties,
                                  const Vector&              rHostStress   = Vector{},
                                  double                     ZOfSecondNode = 1.8)
        : r_model_part(CreateModelPartWithDisplacements(model))
    {
        auto p_host_geometry = CreateHexahedronGeometry(r_model_part);
        p_host               = rHostStress.size() > 0
                                   ? Element::Pointer(Kratos::make_intrusive<StubHostElementWithStress>(
                           100, p_host_geometry, rHostStress))
                                   : Kratos::make_intrusive<Element>(100, p_host_geometry);
        auto p_geometry =
            Kratos::make_shared<Line3D2<Node>>(r_model_part.CreateNewNode(9, 1.0, 1.0, 0.2),
                                               r_model_part.CreateNewNode(10, 1.0, 1.0, ZOfSecondNode));
        AddDisplacementDofs(r_model_part);
        p_element = Kratos::make_intrusive<GeoEmbeddedBeamInterfaceElement<3>>(1, p_geometry, rpProperties);
        AssignHost(*p_element, p_host);
        p_element->Initialize(process_info);
    }

    Model                                       model;
    ModelPart&                                  r_model_part;
    Element::Pointer                            p_host;
    GeoEmbeddedBeamInterfaceElement<3>::Pointer p_element;
    ProcessInfo                                 process_info;
};

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement_ThrowsWhenNoHostElementsAreAssigned, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithDisplacements(model);
    auto  p_geometry   = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0), r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    GeoEmbeddedBeamInterfaceElement<2> element(1, p_geometry, CreateInterfaceProperties(1.0, 1.0));

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        element.Initialize(ProcessInfo{}),
        "No candidate host elements have been assigned to embedded beam interface element 1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement3D2N_ElasticForcesAreProportionalToRelativeDisplacements,
                          KratosGeoMechanicsFastSuite)
{
    const auto shear_stiffness  = 1000.0;
    const auto normal_stiffness = 5000.0;
    auto setup = VerticalBeamIn3DHost(CreateInterfaceProperties(shear_stiffness, normal_stiffness));

    // Dofs: beam nodes 9 and 10, followed by host nodes 1-8 (3 dofs each)
    Element::EquationIdVectorType equation_ids;
    Element::DofsVectorType       dofs;
    setup.p_element->GetDofList(dofs, setup.process_info);
    ASSERT_EQ(dofs.size(), 30u);
    KRATOS_EXPECT_EQ(dofs[0]->Id(), 9u);
    KRATOS_EXPECT_EQ(dofs[3]->Id(), 10u);
    KRATOS_EXPECT_EQ(dofs[6]->Id(), 1u);

    // A tangential (axial) slip of both beam nodes and a normal displacement of beam node 9
    const auto slip = 1.0e-3;
    setup.r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(2.0e-3, 0.0, slip);
    setup.r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, slip);

    Matrix lhs;
    Vector rhs;
    setup.p_element->CalculateLocalSystem(lhs, rhs, setup.process_info);

    const auto half_length = 0.8;
    KRATOS_EXPECT_NEAR(rhs[0], -normal_stiffness * 2.0e-3 * half_length, 1.0e-12);
    KRATOS_EXPECT_NEAR(rhs[1], 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(rhs[2], -shear_stiffness * slip * half_length, 1.0e-12);
    KRATOS_EXPECT_NEAR(rhs[5], -shear_stiffness * slip * half_length, 1.0e-12);

    // The host nodes take the opposite forces
    double total_force_z = 0.0;
    double total_force_x = 0.0;
    for (std::size_t i = 0; i < 10; ++i) {
        total_force_x += rhs[3 * i];
        total_force_z += rhs[3 * i + 2];
    }
    KRATOS_EXPECT_NEAR(total_force_x, 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(total_force_z, 0.0, 1.0e-12);

    // The stiffness of a beam node
    KRATOS_EXPECT_NEAR(lhs(0, 0), normal_stiffness * half_length, 1.0e-9);
    KRATOS_EXPECT_NEAR(lhs(1, 1), normal_stiffness * half_length, 1.0e-9);
    KRATOS_EXPECT_NEAR(lhs(2, 2), shear_stiffness * half_length, 1.0e-9);

    // The local relative displacements [tangential, normal 1, normal 2] at the beam nodes. The
    // local axes of the vertical beam are x' = z, y' = y and z' = -x.
    std::vector<Vector> relative_displacements;
    setup.p_element->CalculateOnIntegrationPoints(GEO_RELATIVE_DISPLACEMENT_VECTOR,
                                                  relative_displacements, setup.process_info);
    ASSERT_EQ(relative_displacements.size(), 2u);
    KRATOS_EXPECT_VECTOR_NEAR(relative_displacements[0], MakeArray(slip, 0.0, -2.0e-3), 1.0e-12)
    KRATOS_EXPECT_VECTOR_NEAR(relative_displacements[1], MakeArray(slip, 0.0, 0.0), 1.0e-12)
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement3D2N_ShearForceIsLimitedByCohesion, KratosGeoMechanicsFastSuite)
{
    const auto shear_stiffness  = 1.0e5;
    const auto normal_stiffness = 5000.0;
    auto       p_properties     = CreateInterfaceProperties(shear_stiffness, normal_stiffness);
    const auto cohesion         = 10.0;
    const auto diameter         = 0.2;
    p_properties->SetValue(GEO_COHESION, cohesion);
    p_properties->SetValue(DIAMETER, diameter);
    auto setup = VerticalBeamIn3DHost(p_properties);

    setup.r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, 1.0e-3);
    setup.r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, 1.0e-3);

    Matrix lhs;
    Vector rhs;
    setup.p_element->CalculateLocalSystem(lhs, rhs, setup.process_info);

    const auto maximum_shear_force = std::numbers::pi * diameter * cohesion;
    KRATOS_EXPECT_NEAR(rhs[2], -maximum_shear_force * 0.8, 1.0e-12);
    KRATOS_EXPECT_NEAR(rhs[5], -maximum_shear_force * 0.8, 1.0e-12);

    // While sliding, there is no tangential stiffness, while the normal stiffness is not affected
    KRATOS_EXPECT_NEAR(lhs(2, 2), 0.0, 1.0e-12);
    KRATOS_EXPECT_NEAR(lhs(0, 0), normal_stiffness * 0.8, 1.0e-9);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement3D2N_FrictionDependsOnNormalStressOfHostElement,
                          KratosGeoMechanicsFastSuite)
{
    auto       p_properties   = CreateInterfaceProperties(1.0e5, 5000.0);
    const auto cohesion       = 5.0;
    const auto friction_angle = 30.0;
    const auto diameter       = 0.2;
    p_properties->SetValue(GEO_COHESION, cohesion);
    p_properties->SetValue(GEO_FRICTION_ANGLE, friction_angle);
    p_properties->SetValue(DIAMETER, diameter);

    // [xx, yy, zz, xy, yz, xz]. The normal stress on the shaft of a vertical beam is the average of xx and yy.
    Vector host_stress = ZeroVector(6);
    host_stress[0]     = -100.0;
    host_stress[1]     = -60.0;
    host_stress[2]     = -300.0;
    auto setup         = VerticalBeamIn3DHost(p_properties, host_stress);
    setup.p_element->InitializeSolutionStep(setup.process_info);

    setup.r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, -1.0e-2);
    setup.r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, -1.0e-2);

    Vector rhs;
    setup.p_element->CalculateRightHandSide(rhs, setup.process_info);

    const auto maximum_shear_force =
        std::numbers::pi * diameter * (cohesion + 80.0 * std::tan(friction_angle * std::numbers::pi / 180.0));
    KRATOS_EXPECT_NEAR(rhs[2], maximum_shear_force * 0.8, 1.0e-9);
    KRATOS_EXPECT_NEAR(rhs[5], maximum_shear_force * 0.8, 1.0e-9);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement3D2N_NodesOutsideHostDoNotContribute, KratosGeoMechanicsFastSuite)
{
    // The second beam node is located above the host
    auto setup = VerticalBeamIn3DHost(CreateInterfaceProperties(1000.0, 5000.0), Vector{}, 2.5);
    KRATOS_EXPECT_TRUE(setup.p_element->GetHostInterpolation().IsEmbedded(0))
    KRATOS_EXPECT_FALSE(setup.p_element->GetHostInterpolation().IsEmbedded(1))

    setup.r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, 0.0, 1.0e-3);
    setup.r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(1.0e-3, 0.0, 1.0e-3);

    Vector rhs;
    setup.p_element->CalculateRightHandSide(rhs, setup.process_info);

    const auto half_length = 0.5 * 2.3;
    KRATOS_EXPECT_NEAR(rhs[2], -1000.0 * 1.0e-3 * half_length, 1.0e-12);
    for (std::size_t i = 3; i < 6; ++i) {
        KRATOS_EXPECT_NEAR(rhs[i], 0.0, 1.0e-12);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement3D2N_DoesNotThrowForUnsupportedOutputVariables,
                          KratosGeoMechanicsFastSuite)
{
    auto setup = VerticalBeamIn3DHost(CreateInterfaceProperties(1000.0, 5000.0));

    // E.g. output processes request stresses from all elements of a model part
    std::vector<Vector> stress_vectors;
    EXPECT_NO_THROW(setup.p_element->CalculateOnIntegrationPoints(
        CAUCHY_STRESS_VECTOR, stress_vectors, setup.process_info));
    KRATOS_EXPECT_TRUE(stress_vectors.empty())
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement2D2N_SlipIsRetainedAfterFinalizingTheSolutionStep,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithDisplacements(model);
    auto  p_host = Kratos::make_intrusive<Element>(100, CreateQuadrilateralGeometry(r_model_part));
    auto  p_geometry = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(9, 1.0, 0.2, 0.0), r_model_part.CreateNewNode(10, 1.0, 1.8, 0.0));
    AddDisplacementDofs(r_model_part);

    const auto shear_stiffness = 1.0e4;
    const auto cohesion        = 10.0;
    auto       p_properties    = CreateInterfaceProperties(shear_stiffness, 5000.0);
    p_properties->SetValue(GEO_COHESION, cohesion);

    GeoEmbeddedBeamInterfaceElement<2> element(1, p_geometry, p_properties);
    AssignHost(element, p_host);
    const auto process_info = ProcessInfo{};
    element.Initialize(process_info);
    element.InitializeSolutionStep(process_info);

    // In 2D, both faces of the beam are in contact with the soil: the maximum shear force per unit
    // length is 2 * cohesion
    const auto maximum_shear_force                                  = 2.0 * cohesion;
    const auto slip                                                 = 1.0e-2;
    r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT)  = MakeArray(0.0, slip, 0.0);
    r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, slip, 0.0);

    std::vector<Vector> tractions;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions, process_info);
    KRATOS_EXPECT_NEAR(tractions[0][0], maximum_shear_force, 1.0e-9);

    element.FinalizeSolutionStep(process_info);
    element.InitializeSolutionStep(process_info);

    // Partial unloading is elastic, starting from the plastic slip reached in the previous step
    const auto unloading = 1.0e-3;
    r_model_part.GetNode(9).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, slip - unloading, 0.0);
    r_model_part.GetNode(10).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, slip - unloading, 0.0);

    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions, process_info);
    const auto plastic_slip = slip - maximum_shear_force / shear_stiffness;
    KRATOS_EXPECT_NEAR(tractions[0][0], shear_stiffness * (slip - unloading - plastic_slip), 1.0e-9);
    KRATOS_EXPECT_NEAR(tractions[1][0], shear_stiffness * (slip - unloading - plastic_slip), 1.0e-9);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement2D2N_BondingForcesAreRetainedInTheNextStage,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithDisplacements(model);
    auto  p_host = Kratos::make_intrusive<Element>(100, CreateQuadrilateralGeometry(r_model_part));
    auto  p_geometry = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(9, 1.0, 0.2, 0.0), r_model_part.CreateNewNode(10, 1.0, 1.8, 0.0));
    AddDisplacementDofs(r_model_part);

    const auto shear_stiffness = 1.0e4;
    const auto cohesion        = 10.0;
    auto       p_properties    = CreateInterfaceProperties(shear_stiffness, 5000.0);
    p_properties->SetValue(GEO_COHESION, cohesion);

    GeoEmbeddedBeamInterfaceElement<2> element(1, p_geometry, p_properties);
    AssignHost(element, p_host);
    const auto process_info = ProcessInfo{};

    // Stage 1: the beam slips along its axis and is displaced normal to its axis
    element.Initialize(process_info);
    element.InitializeSolutionStep(process_info);
    const auto maximum_shear_force = 2.0 * cohesion;
    const auto slip                = 1.0e-2;
    const auto normal_displacement = 1.0e-3;
    for (const auto node_id : {9, 10}) {
        r_model_part.GetNode(node_id).FastGetSolutionStepValue(DISPLACEMENT) =
            MakeArray(normal_displacement, slip, 0.0);
    }

    std::vector<Vector> tractions_stage_1;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions_stage_1, process_info);
    Vector rhs_stage_1;
    element.CalculateRightHandSide(rhs_stage_1, process_info);
    element.FinalizeSolutionStep(process_info);

    // Stage 2: the (stage) displacements are reset and the element is initialized again
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
    }
    element.Initialize(process_info);
    element.InitializeSolutionStep(process_info);

    // Without any displacement in this stage, the bonding forces of stage 1 are retained
    std::vector<Vector> tractions;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions, process_info);
    ASSERT_EQ(tractions.size(), 2u);
    KRATOS_EXPECT_NEAR(tractions[0][0], maximum_shear_force, 1.0e-9);
    KRATOS_EXPECT_VECTOR_NEAR(tractions[0], tractions_stage_1[0], 1.0e-9)
    KRATOS_EXPECT_VECTOR_NEAR(tractions[1], tractions_stage_1[1], 1.0e-9)
    Vector rhs;
    element.CalculateRightHandSide(rhs, process_info);
    KRATOS_EXPECT_VECTOR_NEAR(rhs, rhs_stage_1, 1.0e-9)

    // Partial unloading in stage 2 is elastic, starting from the bonding forces of stage 1
    const auto unloading = 1.0e-3;
    for (const auto node_id : {9, 10}) {
        r_model_part.GetNode(node_id).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(0.0, -unloading, 0.0);
    }
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions, process_info);
    KRATOS_EXPECT_NEAR(tractions[0][0], maximum_shear_force - shear_stiffness * unloading, 1.0e-9);
    KRATOS_EXPECT_NEAR(tractions[1][0], maximum_shear_force - shear_stiffness * unloading, 1.0e-9);
    KRATOS_EXPECT_NEAR(tractions[0][1], tractions_stage_1[0][1], 1.0e-9);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamInterfaceElement2D2N_InactiveElementLosesItsBondingForces,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithDisplacements(model);
    auto  p_host = Kratos::make_intrusive<Element>(100, CreateQuadrilateralGeometry(r_model_part));
    auto  p_geometry = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(9, 1.0, 0.2, 0.0), r_model_part.CreateNewNode(10, 1.0, 1.8, 0.0));
    AddDisplacementDofs(r_model_part);

    GeoEmbeddedBeamInterfaceElement<2> element(1, p_geometry, CreateInterfaceProperties(1.0e4, 5000.0));
    AssignHost(element, p_host);
    const auto process_info = ProcessInfo{};

    element.Initialize(process_info);
    element.InitializeSolutionStep(process_info);
    for (const auto node_id : {9, 10}) {
        r_model_part.GetNode(node_id).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(1.0e-3, 1.0e-3, 0.0);
    }
    element.FinalizeSolutionStep(process_info);

    // The element is deactivated in the next stage
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
    }
    element.Set(ACTIVE, false);
    element.Initialize(process_info);

    // And reactivated in the stage after that
    element.Set(ACTIVE, true);
    element.Initialize(process_info);

    std::vector<Vector> tractions;
    element.CalculateOnIntegrationPoints(GEO_EFFECTIVE_TRACTION_VECTOR, tractions, process_info);
    ASSERT_EQ(tractions.size(), 2u);
    KRATOS_EXPECT_VECTOR_NEAR(tractions[0], ZeroVector(2), 1.0e-12)
    KRATOS_EXPECT_VECTOR_NEAR(tractions[1], ZeroVector(2), 1.0e-12)
}

} // namespace Kratos::Testing

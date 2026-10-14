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

#include "custom_constitutive/timoshenko_beam_elastic_constitutive_law.h"
#include "custom_constitutive/timoshenko_beam_elastic_constitutive_law_3d.h"
#include "custom_elements/beam_elements/timoshenko_beam_element_2D2N.h"
#include "custom_elements/geo_embedded_beam_element.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_2d_4.h"
#include "structural_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "utilities/math_utils.h"

namespace
{

using namespace Kratos;

ModelPart& CreateModelPartWithBeamVariables(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(DISPLACEMENT);
    r_result.AddNodalSolutionStepVariable(ROTATION);
    r_result.AddNodalSolutionStepVariable(VELOCITY);
    r_result.AddNodalSolutionStepVariable(ACCELERATION);
    r_result.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);
    r_result.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
    r_result.AddNodalSolutionStepVariable(REACTION);
    r_result.AddNodalSolutionStepVariable(REACTION_MOMENT);
    return r_result;
}

void AddDisplacementAndRotationDofs(ModelPart& rModelPart)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.AddDof(DISPLACEMENT_X, REACTION_X);
        r_node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        r_node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        r_node.AddDof(ROTATION_X, REACTION_MOMENT_X);
        r_node.AddDof(ROTATION_Y, REACTION_MOMENT_Y);
        r_node.AddDof(ROTATION_Z, REACTION_MOMENT_Z);
    }
}

Properties::Pointer CreateBeamProperties2D()
{
    auto result = std::make_shared<Properties>(1);
    result->SetValue(YOUNG_MODULUS, 2.0e8);
    result->SetValue(POISSON_RATIO, 0.3);
    result->SetValue(CROSS_AREA, 0.01);
    result->SetValue(I33, 2.0e-5);
    result->SetValue(AREA_EFFECTIVE_Y, 0.01 * 5.0 / 6.0);
    result->SetValue(DENSITY, 7.85);
    result->SetValue(CONSTITUTIVE_LAW, std::make_shared<TimoshenkoBeamElasticConstitutiveLaw>());
    return result;
}

Properties::Pointer CreateBeamProperties3D(double ShearArea = 0.01 * 5.0 / 6.0)
{
    auto result = std::make_shared<Properties>(1);
    result->SetValue(YOUNG_MODULUS, 2.0e8);
    result->SetValue(POISSON_RATIO, 0.3);
    result->SetValue(CROSS_AREA, 0.01);
    result->SetValue(I22, 3.0e-5);
    result->SetValue(I33, 2.0e-5);
    result->SetValue(AREA_EFFECTIVE_Y, ShearArea);
    result->SetValue(AREA_EFFECTIVE_Z, ShearArea);
    result->SetValue(DENSITY, 7.85);
    result->SetValue(CONSTITUTIVE_LAW, std::make_shared<TimoshenkoBeamElasticConstitutiveLaw3D>());
    return result;
}

Element::Pointer CreateHostQuadrilateral(ModelPart& rModelPart)
{
    // Nodes 1-4: a square of 2 x 2
    auto p_geometry = Kratos::make_shared<Quadrilateral2D4<Node>>(
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0), rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0),
        rModelPart.CreateNewNode(3, 2.0, 2.0, 0.0), rModelPart.CreateNewNode(4, 0.0, 2.0, 0.0));
    return Kratos::make_intrusive<Element>(100, p_geometry);
}

Element::Pointer CreateHostHexahedron(ModelPart& rModelPart)
{
    // Nodes 1-8: a cube of 2 x 2 x 2
    auto p_geometry = Kratos::make_shared<Hexahedra3D8<Node>>(
        rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0), rModelPart.CreateNewNode(2, 2.0, 0.0, 0.0),
        rModelPart.CreateNewNode(3, 2.0, 2.0, 0.0), rModelPart.CreateNewNode(4, 0.0, 2.0, 0.0),
        rModelPart.CreateNewNode(5, 0.0, 0.0, 2.0), rModelPart.CreateNewNode(6, 2.0, 0.0, 2.0),
        rModelPart.CreateNewNode(7, 2.0, 2.0, 2.0), rModelPart.CreateNewNode(8, 0.0, 2.0, 2.0));
    return Kratos::make_intrusive<Element>(100, p_geometry);
}

void AssignHosts(Element& rElement, const std::vector<Element::Pointer>& rHosts)
{
    GlobalPointersVector<Element> hosts;
    for (const auto& rp_host : rHosts) {
        hosts.push_back(GlobalPointer<Element>(rp_host.get()));
    }
    rElement.SetValue(NEIGHBOUR_ELEMENTS, hosts);
}

void ExpectMatricesNear(const Matrix& rActual, const Matrix& rExpected, double RelativeTolerance = 1.0e-9)
{
    const double tolerance = RelativeTolerance * std::max<double>(norm_frobenius(rExpected), 1.0);
    KRATOS_EXPECT_MATRIX_NEAR(rActual, rExpected, tolerance)
}

void ExpectVectorsNear(const Vector& rActual, const Vector& rExpected, double AbsoluteTolerance){
    KRATOS_EXPECT_VECTOR_NEAR(rActual, rExpected, AbsoluteTolerance)}

array_1d<double, 3> MakeArray(double X, double Y, double Z)
{
    array_1d<double, 3> result;
    result[0] = X;
    result[1] = Y;
    result[2] = Z;
    return result;
}

/// Imposes the rigid body motion u(x) = a + omega x (x - x_0), theta = omega
void ImposeRigidBodyMotion(ModelPart&                 rModelPart,
                           const array_1d<double, 3>& rTranslation,
                           const array_1d<double, 3>& rRotation)
{
    for (auto& r_node : rModelPart.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) =
            rTranslation +
            MathUtils<double>::CrossProduct(rRotation, r_node.GetInitialPosition().Coordinates());
        r_node.FastGetSolutionStepValue(ROTATION) = rRotation;
    }
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement_ThrowsWhenNoHostElementsAreAssigned, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithBeamVariables(model);
    auto  p_geometry   = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(1, 0.0, 0.0, 0.0), r_model_part.CreateNewNode(2, 1.0, 0.0, 0.0));
    GeoEmbeddedBeamElement<2> element(1, p_geometry, CreateBeamProperties2D());

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        element.Initialize(ProcessInfo{}),
        "No candidate host elements have been assigned to embedded beam element 1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement2D2N_DofsConsistOfHostTranslationsAndBeamRotations, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithBeamVariables(model);
    auto  p_host       = CreateHostQuadrilateral(r_model_part);
    auto  p_geometry   = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(5, 0.5, 0.5, 0.0), r_model_part.CreateNewNode(6, 1.5, 1.0, 0.0));
    AddDisplacementAndRotationDofs(r_model_part);

    GeoEmbeddedBeamElement<2> element(1, p_geometry, CreateBeamProperties2D());
    AssignHosts(element, {p_host});
    const auto process_info = ProcessInfo{};
    element.Initialize(process_info);

    Element::DofsVectorType dofs;
    element.GetDofList(dofs, process_info);
    ASSERT_EQ(dofs.size(), 10u);
    const std::vector<std::pair<std::size_t, std::string>> expected_dofs = {
        {1, "DISPLACEMENT_X"}, {1, "DISPLACEMENT_Y"}, {2, "DISPLACEMENT_X"}, {2, "DISPLACEMENT_Y"},
        {3, "DISPLACEMENT_X"}, {3, "DISPLACEMENT_Y"}, {4, "DISPLACEMENT_X"}, {4, "DISPLACEMENT_Y"},
        {5, "ROTATION_Z"},     {6, "ROTATION_Z"}};
    for (std::size_t i = 0; i < dofs.size(); ++i) {
        KRATOS_EXPECT_EQ(dofs[i]->Id(), expected_dofs[i].first);
        KRATOS_EXPECT_EQ(dofs[i]->GetVariable().Name(), expected_dofs[i].second);
    }

    // The translations of the beam nodes are interpolated with the bilinear shape functions of the host
    const auto transformation_matrix = element.CalculateTransformationMatrix();
    ASSERT_EQ(transformation_matrix.size1(), 6u);
    ASSERT_EQ(transformation_matrix.size2(), 10u);
    Matrix     expected_transformation_matrix = ZeroMatrix(6, 10);
    const auto set_row = [&expected_transformation_matrix](std::size_t Row, std::size_t Direction,
                                                           const std::vector<double>& rWeights) {
        for (std::size_t k = 0; k < rWeights.size(); ++k) {
            expected_transformation_matrix(Row, 2 * k + Direction) = rWeights[k];
        }
    };
    set_row(0, 0, {0.5625, 0.1875, 0.0625, 0.1875});
    set_row(1, 1, {0.5625, 0.1875, 0.0625, 0.1875});
    set_row(3, 0, {0.125, 0.375, 0.375, 0.125});
    set_row(4, 1, {0.125, 0.375, 0.375, 0.125});
    expected_transformation_matrix(2, 8) = 1.0;
    expected_transformation_matrix(5, 9) = 1.0;
    ExpectMatricesNear(transformation_matrix, expected_transformation_matrix);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement2D2N_ReproducesTimoshenkoBeamWhenBeamNodesAreHostNodes,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithBeamVariables(model);
    auto  p_host       = CreateHostQuadrilateral(r_model_part);
    AddDisplacementAndRotationDofs(r_model_part);

    // An inclined beam along the diagonal of the host, connecting host nodes 1 and 3
    auto p_geometry =
        Kratos::make_shared<Line2D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(3));
    const auto p_properties = CreateBeamProperties2D();
    const auto process_info = ProcessInfo{};

    GeoEmbeddedBeamElement<2> embedded_beam(1, p_geometry, p_properties);
    AssignHosts(embedded_beam, {p_host});
    embedded_beam.Initialize(process_info);

    LinearTimoshenkoBeamElement2D2N reference_beam(2, p_geometry, p_properties);
    reference_beam.Initialize(process_info);

    r_model_part.GetNode(1).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(1.0e-3, -2.0e-3, 0.0);
    r_model_part.GetNode(1).FastGetSolutionStepValue(ROTATION) = MakeArray(0.0, 0.0, 3.0e-3);
    r_model_part.GetNode(3).FastGetSolutionStepValue(DISPLACEMENT) = MakeArray(-4.0e-3, 5.0e-4, 0.0);
    r_model_part.GetNode(3).FastGetSolutionStepValue(ROTATION) = MakeArray(0.0, 0.0, -1.0e-3);

    Matrix embedded_lhs;
    Vector embedded_rhs;
    embedded_beam.CalculateLocalSystem(embedded_lhs, embedded_rhs, process_info);
    Matrix reference_lhs;
    Vector reference_rhs;
    reference_beam.CalculateLocalSystem(reference_lhs, reference_rhs, process_info);

    // Dofs of the embedded beam: [u1, u2, u3, u4, theta1, theta3], those of the reference beam: [u1, theta1, u3, theta3]
    const std::vector<std::size_t> embedded_indices = {0, 1, 8, 4, 5, 9};
    Matrix                         mapped_lhs(6, 6);
    Vector                         mapped_rhs(6);
    for (std::size_t i = 0; i < 6; ++i) {
        mapped_rhs[i] = embedded_rhs[embedded_indices[i]];
        for (std::size_t j = 0; j < 6; ++j) {
            mapped_lhs(i, j) = embedded_lhs(embedded_indices[i], embedded_indices[j]);
        }
    }
    ExpectMatricesNear(mapped_lhs, reference_lhs);
    ExpectVectorsNear(mapped_rhs, reference_rhs, 1.0e-9 * norm_2(reference_rhs));

    // The nodes of the host that are not connected to the beam do not get any contribution
    for (const std::size_t unconnected_index : {2u, 3u, 6u, 7u}) {
        KRATOS_EXPECT_NEAR(embedded_rhs[unconnected_index], 0.0, 1.0e-9);
        for (std::size_t j = 0; j < embedded_lhs.size2(); ++j) {
            KRATOS_EXPECT_NEAR(embedded_lhs(unconnected_index, j), 0.0, 1.0e-9);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement3D2N_RigidBodyMotionsDoNotGenerateForces, KratosGeoMechanicsFastSuite)
{
    const std::vector<std::pair<array_1d<double, 3>, array_1d<double, 3>>> rigid_body_motions = {
        {MakeArray(1.0e-3, -2.0e-3, 3.0e-3), MakeArray(0.0, 0.0, 0.0)},
        {MakeArray(0.0, 0.0, 0.0), MakeArray(1.0e-3, 0.0, 0.0)},
        {MakeArray(0.0, 0.0, 0.0), MakeArray(0.0, 1.0e-3, 0.0)},
        {MakeArray(0.0, 0.0, 0.0), MakeArray(0.0, 0.0, 1.0e-3)},
        {MakeArray(2.0e-3, 1.0e-3, -1.0e-3), MakeArray(-1.0e-3, 2.0e-3, 0.5e-3)}};

    // Inclined beams with nodes inside the host, including a vertical beam and a beam along -y
    const std::vector<std::pair<array_1d<double, 3>, array_1d<double, 3>>> beam_end_points = {
        {MakeArray(0.3, 0.4, 0.2), MakeArray(1.5, 1.2, 1.7)},
        {MakeArray(1.0, 1.0, 0.1), MakeArray(1.0, 1.0, 1.9)},
        {MakeArray(0.6, 1.8, 0.5), MakeArray(0.6, 0.2, 0.5)}};

    for (const auto& [r_start, r_end] : beam_end_points) {
        for (const auto& [r_translation, r_rotation] : rigid_body_motions) {
            Model model;
            auto& r_model_part = CreateModelPartWithBeamVariables(model);
            auto  p_host       = CreateHostHexahedron(r_model_part);
            auto  p_geometry   = Kratos::make_shared<Line3D2<Node>>(
                r_model_part.CreateNewNode(9, r_start[0], r_start[1], r_start[2]),
                r_model_part.CreateNewNode(10, r_end[0], r_end[1], r_end[2]));
            AddDisplacementAndRotationDofs(r_model_part);

            GeoEmbeddedBeamElement<3> element(1, p_geometry, CreateBeamProperties3D());
            AssignHosts(element, {p_host});
            const auto process_info = ProcessInfo{};
            element.Initialize(process_info);

            ImposeRigidBodyMotion(r_model_part, r_translation, r_rotation);

            Matrix lhs;
            Vector rhs;
            element.CalculateLocalSystem(lhs, rhs, process_info);
            Vector dof_values;
            element.GetValuesVector(dof_values);

            const auto force_scale = norm_frobenius(lhs) * norm_2(dof_values);
            ExpectVectorsNear(rhs, ZeroVector(rhs.size()), 1.0e-9 * force_scale);
            ExpectVectorsNear(prod(lhs, dof_values), ZeroVector(rhs.size()), 1.0e-9 * force_scale);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement3D2N_AxialForceFollowsFromStrainAlongInclinedBeamAxis,
                          KratosGeoMechanicsFastSuite)
{
    Model      model;
    auto&      r_model_part = CreateModelPartWithBeamVariables(model);
    auto       p_host       = CreateHostHexahedron(r_model_part);
    const auto start        = MakeArray(0.3, 0.4, 0.2);
    const auto end          = MakeArray(1.5, 1.2, 1.7);
    auto       p_geometry =
        Kratos::make_shared<Line3D2<Node>>(r_model_part.CreateNewNode(9, start[0], start[1], start[2]),
                                           r_model_part.CreateNewNode(10, end[0], end[1], end[2]));
    AddDisplacementAndRotationDofs(r_model_part);

    const auto                p_properties = CreateBeamProperties3D();
    GeoEmbeddedBeamElement<3> element(1, p_geometry, p_properties);
    AssignHosts(element, {p_host});
    const auto process_info = ProcessInfo{};
    element.Initialize(process_info);

    // A displacement field with a uniform strain along the beam axis: u(x) = strain (d . x) d
    const auto axis   = array_1d<double, 3>((end - start) / norm_2(end - start));
    const auto strain = 1.0e-4;
    for (auto& r_node : r_model_part.Nodes()) {
        r_node.FastGetSolutionStepValue(DISPLACEMENT) =
            strain * inner_prod(axis, r_node.GetInitialPosition().Coordinates()) * axis;
    }

    std::vector<double> axial_forces;
    element.CalculateOnIntegrationPoints(AXIAL_FORCE, axial_forces, process_info);
    ASSERT_FALSE(axial_forces.empty());
    const auto expected_axial_force = (*p_properties)[YOUNG_MODULUS] * (*p_properties)[CROSS_AREA] * strain;
    for (const auto axial_force : axial_forces) {
        KRATOS_EXPECT_RELATIVE_NEAR(axial_force, expected_axial_force, 1.0e-9);
    }

    // Since the beam is straight and the field is linear, the beam is not bent
    std::vector<double> bending_moments;
    element.CalculateOnIntegrationPoints(BENDING_MOMENT_Y, bending_moments, process_info);
    for (const auto bending_moment : bending_moments) {
        KRATOS_EXPECT_NEAR(bending_moment, 0.0, 1.0e-9 * expected_axial_force);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement3D2N_EndMomentsOfVerticalBeamFollowFromBendingStiffnessAboutLocalAxes,
                          KratosGeoMechanicsFastSuite)
{
    // A vertical beam connecting host nodes 1 and 5. The local axes are x' = z, y' = y and z' = -x
    // (see EmbeddedBeamUtilities::CalculateLocalAxes). A rotation about the global x axis bends the
    // beam about its local z' axis (stiffness E I33), a rotation about the global y axis bends it
    // about its local y' axis (stiffness E I22). The shear deformation is negligible.
    const auto shear_area = 1.0e6;
    const std::vector<std::pair<array_1d<double, 3>, double>> rotations_and_second_moments_of_area = {
        {MakeArray(1.0e-3, 0.0, 0.0), (*CreateBeamProperties3D())[I33]},
        {MakeArray(0.0, 1.0e-3, 0.0), (*CreateBeamProperties3D())[I22]}};

    for (const auto& [r_rotation, second_moment_of_area] : rotations_and_second_moments_of_area) {
        Model model;
        auto& r_model_part = CreateModelPartWithBeamVariables(model);
        auto  p_host       = CreateHostHexahedron(r_model_part);
        AddDisplacementAndRotationDofs(r_model_part);
        auto p_geometry =
            Kratos::make_shared<Line3D2<Node>>(r_model_part.pGetNode(1), r_model_part.pGetNode(5));

        const auto                p_properties = CreateBeamProperties3D(shear_area);
        GeoEmbeddedBeamElement<3> element(1, p_geometry, p_properties);
        AssignHosts(element, {p_host});
        const auto process_info = ProcessInfo{};
        element.Initialize(process_info);

        r_model_part.GetNode(1).FastGetSolutionStepValue(ROTATION) = r_rotation;

        Vector rhs;
        element.CalculateRightHandSide(rhs, process_info);

        // Dofs: 8 x 3 host translations, followed by the rotations of nodes 1 and 5
        const auto length    = 2.0;
        const auto stiffness = (*p_properties)[YOUNG_MODULUS] * second_moment_of_area / length;
        const auto theta     = norm_2(r_rotation);
        Vector     expected_moment_node_1 = -4.0 * stiffness * r_rotation;
        Vector     expected_moment_node_5 = -2.0 * stiffness * r_rotation;
        Vector     moment_node_1(3);
        Vector     moment_node_5(3);
        for (std::size_t i = 0; i < 3; ++i) {
            moment_node_1[i] = rhs[24 + i];
            moment_node_5[i] = rhs[27 + i];
        }
        ExpectVectorsNear(moment_node_1, expected_moment_node_1, 1.0e-6 * 4.0 * stiffness * theta);
        ExpectVectorsNear(moment_node_5, expected_moment_node_5, 1.0e-6 * 4.0 * stiffness * theta);

        // The shear forces at the beam ends: |V| = 6 E I theta / L^2, perpendicular to the rotation axis
        const auto                expected_shear_force = 6.0 * stiffness * theta / length;
        const array_1d<double, 3> shear_direction =
            MathUtils<double>::CrossProduct(MakeArray(0.0, 0.0, 1.0), r_rotation) / theta;
        array_1d<double, 3> force_node_1;
        for (std::size_t i = 0; i < 3; ++i) {
            force_node_1[i] = rhs[i]; // node 1 is the first host node
        }
        KRATOS_EXPECT_NEAR(std::abs(inner_prod(force_node_1, shear_direction)),
                           expected_shear_force, 1.0e-6 * expected_shear_force);
    }
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement3D2N_MassAndBodyForcesAreTransferredToHostNodes, KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithBeamVariables(model);
    auto  p_host       = CreateHostHexahedron(r_model_part);
    auto  p_geometry   = Kratos::make_shared<Line3D2<Node>>(
        r_model_part.CreateNewNode(9, 0.3, 0.4, 0.2), r_model_part.CreateNewNode(10, 1.5, 1.2, 1.7));
    AddDisplacementAndRotationDofs(r_model_part);
    const auto p_properties = CreateBeamProperties3D();
    const auto gravity      = MakeArray(0.0, 0.0, -9.81);
    p_properties->SetValue(VOLUME_ACCELERATION, gravity);

    GeoEmbeddedBeamElement<3> element(1, p_geometry, p_properties);
    AssignHosts(element, {p_host});
    const auto process_info = ProcessInfo{};
    element.Initialize(process_info);

    const auto length =
        norm_2(p_geometry->GetPoint(1).Coordinates() - p_geometry->GetPoint(0).Coordinates());
    const auto beam_mass = (*p_properties)[DENSITY] * (*p_properties)[CROSS_AREA] * length;

    // A rigid translation of the host in any direction mobilizes the total mass of the beam
    Matrix mass_matrix;
    element.CalculateMassMatrix(mass_matrix, process_info);
    for (std::size_t direction = 0; direction < 3; ++direction) {
        Vector rigid_translation = ZeroVector(mass_matrix.size1());
        for (std::size_t node = 0; node < 8; ++node) {
            rigid_translation[3 * node + direction] = 1.0;
        }
        KRATOS_EXPECT_RELATIVE_NEAR(
            inner_prod(rigid_translation, prod(mass_matrix, rigid_translation)), beam_mass, 1.0e-12);
    }

    // The total weight of the beam is transferred to the host nodes
    Vector rhs;
    element.CalculateRightHandSide(rhs, process_info);
    array_1d<double, 3> total_force = ZeroVector(3);
    for (std::size_t node = 0; node < 8; ++node) {
        for (std::size_t direction = 0; direction < 3; ++direction) {
            total_force[direction] += rhs[3 * node + direction];
        }
    }
    ExpectVectorsNear(total_force, beam_mass * gravity, 1.0e-9 * beam_mass * 9.81);
}

KRATOS_TEST_CASE_IN_SUITE(GeoEmbeddedBeamElement2D2N_FinalizeSolutionStepInterpolatesBeamNodeDisplacements,
                          KratosGeoMechanicsFastSuite)
{
    Model model;
    auto& r_model_part = CreateModelPartWithBeamVariables(model);
    auto  p_host       = CreateHostQuadrilateral(r_model_part);
    auto  p_geometry   = Kratos::make_shared<Line2D2<Node>>(
        r_model_part.CreateNewNode(5, 0.5, 0.5, 0.0), r_model_part.CreateNewNode(6, 1.5, 1.0, 0.0));
    AddDisplacementAndRotationDofs(r_model_part);

    GeoEmbeddedBeamElement<2> element(1, p_geometry, CreateBeamProperties2D());
    AssignHosts(element, {p_host});
    const auto process_info = ProcessInfo{};
    element.Initialize(process_info);

    // A linear displacement field u = (0.01 x, -0.02 y)
    for (auto& r_node : r_model_part.Nodes()) {
        const auto& r_position = r_node.GetInitialPosition();
        r_node.FastGetSolutionStepValue(DISPLACEMENT) =
            MakeArray(0.01 * r_position.X(), -0.02 * r_position.Y(), 0.0);
        r_node.FastGetSolutionStepValue(VELOCITY) = MakeArray(0.1 * r_position.Y(), 0.0, 0.0);
    }
    r_model_part.GetNode(5).FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);
    r_model_part.GetNode(6).FastGetSolutionStepValue(DISPLACEMENT) = ZeroVector(3);

    element.FinalizeSolutionStep(process_info);

    ExpectVectorsNear(r_model_part.GetNode(5).FastGetSolutionStepValue(DISPLACEMENT),
                      MakeArray(0.005, -0.01, 0.0), 1.0e-12);
    ExpectVectorsNear(r_model_part.GetNode(6).FastGetSolutionStepValue(DISPLACEMENT),
                      MakeArray(0.015, -0.02, 0.0), 1.0e-12);
    ExpectVectorsNear(r_model_part.GetNode(6).FastGetSolutionStepValue(VELOCITY),
                      MakeArray(0.1, 0.0, 0.0), 1.0e-12);
}

} // namespace Kratos::Testing

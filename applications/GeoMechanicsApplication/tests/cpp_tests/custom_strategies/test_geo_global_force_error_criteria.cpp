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
#include "custom_constitutive/local_error_data_provider.h"
#include "custom_strategies/convergence_criteria/current_stiffness_parameter_calculator.h"
#include "custom_strategies/convergence_criteria/geo_global_force_error_criteria.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "geometries/triangle_2d_3.h"
#include "includes/element.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <limits>

using namespace Kratos;

namespace
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;
using CriteriaType    = GeoGlobalForceErrorCriteria<SparseSpaceType, LocalSpaceType>;

class StubLocalErrorLaw : public ConstitutiveLaw, public LocalErrorDataProvider
{
public:
    StubLocalErrorLaw(Geo::StressPointType Type, double TotalStrainEnergyIncrement, double ElasticStrainEnergyIncrement)
        : mType(Type)
    {
        mData.TotalStrainEnergyIncrement   = TotalStrainEnergyIncrement;
        mData.ElasticStrainEnergyIncrement = ElasticStrainEnergyIncrement;
    }

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override
    {
        return std::make_shared<StubLocalErrorLaw>(
            mType, mData.TotalStrainEnergyIncrement, mData.ElasticStrainEnergyIncrement);
    }

    [[nodiscard]] std::optional<Geo::StressPointType> GetStressPointType() const override
    {
        return mType;
    }

    [[nodiscard]] Geo::LocalErrorData CalculateLocalErrorData(Parameters&) override { return mData; }

    void SetStrainEnergyIncrements(double TotalStrainEnergyIncrement, double ElasticStrainEnergyIncrement)
    {
        mData.TotalStrainEnergyIncrement   = TotalStrainEnergyIncrement;
        mData.ElasticStrainEnergyIncrement = ElasticStrainEnergyIncrement;
    }

private:
    Geo::StressPointType mType;
    Geo::LocalErrorData  mData;
};

class ElementWithConstitutiveLaws : public Element
{
public:
    ElementWithConstitutiveLaws(IndexType                             NewId,
                                GeometryType::Pointer                 pGeometry,
                                PropertiesType::Pointer               pProperties,
                                std::vector<ConstitutiveLaw::Pointer> ConstitutiveLaws)
        : Element(NewId, std::move(pGeometry), std::move(pProperties)),
          mConstitutiveLaws(std::move(ConstitutiveLaws))
    {
    }

    using Element::CalculateOnIntegrationPoints;
    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>&    rOutput,
                                      const ProcessInfo&) override
    {
        if (rVariable == CONSTITUTIVE_LAW) rOutput = mConstitutiveLaws;
    }

private:
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws;
};

// Adds a right-angled triangular element with legs of the given length (i.e. an area of half the
// squared leg length)
Element& AddTriangularElement(ModelPart& rModelPart, double LegLength, std::vector<ConstitutiveLaw::Pointer> ConstitutiveLaws)
{
    const auto first_node_id = rModelPart.NumberOfNodes() + 1;
    auto       p_geometry    = Kratos::make_shared<Triangle2D3<Node>>(
        rModelPart.CreateNewNode(first_node_id, 0.0, 0.0, 0.0),
        rModelPart.CreateNewNode(first_node_id + 1, LegLength, 0.0, 0.0),
        rModelPart.CreateNewNode(first_node_id + 2, 0.0, LegLength, 0.0));
    auto p_properties = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0) : rModelPart.CreateNewProperties(0);
    auto p_element    = Kratos::make_intrusive<ElementWithConstitutiveLaws>(
        rModelPart.NumberOfElements() + 1, p_geometry, p_properties, std::move(ConstitutiveLaws));
    rModelPart.AddElement(p_element);
    return *p_element;
}

std::shared_ptr<StubLocalErrorLaw> MakeSoilPoint(double TotalStrainEnergyIncrement, double ElasticStrainEnergyIncrement)
{
    return std::make_shared<StubLocalErrorLaw>(Geo::StressPointType::Soil, TotalStrainEnergyIncrement,
                                               ElasticStrainEnergyIncrement);
}

// Creates a node with free displacement degrees of freedom in x and y direction (equation ids 0 and
// 1), a fixed displacement degree of freedom in z direction (equation id 2) and a free water pressure
// degree of freedom (equation id 3)
ModelPart::DofsArrayType CreateDofSet(ModelPart& rModelPart)
{
    rModelPart.AddNodalSolutionStepVariable(DISPLACEMENT);
    rModelPart.AddNodalSolutionStepVariable(WATER_PRESSURE);
    auto p_node = rModelPart.CreateNewNode(1000, 0.0, 0.0, 0.0);

    auto result      = ModelPart::DofsArrayType{};
    auto equation_id = std::size_t{0};
    for (const auto* p_variable : {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z, &WATER_PRESSURE}) {
        p_node->AddDof(*p_variable);
        auto p_dof = p_node->pGetDof(*p_variable);
        p_dof->SetEquationId(equation_id++);
        result.push_back(p_dof);
    }
    p_node->pGetDof(DISPLACEMENT_Z)->FixDof();

    return result;
}

bool PostCriteria(CriteriaType& rCriteria, ModelPart& rModelPart, ModelPart::DofsArrayType& rDofSet, const Vector& rResidual)
{
    auto A  = CriteriaType::TSystemMatrixType{};
    auto Dx = CriteriaType::TSystemVectorType{};
    return rCriteria.PostCriteria(rModelPart, rDofSet, A, Dx, rResidual);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_HasExpectedDefaultSettings, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto criteria = CriteriaType{};

    KRATOS_EXPECT_DOUBLE_EQ(criteria.GetForceRelativeTolerance(), 0.01);
    KRATOS_EXPECT_DOUBLE_EQ(criteria.GetForceAbsoluteTolerance(), 1.0e-9);
    KRATOS_EXPECT_FALSE(criteria.GetActualizeRHSflag())
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_ThrowsForInvalidTolerances, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CriteriaType{Parameters{R"({"force_relative_tolerance": 0.0})"}},
                                      "The force_relative_tolerance must be positive, got 0")
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(CriteriaType{Parameters{R"({"force_absolute_tolerance": -1.0})"}},
                                      "The force_absolute_tolerance must not be negative, got -1")
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_GlobalErrorIsResidualNormDividedByReferenceForceNorm,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // || r || / (CSP || f_ext^inact || + || f_int^act ||)
    KRATOS_EXPECT_DOUBLE_EQ(CriteriaType::CalculateGlobalForceError(3.0, 0.5, 100.0, 10.0), 0.05);
    KRATOS_EXPECT_DOUBLE_EQ(CriteriaType::CalculateGlobalForceError(3.0, 0.0, 100.0, 10.0), 0.3);
    KRATOS_EXPECT_DOUBLE_EQ(CriteriaType::CalculateGlobalForceError(3.0, 1.0, 100.0, 0.0), 0.03);
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_WithoutReferenceForcesAnyResidualIsAnInfiniteError,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_EQ(CriteriaType::CalculateGlobalForceError(1.0e-6, 1.0, 0.0, 0.0),
                     std::numeric_limits<double>::infinity());
    KRATOS_EXPECT_DOUBLE_EQ(CriteriaType::CalculateGlobalForceError(0.0, 1.0, 0.0, 0.0), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_OnlyFreeDisplacementDofsAreTakenIntoAccount, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    auto  dof_set      = CreateDofSet(r_model_part);
    auto  criteria     = CriteriaType{Parameters{R"({"echo_level": 0})"}};

    // The fixed displacement and the water pressure entries are large, but they must be ignored
    const auto inactive_external_forces = UblasUtilities::CreateVector({600.0, 800.0, 1.0e6, 1.0e6});
    const auto active_internal_forces   = UblasUtilities::CreateVector({0.0, 0.0, 1.0e6, 1.0e6});

    // Act & Assert: without soil stress points, the CSP is one, so the error is 5 / 1000
    criteria.SetReferenceForces(inactive_external_forces, active_internal_forces);
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part, dof_set, UblasUtilities::CreateVector({3.0, 4.0, 1.0e3, 1.0e3})))

    // The error is 50 / 1000
    criteria.SetReferenceForces(inactive_external_forces, active_internal_forces);
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part, dof_set, UblasUtilities::CreateVector({30.0, 40.0, 0.0, 0.0})))

    // The active internal forces are part of the reference as well: the error is 50 / (1000 + 8000)
    criteria.SetReferenceForces(inactive_external_forces, UblasUtilities::CreateVector({4800.0, 6400.0, 0.0, 0.0}));
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part, dof_set, UblasUtilities::CreateVector({30.0, 40.0, 0.0, 0.0})))
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_CurrentStiffnessParameterScalesInactiveExternalForces,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    auto  dof_set      = CreateDofSet(r_model_part);
    auto  p_soil_point = MakeSoilPoint(1.0, 4.0);
    AddTriangularElement(r_model_part, 1.0, {p_soil_point});
    auto       criteria                 = CriteriaType{Parameters{R"({"echo_level": 0})"}};
    const auto inactive_external_forces = UblasUtilities::CreateVector({600.0, 800.0, 0.0, 0.0});
    const auto active_internal_forces   = UblasUtilities::CreateVector({0.0, 0.0, 0.0, 0.0});
    const auto residual                 = UblasUtilities::CreateVector({3.0, 4.0, 0.0, 0.0});

    // Act & Assert: the CSP is 0.25, so the error is 5 / 250
    criteria.SetReferenceForces(inactive_external_forces, active_internal_forces);
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part, dof_set, residual))

    // Once the stress point behaves elastically, the error is 5 / 1000
    p_soil_point->SetStrainEnergyIncrements(4.0, 4.0);
    criteria.SetReferenceForces(inactive_external_forces, active_internal_forces);
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part, dof_set, residual))
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_ResidualBelowAbsoluteToleranceIsConverged, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto       model        = Model{};
    auto&      r_model_part = model.CreateModelPart("Main");
    auto       dof_set      = CreateDofSet(r_model_part);
    auto       criteria     = CriteriaType{Parameters{R"({"echo_level": 0, "force_absolute_tolerance": 1.0e-6})"}};
    const auto no_forces    = Vector{ZeroVector(4)};

    // Act & Assert
    criteria.SetReferenceForces(no_forces, no_forces);
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part, dof_set, UblasUtilities::CreateVector({1.0e-7, 0.0, 0.0, 0.0})))

    criteria.SetReferenceForces(no_forces, no_forces);
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part, dof_set, UblasUtilities::CreateVector({1.0e-5, 0.0, 0.0, 0.0})))
}

KRATOS_TEST_CASE_IN_SUITE(GeoGlobalForceErrorCriteria_ReferenceForcesMustBeSetBeforeEachEvaluation, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto       model        = Model{};
    auto&      r_model_part = model.CreateModelPart("Main");
    auto       dof_set      = CreateDofSet(r_model_part);
    auto       criteria     = CriteriaType{Parameters{R"({"echo_level": 0})"}};
    const auto forces       = UblasUtilities::CreateVector({600.0, 800.0, 0.0, 0.0});
    const auto residual     = Vector{ZeroVector(4)};

    // Act & Assert
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(PostCriteria(criteria, r_model_part, dof_set, residual),
                                      "The reference forces of the GeoGlobalForceErrorCriteria have not been set")

    criteria.SetReferenceForces(forces, forces);
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part, dof_set, residual))
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(PostCriteria(criteria, r_model_part, dof_set, residual),
                                      "The reference forces of the GeoGlobalForceErrorCriteria have not been set")
}

KRATOS_TEST_CASE_IN_SUITE(CurrentStiffnessParameterCalculator_IsRatioOfTotalAndElasticStrainEnergyIncrements,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(1.0, 4.0), 0.25);

    // The CSP is limited to the range [0, 1]
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(5.0, 4.0), 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(-1.0, 4.0), 0.0);

    // Without any elastic strain energy increment, the behaviour is regarded as elastic
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(0.0, 0.0), 1.0);
}

KRATOS_TEST_CASE_IN_SUITE(CurrentStiffnessParameterCalculator_WeighsSoilStressPointsByTheirVolume, KratosGeoMechanicsFastSuite)
{
    // Arrange: an element of area 0.5 with one stress point and an element of area 2 with two stress
    // points, i.e. all stress points represent a volume of 0.5
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    AddTriangularElement(r_model_part, 1.0, {MakeSoilPoint(1.0, 1.0)});
    AddTriangularElement(r_model_part, 2.0, {MakeSoilPoint(0.0, 1.0), MakeSoilPoint(0.0, 1.0)});

    // Act & Assert: (0.5 * 1) / (0.5 * 1 + 1 * 1 + 1 * 1)
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(r_model_part), 0.2);
}

KRATOS_TEST_CASE_IN_SUITE(CurrentStiffnessParameterCalculator_IgnoresInterfacePointsAndInactiveElements,
                          KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    AddTriangularElement(r_model_part, 1.0, {MakeSoilPoint(1.0, 2.0)});
    AddTriangularElement(r_model_part, 1.0,
                         {std::make_shared<StubLocalErrorLaw>(Geo::StressPointType::Interface, 0.0, 1.0)});
    auto& r_inactive_element = AddTriangularElement(r_model_part, 1.0, {MakeSoilPoint(0.0, 1.0)});
    r_inactive_element.Set(ACTIVE, false);

    // Act & Assert
    KRATOS_EXPECT_DOUBLE_EQ(CurrentStiffnessParameterCalculator::Calculate(r_model_part), 0.5);
}

} // namespace Kratos::Testing

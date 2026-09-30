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
#include "custom_strategies/convergence_criteria/geo_local_error_criteria.h"
#include "custom_strategies/convergence_criteria/local_error_evaluator.h"
#include "custom_utilities/ublas_utilities.h"
#include "geometries/point_2d.h"
#include "includes/element.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

#include <sstream>

using namespace Kratos;

namespace
{

using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;
using LocalSpaceType  = UblasSpace<double, Matrix, Vector>;
using CriteriaType    = GeoLocalErrorCriteria<SparseSpaceType, LocalSpaceType>;

class StubLocalErrorLaw : public ConstitutiveLaw, public LocalErrorDataProvider
{
public:
    StubLocalErrorLaw(std::optional<Geo::StressPointType> Type, Geo::LocalErrorData Data)
        : mType(Type), mData(std::move(Data))
    {
    }

    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override
    {
        return std::make_shared<StubLocalErrorLaw>(mType, mData);
    }

    [[nodiscard]] std::optional<Geo::StressPointType> GetStressPointType() const override
    {
        return mType;
    }

    [[nodiscard]] Geo::LocalErrorData CalculateLocalErrorData(Parameters&) override { return mData; }

    void SetElasticPredictorDeviation(const Vector& rDeviation)
    {
        mData.ElasticPredictorDeviation = rDeviation;
    }

private:
    std::optional<Geo::StressPointType> mType;
    Geo::LocalErrorData                 mData;
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

Geo::LocalErrorData MakeLocalErrorData(bool          IsPlastic,
                                       const Vector& rDeviation,
                                       bool          HasStressDependentStiffness = false,
                                       double        MaximumShearStress          = 0.0,
                                       double        Cohesion                    = 0.0)
{
    auto result                        = Geo::LocalErrorData{};
    result.IsPlastic                   = IsPlastic;
    result.HasStressDependentStiffness = HasStressDependentStiffness;
    result.ElasticPredictorDeviation   = rDeviation;
    result.MaximumShearStress          = MaximumShearStress;
    result.Cohesion                    = Cohesion;
    return result;
}

std::shared_ptr<StubLocalErrorLaw> MakeStubLaw(std::optional<Geo::StressPointType> Type,
                                               const Geo::LocalErrorData&          rData)
{
    return std::make_shared<StubLocalErrorLaw>(Type, rData);
}

Element& AddElement(ModelPart& rModelPart, std::vector<ConstitutiveLaw::Pointer> ConstitutiveLaws)
{
    const auto id           = rModelPart.NumberOfElements() + 1;
    auto       p_properties = rModelPart.HasProperties(0) ? rModelPart.pGetProperties(0)
                                                          : rModelPart.CreateNewProperties(0);
    auto       p_geometry = Kratos::make_shared<Point2D<Node>>(rModelPart.CreateNewNode(id, 0.0, 0.0, 0.0));
    auto       p_element  = Kratos::make_intrusive<ElementWithConstitutiveLaws>(
        id, p_geometry, p_properties, std::move(ConstitutiveLaws));
    rModelPart.AddElement(p_element);
    return *p_element;
}

// Adds a single element with the given number of plastic soil points, which all have the given deviation
std::vector<std::shared_ptr<StubLocalErrorLaw>> AddPlasticSoilPoints(ModelPart&    rModelPart,
                                                                     std::size_t   NumberOfPoints,
                                                                     const Vector& rDeviation)
{
    auto stub_laws         = std::vector<std::shared_ptr<StubLocalErrorLaw>>{};
    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    for (std::size_t i = 0; i < NumberOfPoints; ++i) {
        stub_laws.push_back(MakeStubLaw(Geo::StressPointType::Soil, MakeLocalErrorData(true, rDeviation)));
        constitutive_laws.push_back(stub_laws.back());
    }
    AddElement(rModelPart, constitutive_laws);
    return stub_laws;
}

void InitializeSolutionStep(CriteriaType& rCriteria, ModelPart& rModelPart)
{
    auto dof_set = ModelPart::DofsArrayType{};
    auto A       = CriteriaType::TSystemMatrixType{};
    auto Dx      = CriteriaType::TSystemVectorType{};
    auto b       = CriteriaType::TSystemVectorType{};
    rCriteria.InitializeSolutionStep(rModelPart, dof_set, A, Dx, b);
}

bool PostCriteria(CriteriaType& rCriteria, ModelPart& rModelPart)
{
    auto dof_set = ModelPart::DofsArrayType{};
    auto A       = CriteriaType::TSystemMatrixType{};
    auto Dx      = CriteriaType::TSystemVectorType{};
    auto b       = CriteriaType::TSystemVectorType{};
    return rCriteria.PostCriteria(rModelPart, dof_set, A, Dx, b);
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoLocalErrorCriteria_DefaultSettingsFollowPlaxisDefaults, KratosGeoMechanicsFastSuite)
{
    // Act
    const auto criteria = CriteriaType{};

    // Assert
    const auto& r_settings = criteria.GetEvaluator().GetSettings();
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedSoilPlasticLocalError, 0.01);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedSoilElasticLocalError, 0.01);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedInterfacePlasticLocalError, 0.01);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedInaccurateSoilPlasticPointFraction, 0.1);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedInaccurateSoilElasticPointFraction, 0.1);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedInaccurateInterfacePlasticPointFraction, 0.1);
    KRATOS_EXPECT_EQ(r_settings.NumberOfAlwaysToleratedInaccuratePoints, std::size_t{3});
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.MinimumReferenceStress, 1.0e3);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ReferencePressure, 1.0e5);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLocalErrorCriteria_RequestsRebuildOfRightHandSide, KratosGeoMechanicsFastSuite)
{
    auto criteria = CriteriaType{};
    KRATOS_EXPECT_TRUE(criteria.GetActualizeRHSflag())
}

KRATOS_TEST_CASE_IN_SUITE(GeoLocalErrorCriteria_ToleratedErrorIsDefaultOfAllToleratedLocalErrors, KratosGeoMechanicsFastSuite)
{
    // Act
    const auto criteria = CriteriaType{Parameters{R"(
    {
        "tolerated_error"                        : 0.05,
        "tolerated_soil_elastic_local_error"     : 0.2,
        "number_of_always_tolerated_inaccurate_points": 0,
        "minimum_reference_stress"               : 1.0,
        "reference_pressure"                     : 100.0
    })"}};

    // Assert
    const auto& r_settings = criteria.GetEvaluator().GetSettings();
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedSoilPlasticLocalError, 0.05);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedSoilElasticLocalError, 0.2);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ToleratedInterfacePlasticLocalError, 0.05);
    KRATOS_EXPECT_EQ(r_settings.NumberOfAlwaysToleratedInaccuratePoints, std::size_t{0});
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.MinimumReferenceStress, 1.0);
    KRATOS_EXPECT_DOUBLE_EQ(r_settings.ReferencePressure, 100.0);
}

KRATOS_TEST_CASE_IN_SUITE(GeoLocalErrorCriteria_ThrowsForNegativeNumberOfAlwaysToleratedInaccuratePoints,
                          KratosGeoMechanicsFastSuite)
{
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        CriteriaType{Parameters{R"({"number_of_always_tolerated_inaccurate_points": -1})"}},
        "The number_of_always_tolerated_inaccurate_points must not be negative")
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_ToleratesAFractionOfInaccuratePointsPlusThree,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};

    counts.SoilPlastic           = 20;
    counts.InaccurateSoilPlastic = 4; // less than 0.1 * 20 + 3
    KRATOS_EXPECT_TRUE(evaluator.IsConverged(counts))
    counts.InaccurateSoilPlastic = 5;
    KRATOS_EXPECT_FALSE(evaluator.IsConverged(counts))

    counts                            = LocalErrorEvaluator::PointCounts{};
    counts.InaccurateInterfacePlastic = 3; // less than 0.1 * 3 + 3
    counts.InterfacePlastic           = 3;
    KRATOS_EXPECT_TRUE(evaluator.IsConverged(counts))
    counts.InaccurateInterfacePlastic = 4; // not less than 0.1 * 4 + 3
    counts.InterfacePlastic           = 4;
    KRATOS_EXPECT_FALSE(evaluator.IsConverged(counts))

    counts                       = LocalErrorEvaluator::PointCounts{};
    counts.SoilElastic           = 1000;
    counts.InaccurateSoilElastic = 102; // less than 0.1 * 1000 + 3
    KRATOS_EXPECT_TRUE(evaluator.IsConverged(counts))
    counts.InaccurateSoilElastic = 103;
    KRATOS_EXPECT_FALSE(evaluator.IsConverged(counts))
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_PlasticSoilPointsAreInaccurateWhenErrorExceedsTolerance,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};

    // The reference stress is the minimum reference stress (1000), so the local errors are 0.009 and 0.011
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({9.0, 0.0})), Vector{}, counts);
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({11.0, 0.0})), Vector{}, counts);
    // A larger maximum shear stress increases the reference stress (local error = 0.0055)
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({11.0, 0.0}), false, 2000.0),
                             Vector{}, counts);
    // Only the change with respect to the previous iteration counts (local error = 0.001)
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({11.0, 0.0})),
                             UblasUtilities::CreateVector({10.0, 0.0}), counts);

    auto expected_counts                  = LocalErrorEvaluator::PointCounts{};
    expected_counts.SoilPlastic           = 4;
    expected_counts.InaccurateSoilPlastic = 1;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_OnlyStressDependentElasticSoilPointsCanBeInaccurate,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};

    // The reference stress of elastic points is 1/200 of the reference pressure (500), so the local
    // errors are 0.008 and 0.012
    constexpr auto has_stress_dependent_stiffness = true;
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(false, UblasUtilities::CreateVector({4.0, 0.0}), has_stress_dependent_stiffness),
                             Vector{}, counts);
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(false, UblasUtilities::CreateVector({6.0, 0.0}), has_stress_dependent_stiffness),
                             Vector{}, counts);
    // An elastic point with a linear elastic stiffness is never counted as inaccurate
    evaluator.AddStressPoint(Geo::StressPointType::Soil,
                             MakeLocalErrorData(false, UblasUtilities::CreateVector({6.0, 0.0})), Vector{}, counts);

    auto expected_counts                  = LocalErrorEvaluator::PointCounts{};
    expected_counts.SoilElastic           = 3;
    expected_counts.InaccurateSoilElastic = 1;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_OnlyPlasticInterfacePointsAreCounted, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};

    evaluator.AddStressPoint(Geo::StressPointType::Interface,
                             MakeLocalErrorData(false, UblasUtilities::CreateVector({100.0, 0.0})), Vector{}, counts);
    // Local errors are 0.009 and 0.011 (with respect to the cohesion)
    evaluator.AddStressPoint(Geo::StressPointType::Interface,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({0.0, 18.0}), false, 10.0, 2000.0),
                             Vector{}, counts);
    evaluator.AddStressPoint(Geo::StressPointType::Interface,
                             MakeLocalErrorData(true, UblasUtilities::CreateVector({0.0, 22.0}), false, 10.0, 2000.0),
                             Vector{}, counts);

    auto expected_counts                       = LocalErrorEvaluator::PointCounts{};
    expected_counts.InterfacePlastic           = 2;
    expected_counts.InaccurateInterfacePlastic = 1;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_OnlyShearTractionsDetermineAccuracyOfInterfacePoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};

    // The change of the normal traction (e.g. of an opening interface) is not taken into account.
    // The local errors are 0.009 and 0.011 (with respect to the minimum reference stress).
    auto accurate_data                       = MakeLocalErrorData(true, UblasUtilities::CreateVector({500.0, 9.0}));
    accurate_data.IndexOfFirstShearComponent = 1;
    evaluator.AddStressPoint(Geo::StressPointType::Interface, accurate_data, Vector{}, counts);
    auto inaccurate_data = MakeLocalErrorData(true, UblasUtilities::CreateVector({500.0, 11.0}));
    inaccurate_data.IndexOfFirstShearComponent = 1;
    evaluator.AddStressPoint(Geo::StressPointType::Interface, inaccurate_data, Vector{}, counts);

    auto expected_counts                       = LocalErrorEvaluator::PointCounts{};
    expected_counts.InterfacePlastic           = 2;
    expected_counts.InaccurateInterfacePlastic = 1;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_EvaluatesAllStressPointsOfTheModelPart, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    AddPlasticSoilPoints(r_model_part, 10, UblasUtilities::CreateVector({50.0, 0.0}));
    AddElement(r_model_part,
               {MakeStubLaw(Geo::StressPointType::Soil, MakeLocalErrorData(false, ZeroVector(2))),
                MakeStubLaw(Geo::StressPointType::Interface,
                            MakeLocalErrorData(true, UblasUtilities::CreateVector({50.0, 0.0}))),
                MakeStubLaw(std::nullopt, MakeLocalErrorData(true, UblasUtilities::CreateVector({50.0, 0.0})))});
    auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};

    // Act
    const auto counts = evaluator.EvaluateAndStoreIterationState(r_model_part);

    // Assert
    auto expected_counts                       = LocalErrorEvaluator::PointCounts{};
    expected_counts.SoilPlastic                = 10;
    expected_counts.InaccurateSoilPlastic      = 10;
    expected_counts.SoilElastic                = 1;
    expected_counts.InterfacePlastic           = 1;
    expected_counts.InaccurateInterfacePlastic = 1;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
    KRATOS_EXPECT_FALSE(evaluator.IsConverged(counts))
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_CountsStressPointsOfInactiveSoilElementsAsElastic, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    auto& r_element    = AddElement(
        r_model_part, {MakeStubLaw(Geo::StressPointType::Soil,
                                      MakeLocalErrorData(true, UblasUtilities::CreateVector({50.0, 0.0}))),
                          MakeStubLaw(Geo::StressPointType::Soil,
                                      MakeLocalErrorData(true, UblasUtilities::CreateVector({50.0, 0.0}))),
                          MakeStubLaw(Geo::StressPointType::Interface,
                                      MakeLocalErrorData(true, UblasUtilities::CreateVector({50.0, 0.0})))});
    r_element.Set(ACTIVE, false);
    auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};

    // Act
    const auto counts = evaluator.EvaluateAndStoreIterationState(r_model_part);

    // Assert
    auto expected_counts        = LocalErrorEvaluator::PointCounts{};
    expected_counts.SoilElastic = 2;
    KRATOS_EXPECT_TRUE(counts == expected_counts)
}

KRATOS_TEST_CASE_IN_SUITE(GeoLocalErrorCriteria_ComparesWithPreviousIterationWithinAStep, KratosGeoMechanicsFastSuite)
{
    // Arrange
    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main");
    const auto stub_laws = AddPlasticSoilPoints(r_model_part, 10, UblasUtilities::CreateVector({50.0, 0.0}));
    auto       criteria  = CriteriaType{Parameters{R"({"echo_level": 0})"}};

    // Act & Assert: in the first iteration, the deviations are compared with the start of the step
    InitializeSolutionStep(criteria, r_model_part);
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part))

    // The deviations did not change, so all plastic points are accurate
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part))

    // Four of the ten points change significantly, which is not tolerated (4 >= 0.1 * 10 + 3)
    for (std::size_t i = 0; i < 4; ++i) {
        stub_laws[i]->SetElasticPredictorDeviation(UblasUtilities::CreateVector({80.0, 0.0}));
    }
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part))

    // Three of the ten points change significantly, which is tolerated (3 < 0.1 * 10 + 3)
    for (std::size_t i = 0; i < 3; ++i) {
        stub_laws[i]->SetElasticPredictorDeviation(UblasUtilities::CreateVector({50.0, 0.0}));
    }
    KRATOS_EXPECT_TRUE(PostCriteria(criteria, r_model_part))

    // A new step starts from the state at the start of the step again
    InitializeSolutionStep(criteria, r_model_part);
    KRATOS_EXPECT_FALSE(PostCriteria(criteria, r_model_part))
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorEvaluator_PrintsCountsAndToleratedNumbers, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto evaluator = LocalErrorEvaluator{LocalErrorEvaluator::Settings{}};
    auto       counts    = LocalErrorEvaluator::PointCounts{};
    counts.SoilPlastic           = 10;
    counts.InaccurateSoilPlastic = 5;

    auto output = std::ostringstream{};
    evaluator.PrintCounts(output, counts);

    KRATOS_EXPECT_TRUE(output.str().find("Soil plastic points:      10, inaccurate: 5 (tolerated: < 4)") !=
                       std::string::npos)
}

} // namespace Kratos::Testing

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

#include "custom_constitutive/incremental_linear_elastic_interface_law.h"
#include "custom_constitutive/incremental_linear_elastic_law.h"
#include "custom_constitutive/interface_coulomb_with_tension_cut_off.h"
#include "custom_constitutive/interface_plane_strain.h"
#include "custom_constitutive/linear_elastic_2D_beam_law.h"
#include "custom_constitutive/linear_elastic_2D_interface_law.h"
#include "custom_constitutive/linear_elastic_3D_interface_law.h"
#include "custom_constitutive/local_error_data_provider.h"
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_constitutive/plane_strain.h"
#include "custom_constitutive/three_dimensional.h"
#include "custom_utilities/local_error_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace
{

Properties MakeMohrCoulombProperties()
{
    auto result = Properties{};
    result.SetValue(GEO_FRICTION_ANGLE, 30.0);
    result.SetValue(GEO_COHESION, 10.0);
    result.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    result.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    result.SetValue(YOUNG_MODULUS, 1.0e3);
    result.SetValue(POISSON_RATIO, 0.2);
    return result;
}

Properties MakeInterfaceCoulombProperties()
{
    auto result = Properties{};
    result.SetValue(GEO_FRICTION_ANGLE, 30.0);
    result.SetValue(GEO_COHESION, 10.0);
    result.SetValue(GEO_DILATANCY_ANGLE, 0.0);
    result.SetValue(GEO_TENSILE_STRENGTH, 10.0);
    result.SetValue(INTERFACE_NORMAL_STIFFNESS, 1.0e3);
    result.SetValue(INTERFACE_SHEAR_STIFFNESS, 1.0e3);
    return result;
}

// Initializes the law with a zero stress and zero strain at the start of the step
void InitializeLawWithZeroState(ConstitutiveLaw& rLaw, ConstitutiveLaw::Parameters& rParameters, std::size_t StrainSize)
{
    const auto dummy_element_geometry      = Geometry<Node>{};
    const auto dummy_shape_function_values = Vector{};
    rLaw.InitializeMaterial(rParameters.GetMaterialProperties(), dummy_element_geometry,
                            dummy_shape_function_values);

    auto zero_state = Vector{ZeroVector{StrainSize}};
    rParameters.SetStrainVector(zero_state);
    rParameters.SetStressVector(zero_state);
    rLaw.InitializeMaterialResponseCauchy(rParameters);
}

Vector CalculateStress(ConstitutiveLaw& rLaw, ConstitutiveLaw::Parameters& rParameters, Vector StrainVector)
{
    auto stress_vector = Vector{ZeroVector{StrainVector.size()}};
    rParameters.SetStrainVector(StrainVector);
    rParameters.SetStressVector(stress_vector);
    rLaw.CalculateMaterialResponseCauchy(rParameters);
    return rParameters.GetStressVector();
}

LocalErrorDataProvider& AsLocalErrorDataProvider(ConstitutiveLaw& rLaw)
{
    auto* p_provider = dynamic_cast<LocalErrorDataProvider*>(&rLaw);
    KRATOS_ERROR_IF_NOT(p_provider) << "The constitutive law does not provide local error data\n";
    return *p_provider;
}

} // namespace

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_LocalErrorDataOfElasticStressPoint,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto properties = MakeMohrCoulombProperties();
    auto       parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto law = MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()};
    InitializeLawWithZeroState(law, parameters, 4);

    // Act
    const auto stress_vector = CalculateStress(law, parameters, UblasUtilities::CreateVector({-1.0e-3, -2.0e-3, 0.0, 0.0}));
    const auto data = law.CalculateLocalErrorData(parameters);

    // Assert
    KRATOS_EXPECT_TRUE(law.GetStressPointType() == Geo::StressPointType::Soil)
    KRATOS_EXPECT_FALSE(data.IsPlastic)
    KRATOS_EXPECT_FALSE(data.HasStressDependentStiffness)
    KRATOS_EXPECT_VECTOR_NEAR(data.ElasticPredictorDeviation, ZeroVector(4), Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(data.MaximumShearStress, LocalErrorUtilities::CalculateMaximumShearStress(stress_vector),
                       Defaults::absolute_tolerance);
    KRATOS_EXPECT_DOUBLE_EQ(data.Cohesion, 10.0);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_LocalErrorDataOfPlasticStressPoint,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto properties = MakeMohrCoulombProperties();
    auto       parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto law = MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()};
    InitializeLawWithZeroState(law, parameters, 4);

    // Act
    const auto strain_vector = UblasUtilities::CreateVector({0.0, -0.01, 0.0, 0.06});
    const auto stress_vector = CalculateStress(law, parameters, strain_vector);
    const auto data          = law.CalculateLocalErrorData(parameters);

    // Assert: the deviation is the difference between the elastic predictor and the constitutive stress
    const auto   elastic_matrix    = PlaneStrain{}.CalculateElasticMatrix(properties);
    const Vector elastic_predictor = prod(elastic_matrix, strain_vector);
    KRATOS_EXPECT_TRUE(data.IsPlastic)
    KRATOS_EXPECT_VECTOR_NEAR(data.ElasticPredictorDeviation, Vector{elastic_predictor - stress_vector},
                              Defaults::absolute_tolerance);
    KRATOS_EXPECT_GT(norm_2(data.ElasticPredictorDeviation), 1.0);
    KRATOS_EXPECT_NEAR(data.MaximumShearStress, LocalErrorUtilities::CalculateMaximumShearStress(stress_vector),
                       Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_ChangeOfDeviationEqualsDifferenceBetweenEquilibriumAndConstitutiveStress,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto properties = MakeMohrCoulombProperties();
    auto       parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto law = MohrCoulombWithTensionCutOff{std::make_unique<ThreeDimensional>()};
    InitializeLawWithZeroState(law, parameters, 6);

    // Two subsequent iterations of the same step (no finalization in between)
    const auto strain_of_previous_iteration = UblasUtilities::CreateVector({0.0, -0.005, 0.0, 0.03, 0.0, 0.0});
    const auto strain_of_current_iteration = UblasUtilities::CreateVector({0.002, -0.01, 0.0, 0.06, 0.01, 0.0});

    // Act
    const auto constitutive_stress_of_previous_iteration = CalculateStress(law, parameters, strain_of_previous_iteration);
    const auto deviation_of_previous_iteration = law.CalculateLocalErrorData(parameters).ElasticPredictorDeviation;
    const auto constitutive_stress_of_current_iteration = CalculateStress(law, parameters, strain_of_current_iteration);
    const auto deviation_of_current_iteration = law.CalculateLocalErrorData(parameters).ElasticPredictorDeviation;

    // Assert: sigma_eq,j - sigma_c,j = P_j - P_j-1, where sigma_eq,j = sigma_c,j-1 + D^e delta_eps_j
    const auto   elastic_matrix = ThreeDimensional{}.CalculateElasticMatrix(properties);
    const Vector equilibrium_stress =
        constitutive_stress_of_previous_iteration +
        prod(elastic_matrix, Vector{strain_of_current_iteration - strain_of_previous_iteration});
    KRATOS_EXPECT_VECTOR_NEAR(Vector{equilibrium_stress - constitutive_stress_of_current_iteration},
                              Vector{deviation_of_current_iteration - deviation_of_previous_iteration}, 1.0e-9);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_ImposedStressIsRegardedAsElastic, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto properties = MakeMohrCoulombProperties();
    auto       parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto law = MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()};
    InitializeLawWithZeroState(law, parameters, 4);
    [[maybe_unused]] const auto stress_vector =
        CalculateStress(law, parameters, UblasUtilities::CreateVector({0.0, -0.01, 0.0, 0.06}));

    // Act
    law.SetValue(CAUCHY_STRESS_VECTOR, UblasUtilities::CreateVector({-10.0, -20.0, -10.0, 1.0}), ProcessInfo{});
    const auto data = law.CalculateLocalErrorData(parameters);

    // Assert
    KRATOS_EXPECT_FALSE(data.IsPlastic)
    KRATOS_EXPECT_VECTOR_NEAR(data.ElasticPredictorDeviation, ZeroVector(4), Defaults::absolute_tolerance);
}

KRATOS_TEST_CASE_IN_SUITE(MohrCoulombWithTensionCutOff_ClonedLawProvidesLocalErrorData, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto law      = MohrCoulombWithTensionCutOff{std::make_unique<PlaneStrain>()};
    const auto p_cloned = law.Clone();
    KRATOS_EXPECT_NE(dynamic_cast<const LocalErrorDataProvider*>(p_cloned.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(InterfaceCoulombWithTensionCutOff_LocalErrorDataOfElasticAndPlasticStressPoints,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Arrange
    const auto properties = MakeInterfaceCoulombProperties();
    auto       parameters = ConstitutiveLaw::Parameters{};
    parameters.SetMaterialProperties(properties);
    parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
    auto law = InterfaceCoulombWithTensionCutOff{std::make_unique<InterfacePlaneStrain>()};
    InitializeLawWithZeroState(law, parameters, 2);

    // Act: elastic state
    const auto elastic_traction = CalculateStress(law, parameters, UblasUtilities::CreateVector({-0.01, 0.002}));
    const auto elastic_data     = law.CalculateLocalErrorData(parameters);

    // Assert
    KRATOS_EXPECT_TRUE(law.GetStressPointType() == Geo::StressPointType::Interface)
    KRATOS_EXPECT_FALSE(elastic_data.IsPlastic)
    KRATOS_EXPECT_VECTOR_NEAR(elastic_data.ElasticPredictorDeviation, ZeroVector(2), Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(elastic_data.MaximumShearStress, std::abs(elastic_traction[1]), Defaults::absolute_tolerance);
    KRATOS_EXPECT_DOUBLE_EQ(elastic_data.Cohesion, 10.0);

    // Act: plastic state (the trial traction is (-10, 50))
    const auto plastic_traction = CalculateStress(law, parameters, UblasUtilities::CreateVector({-0.01, 0.05}));
    const auto plastic_data     = law.CalculateLocalErrorData(parameters);

    // Assert
    KRATOS_EXPECT_TRUE(plastic_data.IsPlastic)
    KRATOS_EXPECT_VECTOR_NEAR(plastic_data.ElasticPredictorDeviation,
                              Vector{UblasUtilities::CreateVector({-10.0, 50.0}) - plastic_traction}, 1.0e-9);
    KRATOS_EXPECT_NEAR(plastic_data.MaximumShearStress, std::abs(plastic_traction[1]), Defaults::absolute_tolerance);
    // Only the shear traction is used by the local error of an interface point
    KRATOS_EXPECT_EQ(plastic_data.IndexOfFirstShearComponent, std::size_t{1});
}

KRATOS_TEST_CASE_IN_SUITE(LinearElasticLaws_LocalErrorDataIsAlwaysElasticWithoutDeviation,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto parameters = ConstitutiveLaw::Parameters{};

    auto soil_law = GeoIncrementalLinearElasticLaw{std::make_unique<PlaneStrain>()};
    KRATOS_EXPECT_TRUE(soil_law.GetStressPointType() == Geo::StressPointType::Soil)
    const auto soil_data = soil_law.CalculateLocalErrorData(parameters);
    KRATOS_EXPECT_FALSE(soil_data.IsPlastic)
    KRATOS_EXPECT_FALSE(soil_data.HasStressDependentStiffness)
    KRATOS_EXPECT_VECTOR_NEAR(soil_data.ElasticPredictorDeviation, ZeroVector(4), Defaults::absolute_tolerance);

    auto interface_law = GeoIncrementalLinearElasticInterfaceLaw{std::make_unique<InterfacePlaneStrain>()};
    KRATOS_EXPECT_TRUE(interface_law.GetStressPointType() == Geo::StressPointType::Interface)
    const auto interface_data = interface_law.CalculateLocalErrorData(parameters);
    KRATOS_EXPECT_FALSE(interface_data.IsPlastic)
    KRATOS_EXPECT_VECTOR_NEAR(interface_data.ElasticPredictorDeviation, ZeroVector(2), Defaults::absolute_tolerance);

    KRATOS_EXPECT_TRUE(LinearElastic2DInterfaceLaw{}.GetStressPointType() == Geo::StressPointType::Interface)
    KRATOS_EXPECT_TRUE(LinearElastic3DInterfaceLaw{}.GetStressPointType() == Geo::StressPointType::Interface)
}

KRATOS_TEST_CASE_IN_SUITE(LinearElastic2DBeamLaw_DoesNotTakePartInLocalErrorCriteria, KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_FALSE(LinearElastic2DBeamLaw{}.GetStressPointType().has_value())
}

} // namespace Kratos::Testing

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

#include "custom_utilities/local_error_utilities.h"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"
#include "tests/cpp_tests/test_utilities.h"

using namespace Kratos;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_MaximumShearStressIsRadiusOfLargestMohrCircle,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // Plane strain stress vectors (xx, yy, zz, xy)
    KRATOS_EXPECT_NEAR(LocalErrorUtilities::CalculateMaximumShearStress(
                           UblasUtilities::CreateVector({-10.0, -30.0, -20.0, 0.0})),
                       10.0, Defaults::absolute_tolerance);
    KRATOS_EXPECT_NEAR(LocalErrorUtilities::CalculateMaximumShearStress(
                           UblasUtilities::CreateVector({0.0, 0.0, 0.0, 5.0})),
                       5.0, Defaults::absolute_tolerance);

    // Three-dimensional stress vector (xx, yy, zz, xy, yz, xz)
    KRATOS_EXPECT_NEAR(LocalErrorUtilities::CalculateMaximumShearStress(
                           UblasUtilities::CreateVector({1.0, 2.0, 3.0, 0.0, 0.0, 0.0})),
                       1.0, Defaults::absolute_tolerance);

    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateMaximumShearStress(Vector{}), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_ShearTractionMagnitudeOnlyAccountsForShearComponents,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateShearTractionMagnitude(
                                UblasUtilities::CreateVector({-10.0, -3.0}), 1),
                            3.0);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateShearTractionMagnitude(
                                UblasUtilities::CreateVector({-10.0, 3.0, 4.0}), 1),
                            5.0);
    // Interface components that are mapped onto three-dimensional stress components
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateShearTractionMagnitude(
                                UblasUtilities::CreateVector({0.0, 0.0, -10.0, 1.0, 2.0, 2.0}), 3),
                            3.0);
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_CohesionIsTakenFromGeoCohesionOrUMatParameters,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    auto properties = Properties{};
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::GetCohesionIfAvailable(properties), 0.0);

    properties.SetValue(UMAT_PARAMETERS, UblasUtilities::CreateVector({1.0e4, 7.0, 30.0}));
    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 2);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::GetCohesionIfAvailable(properties), 7.0);

    properties.SetValue(INDEX_OF_UMAT_C_PARAMETER, 4);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::GetCohesionIfAvailable(properties), 0.0);

    properties.SetValue(GEO_COHESION, 12.0);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::GetCohesionIfAvailable(properties), 12.0);
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_LocalErrorIsChangeOfDeviationOverReferenceStress,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto deviation = UblasUtilities::CreateVector({3.0, 4.0});

    // Without a previous deviation, the change equals the deviation itself (norm = 5)
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(deviation, Vector{}, 2.0, 1.0, 0.5), 2.5);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(deviation, Vector{}, 2.0, 10.0, 0.5), 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(deviation, Vector{}, 2.0, 1.0, 100.0), 0.05);

    // Only the change with respect to the previous deviation counts (norm = 4)
    KRATOS_EXPECT_DOUBLE_EQ(
        LocalErrorUtilities::CalculateLocalError(deviation, UblasUtilities::CreateVector({3.0, 0.0}), 2.0, 1.0, 0.5),
        2.0);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(deviation, deviation, 2.0, 1.0, 0.5), 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_LocalErrorOnlyAccountsForComponentsFromGivenIndex,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // e.g. an interface traction (normal, shear), of which only the shear component is considered
    const auto deviation = UblasUtilities::CreateVector({30.0, 4.0});

    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(deviation, Vector{}, 2.0, 1.0, 0.5, 1), 2.0);
    KRATOS_EXPECT_DOUBLE_EQ(LocalErrorUtilities::CalculateLocalError(
                                deviation, UblasUtilities::CreateVector({-30.0, 3.0}), 2.0, 1.0, 0.5, 1),
                            0.5);
}

KRATOS_TEST_CASE_IN_SUITE(LocalErrorUtilities_DeviationFromElasticPredictorIgnoresRoundOff,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    const auto stress_vector = UblasUtilities::CreateVector({-100.0, -50.0, -60.0, 10.0});

    KRATOS_EXPECT_FALSE(LocalErrorUtilities::DeviatesFromElasticPredictor(stress_vector, stress_vector))
    KRATOS_EXPECT_FALSE(LocalErrorUtilities::DeviatesFromElasticPredictor(
        stress_vector, Vector{stress_vector * (1.0 + 1.0e-13)}))
    KRATOS_EXPECT_TRUE(LocalErrorUtilities::DeviatesFromElasticPredictor(
        stress_vector, UblasUtilities::CreateVector({-90.0, -50.0, -60.0, 10.0})))
    KRATOS_EXPECT_FALSE(LocalErrorUtilities::DeviatesFromElasticPredictor(ZeroVector(4), ZeroVector(4)))
}

} // namespace Kratos::Testing

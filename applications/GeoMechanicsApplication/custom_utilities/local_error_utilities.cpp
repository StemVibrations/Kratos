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
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Kratos
{

double LocalErrorUtilities::CalculateMaximumShearStress(const Vector& rStressVector)
{
    if (rStressVector.empty()) return 0.0;

    Vector principal_stresses;
    Matrix eigen_vectors;
    StressStrainUtilities::CalculatePrincipalStresses(rStressVector, principal_stresses, eigen_vectors);

    // The principal stresses are sorted from the largest to the smallest one
    return 0.5 * (principal_stresses[0] - principal_stresses[principal_stresses.size() - 1]);
}

double LocalErrorUtilities::CalculateShearTractionMagnitude(const Vector& rTractionVector, std::size_t IndexOfFirstShearComponent)
{
    auto sum_of_squares = 0.0;
    for (auto i = IndexOfFirstShearComponent; i < rTractionVector.size(); ++i) {
        sum_of_squares += rTractionVector[i] * rTractionVector[i];
    }
    return std::sqrt(sum_of_squares);
}

double LocalErrorUtilities::GetCohesionIfAvailable(const Properties& rProperties)
{
    if (rProperties.Has(GEO_COHESION)) return rProperties[GEO_COHESION];

    if (rProperties.Has(UMAT_PARAMETERS) && rProperties.Has(INDEX_OF_UMAT_C_PARAMETER)) {
        const auto& r_umat_parameters = rProperties[UMAT_PARAMETERS];
        const auto  index             = rProperties[INDEX_OF_UMAT_C_PARAMETER]; // 1-based index
        if (index >= 1 && static_cast<std::size_t>(index) <= r_umat_parameters.size()) {
            return r_umat_parameters[index - 1];
        }
    }

    return 0.0;
}

double LocalErrorUtilities::CalculateLocalError(const Vector& rElasticPredictorDeviation,
                                                const Vector& rPreviousElasticPredictorDeviation,
                                                double        MaximumShearStress,
                                                double        Cohesion,
                                                double        MinimumReferenceStress,
                                                std::size_t   IndexOfFirstComponent)
{
    const auto has_previous_deviation =
        rPreviousElasticPredictorDeviation.size() == rElasticPredictorDeviation.size();
    auto sum_of_squares = 0.0;
    for (auto i = IndexOfFirstComponent; i < rElasticPredictorDeviation.size(); ++i) {
        const auto change = has_previous_deviation
                                ? rElasticPredictorDeviation[i] - rPreviousElasticPredictorDeviation[i]
                                : rElasticPredictorDeviation[i];
        sum_of_squares += change * change;
    }
    const auto numerator   = std::sqrt(sum_of_squares);
    const auto denominator = std::max({MaximumShearStress, Cohesion, MinimumReferenceStress});

    if (denominator <= 0.0) return numerator > 0.0 ? std::numeric_limits<double>::infinity() : 0.0;

    return numerator / denominator;
}

bool LocalErrorUtilities::DeviatesFromElasticPredictor(const Vector& rElasticPredictorStress,
                                                       const Vector& rConstitutiveStress)
{
    constexpr auto relative_tolerance = 1.0e-8;
    const auto     stress_scale = std::max(norm_2(rElasticPredictorStress), norm_2(rConstitutiveStress));
    return norm_2(rElasticPredictorStress - rConstitutiveStress) > relative_tolerance * stress_scale;
}

} // namespace Kratos

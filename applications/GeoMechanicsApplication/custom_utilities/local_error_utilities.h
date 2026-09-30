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

#pragma once

#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <cstddef>

namespace Kratos
{

class Properties;

class KRATOS_API(GEO_MECHANICS_APPLICATION) LocalErrorUtilities
{
public:
    // Returns the radius of the largest Mohr circle, i.e. (sigma_1 - sigma_3) / 2
    [[nodiscard]] static double CalculateMaximumShearStress(const Vector& rStressVector);

    // Returns the magnitude of the shear components of a traction vector, which are the components
    // starting at the given index
    [[nodiscard]] static double CalculateShearTractionMagnitude(const Vector& rTractionVector,
                                                                std::size_t   IndexOfFirstShearComponent);

    // Returns GEO_COHESION or the cohesion that is stored in the UMAT_PARAMETERS (at
    // INDEX_OF_UMAT_C_PARAMETER). When neither is available, zero is returned.
    [[nodiscard]] static double GetCohesionIfAvailable(const Properties& rProperties);

    // Returns || P_j - P_j-1 || / max(tau_max, c, MinimumReferenceStress), where P is the deviation
    // of the constitutive stress from the elastic predictor. An empty previous deviation (or one of
    // a different size) is treated as a zero vector. Only the components starting at the given
    // index are taken into account (e.g. the shear components of an interface traction).
    [[nodiscard]] static double CalculateLocalError(const Vector& rElasticPredictorDeviation,
                                                    const Vector& rPreviousElasticPredictorDeviation,
                                                    double        MaximumShearStress,
                                                    double        Cohesion,
                                                    double        MinimumReferenceStress,
                                                    std::size_t   IndexOfFirstComponent = 0);

    // Returns whether a constitutive stress differs (beyond round-off) from its elastic predictor
    [[nodiscard]] static bool DeviatesFromElasticPredictor(const Vector& rElasticPredictorStress,
                                                           const Vector& rConstitutiveStress);
};

} // namespace Kratos

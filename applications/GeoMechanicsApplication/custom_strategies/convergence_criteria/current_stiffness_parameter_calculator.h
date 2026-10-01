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

namespace Kratos
{

class ModelPart;

// Calculates the current stiffness parameter (CSP) of the global error criterion The CSP is the
// ratio of the total and the elastic strain energy
// increments of the soil stress points of the active elements:
//
//     CSP = sum_j (delta_eps_j . delta_sigma_j) dV_j / sum_j (delta_eps_j . D^e delta_eps_j) dV_j
//
// where the increments are taken with respect to the start of the step. The CSP equals one when all
// stress points behave elastically and it approaches zero at failure. Note that the ratio is sometimes
// defined the other way around, which would give a CSP that grows without bounds at failure.
// The volume dV_j of a stress point is taken as the domain size of its element divided by the
// number of stress points of that element.
class KRATOS_API(GEO_MECHANICS_APPLICATION) CurrentStiffnessParameterCalculator
{
public:
    // Evaluates the CSP of all soil stress points of the model part, based on the most recently
    // calculated material responses. Stress points of constitutive laws that don't implement
    // LocalErrorDataProvider are ignored.
    [[nodiscard]] static double Calculate(ModelPart& rModelPart);

    // Returns the CSP for the given (volume integrated) strain energy increments, limited to the
    // range [0, 1]. Without any elastic strain energy increment (e.g. when nothing deforms), the
    // behaviour is regarded as elastic.
    [[nodiscard]] static double Calculate(double TotalStrainEnergyIncrement, double ElasticStrainEnergyIncrement);
};

} // namespace Kratos

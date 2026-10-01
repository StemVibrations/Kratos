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

#include "custom_strategies/convergence_criteria/current_stiffness_parameter_calculator.h"
#include "custom_constitutive/local_error_data_provider.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

#include <algorithm>
#include <tuple>

namespace
{

using namespace Kratos;

// Returns the total and elastic strain energy increments of the soil stress points of an element,
// integrated over the element
std::tuple<double, double> CalculateStrainEnergyIncrements(Element& rElement, const ProcessInfo& rProcessInfo)
{
    if (!rElement.IsActive()) return {0.0, 0.0};

    auto constitutive_laws = std::vector<ConstitutiveLaw::Pointer>{};
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, rProcessInfo);
    if (constitutive_laws.empty()) return {0.0, 0.0};

    const auto volume_per_stress_point =
        rElement.GetGeometry().DomainSize() / static_cast<double>(constitutive_laws.size());

    auto total_strain_energy_increment   = 0.0;
    auto elastic_strain_energy_increment = 0.0;
    auto parameters = ConstitutiveLaw::Parameters{rElement.GetGeometry(), rElement.GetProperties(), rProcessInfo};
    for (const auto& rp_constitutive_law : constitutive_laws) {
        auto* p_provider = dynamic_cast<LocalErrorDataProvider*>(rp_constitutive_law.get());
        if (!p_provider || p_provider->GetStressPointType() != Geo::StressPointType::Soil) continue;

        const auto data = p_provider->CalculateLocalErrorData(parameters);
        total_strain_energy_increment += data.TotalStrainEnergyIncrement * volume_per_stress_point;
        elastic_strain_energy_increment += data.ElasticStrainEnergyIncrement * volume_per_stress_point;
    }

    return {total_strain_energy_increment, elastic_strain_energy_increment};
}

} // namespace

namespace Kratos
{

double CurrentStiffnessParameterCalculator::Calculate(ModelPart& rModelPart)
{
    const auto& r_process_info = rModelPart.GetProcessInfo();
    const auto [total_strain_energy_increment, elastic_strain_energy_increment] =
        block_for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(
            rModelPart.Elements(), [&r_process_info](Element& rElement) {
        return CalculateStrainEnergyIncrements(rElement, r_process_info);
    });

    return Calculate(total_strain_energy_increment, elastic_strain_energy_increment);
}

double CurrentStiffnessParameterCalculator::Calculate(double TotalStrainEnergyIncrement, double ElasticStrainEnergyIncrement)
{
    if (ElasticStrainEnergyIncrement <= 0.0) return 1.0;

    return std::clamp(TotalStrainEnergyIncrement / ElasticStrainEnergyIncrement, 0.0, 1.0);
}

} // namespace Kratos

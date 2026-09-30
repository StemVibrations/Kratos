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

#include "custom_strategies/convergence_criteria/local_error_evaluator.h"
#include "custom_utilities/local_error_utilities.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include <ostream>

namespace
{

using namespace Kratos;

class PointCountsReduction
{
public:
    using value_type  = LocalErrorEvaluator::PointCounts;
    using return_type = LocalErrorEvaluator::PointCounts;

    return_type mValue;

    [[nodiscard]] return_type GetValue() const { return mValue; }

    void LocalReduce(const value_type& rValue) { mValue += rValue; }

    void ThreadSafeReduce(const PointCountsReduction& rOther)
    {
        KRATOS_CRITICAL_SECTION
        mValue += rOther.mValue;
    }
};

std::vector<ConstitutiveLaw::Pointer> GetConstitutiveLaws(Element& rElement, const ProcessInfo& rProcessInfo)
{
    auto result = std::vector<ConstitutiveLaw::Pointer>{};
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, result, rProcessInfo);
    return result;
}

LocalErrorDataProvider* GetLocalErrorDataProvider(const ConstitutiveLaw::Pointer& rpConstitutiveLaw)
{
    return dynamic_cast<LocalErrorDataProvider*>(rpConstitutiveLaw.get());
}

} // namespace

namespace Kratos
{

LocalErrorEvaluator::PointCounts& LocalErrorEvaluator::PointCounts::operator+=(const PointCounts& rOther)
{
    SoilPlastic += rOther.SoilPlastic;
    InaccurateSoilPlastic += rOther.InaccurateSoilPlastic;
    SoilElastic += rOther.SoilElastic;
    InaccurateSoilElastic += rOther.InaccurateSoilElastic;
    InterfacePlastic += rOther.InterfacePlastic;
    InaccurateInterfacePlastic += rOther.InaccurateInterfacePlastic;
    return *this;
}

LocalErrorEvaluator::LocalErrorEvaluator(const Settings& rSettings) : mSettings(rSettings) {}

const LocalErrorEvaluator::Settings& LocalErrorEvaluator::GetSettings() const { return mSettings; }

void LocalErrorEvaluator::ResetPreviousIterationState() { mPreviousElasticPredictorDeviations.clear(); }

LocalErrorEvaluator::PointCounts LocalErrorEvaluator::EvaluateAndStoreIterationState(ModelPart& rModelPart)
{
    const auto number_of_elements = rModelPart.NumberOfElements();
    if (mPreviousElasticPredictorDeviations.size() != number_of_elements) {
        mPreviousElasticPredictorDeviations.assign(number_of_elements, std::vector<Vector>{});
    }

    const auto& r_process_info    = rModelPart.GetProcessInfo();
    const auto  it_elements_begin = rModelPart.ElementsBegin();

    return IndexPartition<std::size_t>(number_of_elements).for_each<PointCountsReduction>([&](std::size_t Index) {
        auto& r_element = *(it_elements_begin + Index);
        return r_element.IsActive()
                   ? EvaluateActiveElement(r_element, mPreviousElasticPredictorDeviations[Index], r_process_info)
                   : CountInactiveElement(r_element, r_process_info);
    });
}

LocalErrorEvaluator::PointCounts LocalErrorEvaluator::EvaluateActiveElement(Element& rElement,
                                                                            std::vector<Vector>& rPreviousElasticPredictorDeviations,
                                                                            const ProcessInfo& rProcessInfo) const
{
    auto result = PointCounts{};

    const auto constitutive_laws = GetConstitutiveLaws(rElement, rProcessInfo);
    rPreviousElasticPredictorDeviations.resize(constitutive_laws.size());

    auto parameters = ConstitutiveLaw::Parameters{rElement.GetGeometry(), rElement.GetProperties(), rProcessInfo};
    for (std::size_t i = 0; i < constitutive_laws.size(); ++i) {
        auto* p_provider = GetLocalErrorDataProvider(constitutive_laws[i]);
        if (!p_provider) continue;

        const auto stress_point_type = p_provider->GetStressPointType();
        if (!stress_point_type) continue;

        const auto data = p_provider->CalculateLocalErrorData(parameters);
        AddStressPoint(*stress_point_type, data, rPreviousElasticPredictorDeviations[i], result);
        rPreviousElasticPredictorDeviations[i] = data.ElasticPredictorDeviation;
    }

    return result;
}

LocalErrorEvaluator::PointCounts LocalErrorEvaluator::CountInactiveElement(Element& rElement, const ProcessInfo& rProcessInfo)
{
    // The stress points of inactive soil elements are counted as elastic points. They are not
    // evaluated, since their material state is not updated.
    auto result = PointCounts{};
    for (const auto& rp_constitutive_law : GetConstitutiveLaws(rElement, rProcessInfo)) {
        const auto* p_provider = GetLocalErrorDataProvider(rp_constitutive_law);
        if (p_provider && p_provider->GetStressPointType() == Geo::StressPointType::Soil) {
            ++result.SoilElastic;
        }
    }
    return result;
}

void LocalErrorEvaluator::AddStressPoint(Geo::StressPointType      Type,
                                         const Geo::LocalErrorData& rData,
                                         const Vector&              rPreviousElasticPredictorDeviation,
                                         PointCounts&               rCounts) const
{
    const auto local_error_for = [&rData, &rPreviousElasticPredictorDeviation](
                                     double MinimumReferenceStress, std::size_t IndexOfFirstComponent = 0) {
        return LocalErrorUtilities::CalculateLocalError(
            rData.ElasticPredictorDeviation, rPreviousElasticPredictorDeviation, rData.MaximumShearStress,
            rData.Cohesion, MinimumReferenceStress, IndexOfFirstComponent);
    };

    if (Type == Geo::StressPointType::Interface) {
        // Only plastic interface points are monitored, based on their shear tractions (section 9.1.2.3)
        if (!rData.IsPlastic) return;

        ++rCounts.InterfacePlastic;
        if (local_error_for(mSettings.MinimumReferenceStress, rData.IndexOfFirstShearComponent) >
            mSettings.ToleratedInterfacePlasticLocalError) {
            ++rCounts.InaccurateInterfacePlastic;
        }
        return;
    }

    if (rData.IsPlastic) {
        // Section 9.1.2.1
        ++rCounts.SoilPlastic;
        if (local_error_for(mSettings.MinimumReferenceStress) > mSettings.ToleratedSoilPlasticLocalError) {
            ++rCounts.InaccurateSoilPlastic;
        }
        return;
    }

    // Section 9.1.2.2: all elastic points count towards the total number of elastic points, but only
    // non-linear elastic points (i.e. those with a stress-dependent stiffness) can be inaccurate
    ++rCounts.SoilElastic;
    constexpr auto reference_pressure_divisor = 200.0;
    if (rData.HasStressDependentStiffness &&
        local_error_for(mSettings.ReferencePressure / reference_pressure_divisor) >
            mSettings.ToleratedSoilElasticLocalError) {
        ++rCounts.InaccurateSoilElastic;
    }
}

double LocalErrorEvaluator::GetToleratedNumberOfInaccuratePoints(std::size_t NumberOfPoints,
                                                                 double ToleratedInaccuratePointFraction) const
{
    return ToleratedInaccuratePointFraction * static_cast<double>(NumberOfPoints) +
           static_cast<double>(mSettings.NumberOfAlwaysToleratedInaccuratePoints);
}

bool LocalErrorEvaluator::IsConverged(const PointCounts& rCounts) const
{
    const auto is_accurate_enough = [this](std::size_t NumberOfInaccuratePoints, std::size_t NumberOfPoints,
                                           double ToleratedInaccuratePointFraction) {
        return static_cast<double>(NumberOfInaccuratePoints) <
               GetToleratedNumberOfInaccuratePoints(NumberOfPoints, ToleratedInaccuratePointFraction);
    };

    return is_accurate_enough(rCounts.InaccurateSoilPlastic, rCounts.SoilPlastic,
                              mSettings.ToleratedInaccurateSoilPlasticPointFraction) &&
           is_accurate_enough(rCounts.InaccurateInterfacePlastic, rCounts.InterfacePlastic,
                              mSettings.ToleratedInaccurateInterfacePlasticPointFraction) &&
           is_accurate_enough(rCounts.InaccurateSoilElastic, rCounts.SoilElastic,
                              mSettings.ToleratedInaccurateSoilElasticPointFraction);
}

void LocalErrorEvaluator::PrintCounts(std::ostream& rOStream, const PointCounts& rCounts) const
{
    const auto print_line = [this, &rOStream](const std::string& rLabel, std::size_t NumberOfPoints,
                                              std::size_t NumberOfInaccuratePoints, double Fraction) {
        rOStream << rLabel << NumberOfPoints << ", inaccurate: " << NumberOfInaccuratePoints
                 << " (tolerated: < " << GetToleratedNumberOfInaccuratePoints(NumberOfPoints, Fraction) << ")\n";
    };

    print_line("Soil plastic points:      ", rCounts.SoilPlastic, rCounts.InaccurateSoilPlastic,
               mSettings.ToleratedInaccurateSoilPlasticPointFraction);
    print_line("Soil elastic points:      ", rCounts.SoilElastic, rCounts.InaccurateSoilElastic,
               mSettings.ToleratedInaccurateSoilElasticPointFraction);
    print_line("Interface plastic points: ", rCounts.InterfacePlastic, rCounts.InaccurateInterfacePlastic,
               mSettings.ToleratedInaccurateInterfacePlasticPointFraction);
}

} // namespace Kratos

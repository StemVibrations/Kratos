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

#include "custom_constitutive/local_error_data_provider.h"
#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <cstddef>
#include <iosfwd>
#include <vector>

namespace Kratos
{

class Element;
class ModelPart;
class ProcessInfo;

// Evaluates the local error criteria of section 9.1.2 of the PLAXIS Scientific Manual:
//  - inaccurate plastic points for soil elements (9.1.2.1),
//  - non-linear inaccurate elastic points for soil elements (9.1.2.2),
//  - inaccurate plastic points for interfaces (9.1.2.3).
//
// A stress point j is inaccurate when || sigma_eq,j - sigma_c,j || / max(tau_max,j, c, sigma_ref)
// exceeds the tolerated local error, where sigma_eq,j = sigma_c,j-1 + D^e delta_eps_j is the
// equilibrium stress, sigma_c,j is the constitutive stress and sigma_ref is either the minimum
// reference stress (plastic points) or 1/200 of the reference pressure (elastic points).
// The equilibrium stress follows from the deviation P = sigma_0 + D^e (eps - eps_0) - sigma_c of the
// constitutive stress from the elastic predictor of the step: sigma_eq,j - sigma_c,j = P_j - P_j-1.
// Therefore, only the deviations of the previous iteration need to be stored.
class KRATOS_API(GEO_MECHANICS_APPLICATION) LocalErrorEvaluator
{
public:
    struct Settings {
        double      ToleratedSoilPlasticLocalError                   = 0.01;
        double      ToleratedSoilElasticLocalError                   = 0.01;
        double      ToleratedInterfacePlasticLocalError              = 0.01;
        double      ToleratedInaccurateSoilPlasticPointFraction      = 0.1;
        double      ToleratedInaccurateSoilElasticPointFraction      = 0.1;
        double      ToleratedInaccurateInterfacePlasticPointFraction = 0.1;
        std::size_t NumberOfAlwaysToleratedInaccuratePoints          = 3;
        double      MinimumReferenceStress                           = 1.0e3; // 1 kPa, in Pa
        double      ReferencePressure                                = 1.0e5; // 100 kPa, in Pa
    };

    struct PointCounts {
        std::size_t SoilPlastic                = 0;
        std::size_t InaccurateSoilPlastic      = 0;
        std::size_t SoilElastic                = 0;
        std::size_t InaccurateSoilElastic      = 0;
        std::size_t InterfacePlastic           = 0;
        std::size_t InaccurateInterfacePlastic = 0;

        PointCounts& operator+=(const PointCounts& rOther);
        bool         operator==(const PointCounts& rOther) const = default;
    };

    LocalErrorEvaluator() = default;
    explicit LocalErrorEvaluator(const Settings& rSettings);

    [[nodiscard]] const Settings& GetSettings() const;

    // Forgets the state of the previous iteration, such that the next evaluation is carried out
    // with respect to the state at the start of the step (i.e. sigma_c,0 = sigma_0).
    void ResetPreviousIterationState();

    // Counts the (inaccurate) stress points of all elements of the model part, based on the most
    // recently calculated material responses. Afterwards, the current state is stored as the
    // state of the previous iteration.
    PointCounts EvaluateAndStoreIterationState(ModelPart& rModelPart);

    // Adds a single stress point to the given counts
    void AddStressPoint(Geo::StressPointType      Type,
                        const Geo::LocalErrorData& rData,
                        const Vector&              rPreviousElasticPredictorDeviation,
                        PointCounts&               rCounts) const;

    [[nodiscard]] double GetToleratedNumberOfInaccuratePoints(std::size_t NumberOfPoints,
                                                              double ToleratedInaccuratePointFraction) const;
    [[nodiscard]] bool   IsConverged(const PointCounts& rCounts) const;
    void                 PrintCounts(std::ostream& rOStream, const PointCounts& rCounts) const;

private:
    Settings mSettings;

    // The elastic predictor deviations of the previous iteration, per element (in the order of the
    // elements of the model part) and per integration point. An empty vector means zero deviation.
    std::vector<std::vector<Vector>> mPreviousElasticPredictorDeviations;

    PointCounts EvaluateActiveElement(Element&           rElement,
                                      std::vector<Vector>& rPreviousElasticPredictorDeviations,
                                      const ProcessInfo&   rProcessInfo) const;
    static PointCounts CountInactiveElement(Element& rElement, const ProcessInfo& rProcessInfo);
};

} // namespace Kratos

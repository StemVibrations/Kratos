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

#include "includes/constitutive_law.h"
#include "includes/kratos_export_api.h"
#include "includes/ublas_interface.h"

#include <optional>

namespace Kratos
{

namespace Geo
{

enum class StressPointType { Soil, Interface };

// The state of a single stress point that is required to evaluate its local error indicator
struct LocalErrorData {
    // Whether plastic flow occurred in the current step, i.e. whether the constitutive stress
    // deviates from the elastic predictor of the step
    bool IsPlastic = false;

    // Whether the elastic stiffness depends on the stress state (e.g. hardening soil type models)
    bool HasStressDependentStiffness = false;

    // The deviation of the constitutive stress from the elastic predictor of the current step:
    //     sigma_0 + D^e (eps - eps_0) - sigma_c
    // where sigma_0 and eps_0 are the stress and strain at the start of the step, D^e is the
    // elastic stiffness matrix and sigma_c is the constitutive stress. Since sigma_0 and D^e are
    // constant within a step, the change of this quantity between two subsequent iterations
    // j-1 and j equals sigma_eq,j - sigma_c,j, where sigma_eq,j = sigma_c,j-1 + D^e delta_eps_j is
    // the equilibrium stress.
    Vector ElasticPredictorDeviation;

    // The index of the first shear component of the stress vector. For interfaces, only the shear
    // components (i.e. the components from this index onwards) are taken into account by the local
    // error. Zero means all components.
    std::size_t IndexOfFirstShearComponent = 0;

    // The maximum shear stress of the constitutive stress (soil: the radius of the largest Mohr
    // circle, interface: the magnitude of the shear traction)
    double MaximumShearStress = 0.0;

    double Cohesion = 0.0;

    // The strain energy increments (per unit volume) of the current step, which define the current
    // stiffness parameter of the global error criterion:
    //     total:   delta_eps . delta_sigma
    //     elastic: delta_eps . D^e delta_eps
    // where delta_eps and delta_sigma are the strain and stress increments with respect to the start
    // of the step. Both are zero when the law does not provide them.
    double TotalStrainEnergyIncrement   = 0.0;
    double ElasticStrainEnergyIncrement = 0.0;
};

} // namespace Geo

// Interface of constitutive laws that can provide the data needed by the local error criteria
class KRATOS_API(GEO_MECHANICS_APPLICATION) LocalErrorDataProvider
{
public:
    virtual ~LocalErrorDataProvider();

    // Returns the type of stress point that is represented by the constitutive law. An empty
    // optional means that the local error criteria don't apply (e.g. for structural elements).
    [[nodiscard]] virtual std::optional<Geo::StressPointType> GetStressPointType() const = 0;

    // Returns the local error data that belongs to the most recently calculated material
    // response. Calling this function must not alter the material state of the law.
    [[nodiscard]] virtual Geo::LocalErrorData CalculateLocalErrorData(ConstitutiveLaw::Parameters& rParameters) = 0;

protected:
    LocalErrorDataProvider()                                             = default;
    LocalErrorDataProvider(const LocalErrorDataProvider&)                = default;
    LocalErrorDataProvider& operator=(const LocalErrorDataProvider&)     = default;
    LocalErrorDataProvider(LocalErrorDataProvider&&) noexcept            = default;
    LocalErrorDataProvider& operator=(LocalErrorDataProvider&&) noexcept = default;
};

} // namespace Kratos

// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Mohamed Nabi,
//                   Wijtze Pieter Kikstra
//

// Application includes
#include "custom_constitutive/interface_coulomb_with_tension_cut_off.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/sigma_tau.hpp"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

#include <algorithm>
#include <cmath>

namespace
{
using namespace Kratos;

// Estimates how many sub-steps are needed to integrate the plastic remainder accurately.
// It performs a single "probe" return mapping of the elastic predictor to measure how far
// that predictor overshoots the yield surface, and subdivides such that each sub-step
// overshoots by at most a target fraction of the traction magnitude. The kappa (hardening)
// update caused by the probe is rolled back. The result is clamped to [1, MaxNumberOfSubSteps].
std::size_t CalculateAdaptiveNumberOfSubSteps(CoulombWithTensionCutOffImpl& rImpl,
                                              const Geo::SigmaTau&          rTrialTraction,
                                              const Matrix&                 rElasticMatrix,
                                              std::size_t                   MaxNumberOfSubSteps)
{
    rImpl.SaveKappaOfCoulombYieldSurface();
    const auto mapped_traction = rImpl.DoReturnMapping(
        rTrialTraction, rElasticMatrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
    rImpl.RestoreKappaOfCoulombYieldSurface();

    const Vector trial_values  = rTrialTraction.CopyTo<Vector>();
    const Vector mapped_values = mapped_traction.CopyTo<Vector>();

    const auto overshoot          = norm_2(trial_values - mapped_values);
    const auto stress_scale       = std::max(norm_2(trial_values), 1.0e-12);
    const auto relative_overshoot = overshoot / stress_scale;

    // Each sub-step's elastic predictor is allowed to overshoot the yield surface by at most
    // this fraction of the traction magnitude before we subdivide further.
    constexpr auto target_relative_overshoot_per_sub_step = 0.1;

    const auto number_of_sub_steps = static_cast<std::size_t>(
        std::ceil(relative_overshoot / target_relative_overshoot_per_sub_step));

    return std::clamp(number_of_sub_steps, std::size_t{1}, MaxNumberOfSubSteps);
}

// Pegasus algorithm (an accelerated regula-falsi / false-position method with guaranteed
// convergence). It locates the factor alpha in [0, 1] at which the elastic traction path
//     traction(alpha) = traction_start + alpha * elastic_traction_increment
// first touches the yield surface, i.e. where the combined yield function F(alpha) == 0.
// The yield functions are expressed in terms of the shear magnitude, so |tau| is used.
double CalculateYieldSurfaceIntersectionFactor(CoulombWithTensionCutOffImpl& rImpl,
                                               const Geo::SigmaTau& rStartTraction,
                                               const Geo::SigmaTau& rElasticTractionIncrement)
{
    // Hard-coded control parameters (no extra material input required).
    constexpr auto        yield_function_tolerance = 1.0e-10;
    constexpr std::size_t max_number_of_iterations = 100;

    const auto yield_function_value_at = [&](double Alpha) {
        auto traction =
            Geo::SigmaTau{rStartTraction.Values() + Alpha * rElasticTractionIncrement.Values()};
        traction.Tau() = std::abs(traction.Tau());
        return rImpl.YieldFunctionValue(traction);
    };

    auto alpha_0 = 0.0;
    auto alpha_1 = 1.0;
    auto f_0     = yield_function_value_at(alpha_0);
    auto f_1     = yield_function_value_at(alpha_1);

    // If the bracket is not valid (e.g. the start state is already on/outside the surface),
    // there is no elastic portion to skip.
    if (f_0 >= 0.0) return 0.0;
    if (f_1 <= 0.0) return 1.0;

    auto alpha = alpha_1;
    for (std::size_t iteration = 0; iteration < max_number_of_iterations; ++iteration) {
        alpha            = alpha_1 - f_1 * (alpha_1 - alpha_0) / (f_1 - f_0);
        const auto f_new = yield_function_value_at(alpha);
        if (std::abs(f_new) < yield_function_tolerance) break;

        if (f_new * f_1 < 0.0) {
            alpha_0 = alpha_1;
            f_0     = f_1;
        } else {
            // Pegasus acceleration to avoid one-sided stagnation of classic regula-falsi.
            f_0 *= f_1 / (f_1 + f_new);
        }
        alpha_1 = alpha;
        f_1     = f_new;
    }

    return std::clamp(alpha, 0.0, 1.0);
}

} // namespace

namespace Kratos
{

InterfaceCoulombWithTensionCutOff::InterfaceCoulombWithTensionCutOff(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveDimension)
    : mpConstitutiveDimension(std::move(pConstitutiveDimension)),
      mTractionVector(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mTractionVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize())),
      mRelativeDisplacementVectorFinalized(ZeroVector(mpConstitutiveDimension->GetStrainSize()))
{
}

ConstitutiveLaw::Pointer InterfaceCoulombWithTensionCutOff::Clone() const
{
    auto p_result = std::make_shared<InterfaceCoulombWithTensionCutOff>(mpConstitutiveDimension->Clone());
    p_result->mTractionVector                      = mTractionVector;
    p_result->mTractionVectorFinalized             = mTractionVectorFinalized;
    p_result->mRelativeDisplacementVectorFinalized = mRelativeDisplacementVectorFinalized;
    p_result->mCoulombWithTensionCutOffImpl        = mCoulombWithTensionCutOffImpl;
    p_result->mIsModelInitialized                  = mIsModelInitialized;
    return p_result;
}

Vector& InterfaceCoulombWithTensionCutOff::GetValue(const Variable<Vector>& rVariable, Vector& rValue)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        rValue = mTractionVector;
    } else {
        rValue = ConstitutiveLaw::GetValue(rVariable, rValue);
    }
    return rValue;
}

int& InterfaceCoulombWithTensionCutOff::GetValue(const Variable<int>& rVariable, int& rValue)
{
    if (rVariable == GEO_PLASTICITY_STATUS) {
        rValue = static_cast<int>(mCoulombWithTensionCutOffImpl.GetPlasticityStatus());
    }
    return rValue;
}

void InterfaceCoulombWithTensionCutOff::SetValue(const Variable<Vector>& rVariable,
                                                 const Vector&           rValue,
                                                 const ProcessInfo&      rCurrentProcessInfo)
{
    if (rVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        mTractionVector = rValue;
    } else {
        KRATOS_ERROR << "Can't set value of " << rVariable.Name() << ": unsupported variable\n";
    }
}

SizeType InterfaceCoulombWithTensionCutOff::WorkingSpaceDimension()
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return N_DIM_2D;
}

int InterfaceCoulombWithTensionCutOff::Check(const Properties&   rMaterialProperties,
                                             const GeometryType& rElementGeometry,
                                             const ProcessInfo&  rCurrentProcessInfo) const
{
    const auto result = ConstitutiveLaw::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "property", CheckProperties::Bounds::AllInclusive);
    check_properties.Check(GEO_COHESION);
    constexpr auto max_value_angle = 90.0;
    check_properties.SingleUseBounds(CheckProperties::Bounds::AllExclusive).Check(GEO_FRICTION_ANGLE, 0.0, max_value_angle);
    check_properties.Check(GEO_DILATANCY_ANGLE, rMaterialProperties[GEO_FRICTION_ANGLE]);
   /* check_properties.Check(
        GEO_TENSILE_STRENGTH,
        rMaterialProperties[GEO_COHESION] /
            std::tan(MathUtils<>::DegreesToRadians(rMaterialProperties[GEO_FRICTION_ANGLE])));*/
    check_properties.Check(INTERFACE_NORMAL_STIFFNESS);
    check_properties.Check(INTERFACE_SHEAR_STIFFNESS);
    return result;
}

ConstitutiveLaw::StressMeasure InterfaceCoulombWithTensionCutOff::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

SizeType InterfaceCoulombWithTensionCutOff::GetStrainSize() const
{
    // Note that this implementation assumes line interface elements. It needs to be modified when planar interface elements become available.
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StrainMeasure InterfaceCoulombWithTensionCutOff::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

bool InterfaceCoulombWithTensionCutOff::IsIncremental() { return true; }

bool InterfaceCoulombWithTensionCutOff::RequiresInitializeMaterialResponse() { return true; }

void InterfaceCoulombWithTensionCutOff::InitializeMaterial(const Properties& rMaterialProperties,
                                                           const Geometry<Node>&,
                                                           const Vector&)
{
    mCoulombWithTensionCutOffImpl = CoulombWithTensionCutOffImpl{rMaterialProperties};

    mRelativeDisplacementVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mTractionVectorFinalized =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void InterfaceCoulombWithTensionCutOff::InitializeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!mIsModelInitialized) {
        mTractionVectorFinalized             = rConstitutiveLawParameters.GetStressVector();
        mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
        mIsModelInitialized                  = true;
    }
}

void InterfaceCoulombWithTensionCutOff::CalculateMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    if (!rConstitutiveLawParameters.GetOptions().Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        return;
    }

    const auto& r_properties   = rConstitutiveLawParameters.GetMaterialProperties();
    const auto  elastic_matrix = mpConstitutiveDimension->CalculateElasticMatrix(r_properties);

    const Vector& r_relative_displacement_vector = rConstitutiveLawParameters.GetStrainVector();

    // Full elastic predictor over the entire relative-displacement increment.
    const auto full_trial_sigma_tau = CalculateTrialTractionVector(
        r_relative_displacement_vector, r_properties[INTERFACE_NORMAL_STIFFNESS],
        r_properties[INTERFACE_SHEAR_STIFFNESS]);

    // Admissibility (and the yield functions) are evaluated on the shear magnitude.
    auto full_trial_for_check  = full_trial_sigma_tau;
    full_trial_for_check.Tau() = std::abs(full_trial_for_check.Tau());

    // If the whole step stays elastic, there is nothing to integrate and no need to sub-step.
    if (mCoulombWithTensionCutOffImpl.IsAdmissibleStressState(full_trial_for_check)) {
        mTractionVector                              = full_trial_sigma_tau.CopyTo<Vector>();
        rConstitutiveLawParameters.GetStressVector() = mTractionVector;
        return;
    }

    // ----- Pegasus yield-surface intersection -----
    // Locate the fraction of the increment that remains purely elastic before the traction
    // path first touches the yield surface. Only the remaining (plastic) part is integrated
    // with the return mapping, which improves both accuracy and robustness.
    const auto finalized_traction = Geo::SigmaTau{mTractionVectorFinalized};
    const auto elastic_traction_increment =
        Geo::SigmaTau{full_trial_sigma_tau.Values() - finalized_traction.Values()};
    const double intersection_factor = CalculateYieldSurfaceIntersectionFactor(
        mCoulombWithTensionCutOffImpl, finalized_traction, elastic_traction_increment);

    const Vector total_relative_displacement_increment =
        r_relative_displacement_vector - mRelativeDisplacementVectorFinalized;

    // Advance elastically up to the intersection point; this becomes the starting (committed)
    // state for the sub-stepped plastic integration.
    const auto plastic_start_traction = Geo::SigmaTau{
        finalized_traction.Values() + intersection_factor * elastic_traction_increment.Values()};
    const Vector plastic_start_relative_displacement =
        mRelativeDisplacementVectorFinalized + intersection_factor * total_relative_displacement_increment;
    const Vector plastic_relative_displacement_increment =
        (1.0 - intersection_factor) * total_relative_displacement_increment;

    mTractionVector = plastic_start_traction.CopyTo<Vector>();

    // ----- adaptive sub-stepping (with an upper bound) of the plastic remainder -----
    constexpr std::size_t max_number_of_sub_steps = 5000;
    const std::size_t     number_of_sub_steps     = CalculateAdaptiveNumberOfSubSteps(
        mCoulombWithTensionCutOffImpl, full_trial_for_check, elastic_matrix, max_number_of_sub_steps);

    // Running committed state for the sub-stepping (start from the Pegasus intersection point)
    Vector committed_traction              = mTractionVector;
    Vector committed_relative_displacement = plastic_start_relative_displacement;

    for (std::size_t sub = 1; sub <= number_of_sub_steps; ++sub) {
        // relative displacement at the end of this sub-step
        const Vector sub_relative_displacement =
            plastic_start_relative_displacement +
            (static_cast<double>(sub) / static_cast<double>(number_of_sub_steps)) *
                plastic_relative_displacement_increment;

        // elastic predictor from the *committed* sub-step state
        auto trial_sigma_tau = Geo::SigmaTau{
            committed_traction +
            prod(elastic_matrix, sub_relative_displacement - committed_relative_displacement)};

        const auto negative   = std::signbit(trial_sigma_tau.Tau());
        trial_sigma_tau.Tau() = std::abs(trial_sigma_tau.Tau());

        auto mapped_sigma_tau = trial_sigma_tau;
        if (!mCoulombWithTensionCutOffImpl.IsAdmissibleStressState(trial_sigma_tau)) {
            mapped_sigma_tau = mCoulombWithTensionCutOffImpl.DoReturnMapping(
                trial_sigma_tau, elastic_matrix, Geo::PrincipalStresses::AveragingType::NO_AVERAGING);
        }
        if (negative) mapped_sigma_tau.Tau() *= -1.0;

        // commit this sub-step (kappa is already updated inside DoReturnMapping)
        mTractionVector                 = mapped_sigma_tau.CopyTo<Vector>();
        committed_traction              = mTractionVector;
        committed_relative_displacement = sub_relative_displacement;
    }

    rConstitutiveLawParameters.GetStressVector() = mTractionVector;
}

Geo::SigmaTau InterfaceCoulombWithTensionCutOff::CalculateTrialTractionVector(const Vector& rRelativeDisplacementVector,
                                                                              double NormalStiffness,
                                                                              double ShearStiffness) const
{
    constexpr auto number_of_normal_components = std::size_t{1};
    return Geo::SigmaTau{mTractionVectorFinalized +
                         prod(ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
                                  NormalStiffness, ShearStiffness, GetStrainSize(), number_of_normal_components),
                              rRelativeDisplacementVector - mRelativeDisplacementVectorFinalized)};
}

void InterfaceCoulombWithTensionCutOff::FinalizeMaterialResponseCauchy(Parameters& rConstitutiveLawParameters)
{
    mRelativeDisplacementVectorFinalized = rConstitutiveLawParameters.GetStrainVector();
    mTractionVectorFinalized             = mTractionVector;
}

Matrix& InterfaceCoulombWithTensionCutOff::CalculateValue(Parameters& rConstitutiveLawParameters,
                                                          const Variable<Matrix>& rVariable,
                                                          Matrix&                 rValue)
{
    if (rVariable == CONSTITUTIVE_MATRIX) {
        const auto&    r_properties = rConstitutiveLawParameters.GetMaterialProperties();
        constexpr auto number_of_normal_components = std::size_t{1};
        rValue = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
            r_properties[INTERFACE_NORMAL_STIFFNESS], r_properties[INTERFACE_SHEAR_STIFFNESS],
            GetStrainSize(), number_of_normal_components);
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

void InterfaceCoulombWithTensionCutOff::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.save("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.save("TractionVector", mTractionVector);
    rSerializer.save("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.save("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.save("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.save("IsModelInitialized", mIsModelInitialized);
}

void InterfaceCoulombWithTensionCutOff::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    rSerializer.load("ConstitutiveDimension", mpConstitutiveDimension);
    rSerializer.load("TractionVector", mTractionVector);
    rSerializer.load("TractionVectorFinalized", mTractionVectorFinalized);
    rSerializer.load("RelativeDisplacementVectorFinalized", mRelativeDisplacementVectorFinalized);
    rSerializer.load("CoulombWithTensionCutOffImpl", mCoulombWithTensionCutOffImpl);
    rSerializer.load("IsModelInitialized", mIsModelInitialized);
}

} // Namespace Kratos
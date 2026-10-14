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

#include "custom_strategies/convergence_criteria/local_error_evaluator.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

#include <sstream>
#include <string>

namespace Kratos
{

/**
 * @class GeoLocalErrorCriteria
 * @ingroup GeoMechanicsApplication
 * @brief Convergence criterion based on the local error indicators of the stress points
 * @details Implements the local error criteria:
 * the number of inaccurate plastic soil points, inaccurate non-linear elastic soil points and
 * inaccurate plastic interface points must each remain below a tolerated number. This criterion
 * does not check global equilibrium, so it is meant to be combined with a global criterion (e.g.
 * by means of an AndCriteria). The constitutive laws need to implement LocalErrorDataProvider; the
 * stress points of other laws are ignored.
 * The evaluation requires the material responses to be up to date with the latest solution, hence
 * this criterion requests the right-hand side to be rebuilt before calling PostCriteria.
 */
template <class TSparseSpace, class TDenseSpace>
class GeoLocalErrorCriteria : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoLocalErrorCriteria);

    using BaseType          = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using ClassType         = GeoLocalErrorCriteria<TSparseSpace, TDenseSpace>;
    using DofsArrayType     = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    GeoLocalErrorCriteria() : GeoLocalErrorCriteria(Parameters{"{}"}) {}

    explicit GeoLocalErrorCriteria(Parameters ThisParameters) : BaseType()
    {
        ThisParameters = this->ValidateAndAssignParameters(
            UseToleratedErrorAsDefaultForLocalErrors(ThisParameters), this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
        this->mActualizeRHSIsNeeded = true;
    }

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    void InitializeSolutionStep(ModelPart&               rModelPart,
                                DofsArrayType&           rDofSet,
                                const TSystemMatrixType& rA,
                                const TSystemVectorType& rDx,
                                const TSystemVectorType& rb) override
    {
        BaseType::InitializeSolutionStep(rModelPart, rDofSet, rA, rDx, rb);

        // The first iteration of a step is evaluated with respect to the state at the start of the step
        mEvaluator.ResetPreviousIterationState();
    }

    bool PostCriteria(ModelPart& rModelPart,
                      DofsArrayType&,
                      const TSystemMatrixType&,
                      const TSystemVectorType&,
                      const TSystemVectorType&) override
    {
        const auto counts       = mEvaluator.EvaluateAndStoreIterationState(rModelPart);
        const auto is_converged = mEvaluator.IsConverged(counts);

        if (this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0) {
            std::ostringstream message;
            mEvaluator.PrintCounts(message, counts);
            message << "Local error criteria " << (is_converged ? "satisfied" : "not satisfied");
            KRATOS_INFO("GeoLocalErrorCriteria") << message.str() << std::endl;
        }

        return is_converged;
    }

    [[nodiscard]] const LocalErrorEvaluator& GetEvaluator() const { return mEvaluator; }

    [[nodiscard]] Parameters GetDefaultParameters() const override
    {
        auto default_parameters = Parameters(R"(
        {
            "name"                                                 : "geo_local_error_criteria",
            "tolerated_error"                                      : 0.01,
            "tolerated_soil_plastic_local_error"                   : 0.01,
            "tolerated_soil_elastic_local_error"                   : 0.01,
            "tolerated_interface_plastic_local_error"              : 0.01,
            "tolerated_inaccurate_soil_plastic_point_fraction"     : 0.1,
            "tolerated_inaccurate_soil_elastic_point_fraction"     : 0.1,
            "tolerated_inaccurate_interface_plastic_point_fraction": 0.1,
            "number_of_always_tolerated_inaccurate_points"         : 3,
            "minimum_reference_stress"                             : 1.0e3,
            "reference_pressure"                                   : 1.0e5
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());
        return default_parameters;
    }

    static std::string Name() { return "geo_local_error_criteria"; }

    [[nodiscard]] std::string Info() const override { return "GeoLocalErrorCriteria"; }

protected:
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        KRATOS_ERROR_IF(ThisParameters["number_of_always_tolerated_inaccurate_points"].GetInt() < 0)
            << "The number_of_always_tolerated_inaccurate_points must not be negative" << std::endl;

        auto settings                                = LocalErrorEvaluator::Settings{};
        settings.ToleratedSoilPlasticLocalError      = ThisParameters["tolerated_soil_plastic_local_error"].GetDouble();
        settings.ToleratedSoilElasticLocalError      = ThisParameters["tolerated_soil_elastic_local_error"].GetDouble();
        settings.ToleratedInterfacePlasticLocalError = ThisParameters["tolerated_interface_plastic_local_error"].GetDouble();
        settings.ToleratedInaccurateSoilPlasticPointFraction =
            ThisParameters["tolerated_inaccurate_soil_plastic_point_fraction"].GetDouble();
        settings.ToleratedInaccurateSoilElasticPointFraction =
            ThisParameters["tolerated_inaccurate_soil_elastic_point_fraction"].GetDouble();
        settings.ToleratedInaccurateInterfacePlasticPointFraction =
            ThisParameters["tolerated_inaccurate_interface_plastic_point_fraction"].GetDouble();
        settings.NumberOfAlwaysToleratedInaccuratePoints = static_cast<std::size_t>(
            ThisParameters["number_of_always_tolerated_inaccurate_points"].GetInt());
        settings.MinimumReferenceStress = ThisParameters["minimum_reference_stress"].GetDouble();
        settings.ReferencePressure      = ThisParameters["reference_pressure"].GetDouble();
        mEvaluator                      = LocalErrorEvaluator{settings};
    }

private:
    LocalErrorEvaluator mEvaluator;

    // By default, the tolerated local errors are equal to the tolerated error
    static Parameters UseToleratedErrorAsDefaultForLocalErrors(Parameters ThisParameters)
    {
        if (!ThisParameters.Has("tolerated_error")) return ThisParameters;

        const auto tolerated_error = ThisParameters["tolerated_error"].GetDouble();
        for (const auto& r_name :
             {"tolerated_soil_plastic_local_error", "tolerated_soil_elastic_local_error",
              "tolerated_interface_plastic_local_error"}) {
            if (!ThisParameters.Has(r_name)) ThisParameters.AddDouble(r_name, tolerated_error);
        }
        return ThisParameters;
    }
}; // Class GeoLocalErrorCriteria

} // namespace Kratos

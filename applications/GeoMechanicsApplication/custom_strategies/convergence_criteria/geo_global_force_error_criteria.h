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

#include "custom_strategies/convergence_criteria/current_stiffness_parameter_calculator.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace Kratos
{

/**
 * @class GeoGlobalForceErrorCriteria
 * @ingroup GeoMechanicsApplication
 * @brief Convergence criterion based on the global error of the nodal forces
 *
 *     GlobalErrorForce = || r || / (CSP || f_ext^inact || + || f_int^act ||) < ForceRelativeTolerance
 *
 * where r is the out-of-balance force vector, f_ext^inact are the external forces that were already
 * in equilibrium at the start of the stage, f_int^act are the internal forces that have developed in
 * the stage so far and CSP is the current stiffness parameter (see CurrentStiffnessParameterCalculator).
 * Only the free displacement degrees of freedom are taken into account, i.e. the residuals of other
 * degrees of freedom (e.g. water pressures and rotations) need to be checked by another criterion.
 * The residual is taken from the strategy, which also needs to provide the reference forces (by means
 * of SetReferenceForces) before each call of PostCriteria.
 */
template <class TSparseSpace, class TDenseSpace>
class GeoGlobalForceErrorCriteria : public ConvergenceCriteria<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoGlobalForceErrorCriteria);

    using BaseType          = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using ClassType         = GeoGlobalForceErrorCriteria<TSparseSpace, TDenseSpace>;
    using DofsArrayType     = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    GeoGlobalForceErrorCriteria() : GeoGlobalForceErrorCriteria(Parameters{"{}"}) {}

    explicit GeoGlobalForceErrorCriteria(Parameters ThisParameters) : BaseType()
    {
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

    typename BaseType::Pointer Create(Parameters ThisParameters) const override
    {
        return Kratos::make_shared<ClassType>(ThisParameters);
    }

    // Sets the forces that the residual is compared with, expressed in the same space as the residual:
    //  - the inactive external forces f_ext^inact, i.e. the external forces that were already in
    //    equilibrium at the start of the stage,
    //  - the active internal forces f_int^act, i.e. the internal forces that have developed in the
    //    stage so far.
    // Both need to be set again before each evaluation of the criterion.
    void SetReferenceForces(const TSystemVectorType& rInactiveExternalForces, const TSystemVectorType& rActiveInternalForces)
    {
        mInactiveExternalForces = rInactiveExternalForces;
        mActiveInternalForces   = rActiveInternalForces;
        mHasReferenceForces     = true;
    }

    bool PostCriteria(ModelPart&     rModelPart,
                      DofsArrayType& rDofSet,
                      const TSystemMatrixType&,
                      const TSystemVectorType&,
                      const TSystemVectorType& rb) override
    {
        KRATOS_ERROR_IF_NOT(mHasReferenceForces)
            << "The reference forces of the " << Info() << " have not been set. This criterion "
            << "requires a strategy that provides them (e.g. GeoMechanicsNewtonRaphsonStrategy)." << std::endl;
        mHasReferenceForces = false;

        const auto equation_ids = GetEquationIdsOfFreeDisplacementDofs(rDofSet, TSparseSpace::Size(rb));
        const auto residual_norm                = CalculateNorm(rb, equation_ids);
        const auto inactive_external_force_norm = CalculateNorm(mInactiveExternalForces, equation_ids);
        const auto active_internal_force_norm   = CalculateNorm(mActiveInternalForces, equation_ids);
        const auto current_stiffness_parameter  = CurrentStiffnessParameterCalculator::Calculate(rModelPart);
        const auto global_error = CalculateGlobalForceError(residual_norm, current_stiffness_parameter,
                                                            inactive_external_force_norm, active_internal_force_norm);

        const auto is_converged = global_error < mForceRelativeTolerance || residual_norm <= mForceAbsoluteTolerance;

        KRATOS_INFO_IF("GeoGlobalForceErrorCriteria",
                       this->GetEchoLevel() > 0 && rModelPart.GetCommunicator().MyPID() == 0)
            << "Global force error: " << global_error << " (tolerated: < " << mForceRelativeTolerance
            << "), CSP: " << current_stiffness_parameter << ", |r|: " << residual_norm
            << ", |f_ext^inact|: " << inactive_external_force_norm
            << ", |f_int^act|: " << active_internal_force_norm << ". Global force error criterion "
            << (is_converged ? "satisfied" : "not satisfied") << std::endl;

        return is_converged;
    }

    [[nodiscard]] static double CalculateGlobalForceError(double ResidualNorm,
                                                          double CurrentStiffnessParameter,
                                                          double InactiveExternalForceNorm,
                                                          double ActiveInternalForceNorm)
    {
        const auto reference_force_norm =
            CurrentStiffnessParameter * InactiveExternalForceNorm + ActiveInternalForceNorm;
        if (reference_force_norm > 0.0) return ResidualNorm / reference_force_norm;

        return ResidualNorm > 0.0 ? std::numeric_limits<double>::infinity() : 0.0;
    }

    [[nodiscard]] double GetForceRelativeTolerance() const { return mForceRelativeTolerance; }

    [[nodiscard]] double GetForceAbsoluteTolerance() const { return mForceAbsoluteTolerance; }

    [[nodiscard]] Parameters GetDefaultParameters() const override
    {
        auto default_parameters = Parameters(R"(
        {
            "name"                     : "geo_global_force_error_criteria",
            "force_relative_tolerance" : 0.01,
            "force_absolute_tolerance" : 1.0e-9
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());
        return default_parameters;
    }

    static std::string Name() { return "geo_global_force_error_criteria"; }

    [[nodiscard]] std::string Info() const override { return "GeoGlobalForceErrorCriteria"; }

protected:
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);

        mForceRelativeTolerance = ThisParameters["force_relative_tolerance"].GetDouble();
        mForceAbsoluteTolerance = ThisParameters["force_absolute_tolerance"].GetDouble();
        KRATOS_ERROR_IF(mForceRelativeTolerance <= 0.0)
            << "The force_relative_tolerance must be positive, got " << mForceRelativeTolerance << std::endl;
        KRATOS_ERROR_IF(mForceAbsoluteTolerance < 0.0)
            << "The force_absolute_tolerance must not be negative, got " << mForceAbsoluteTolerance << std::endl;
    }

private:
    double            mForceRelativeTolerance = 0.01;
    double            mForceAbsoluteTolerance = 1.0e-9;
    TSystemVectorType mInactiveExternalForces;
    TSystemVectorType mActiveInternalForces;
    bool              mHasReferenceForces     = false;

    static std::vector<std::size_t> GetEquationIdsOfFreeDisplacementDofs(const DofsArrayType& rDofSet,
                                                                         std::size_t SystemSize)
    {
        auto result = std::vector<std::size_t>{};
        for (const auto& r_dof : rDofSet) {
            const auto key = r_dof.GetVariable().Key();
            const auto is_displacement_dof =
                key == DISPLACEMENT_X.Key() || key == DISPLACEMENT_Y.Key() || key == DISPLACEMENT_Z.Key();
            // Note that the elimination builder numbers the fixed degrees of freedom beyond the system size
            if (is_displacement_dof && r_dof.IsFree() && r_dof.EquationId() < SystemSize) {
                result.push_back(r_dof.EquationId());
            }
        }
        return result;
    }

    // Returns the quadratic norm of the given entries. Entries beyond the size of the vector (e.g. of
    // an empty vector) are regarded as zeros.
    static double CalculateNorm(const TSystemVectorType& rVector, const std::vector<std::size_t>& rEquationIds)
    {
        const auto size           = TSparseSpace::Size(rVector);
        auto       sum_of_squares = 0.0;
        for (const auto equation_id : rEquationIds) {
            if (equation_id < size) sum_of_squares += rVector[equation_id] * rVector[equation_id];
        }
        return std::sqrt(sum_of_squares);
    }
}; // Class GeoGlobalForceErrorCriteria

} // namespace Kratos

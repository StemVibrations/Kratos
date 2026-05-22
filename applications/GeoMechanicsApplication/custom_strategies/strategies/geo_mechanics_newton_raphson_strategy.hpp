// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "custom_strategies/strategies/residualbased_newton_raphson_strategy_two.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class GeoMechanicsNewtonRaphsonStrategy
    : public ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(GeoMechanicsNewtonRaphsonStrategy);

    using BaseType   = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
    using MotherType = ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver>;
    using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
    using TBuilderAndSolverType    = typename BaseType::TBuilderAndSolverType;
    using TSchemeType              = typename BaseType::TSchemeType;
    using DofsArrayType            = typename BaseType::DofsArrayType;
    using TSystemMatrixType        = typename BaseType::TSystemMatrixType;
    using TSystemVectorType        = typename BaseType::TSystemVectorType;
    using MotherType::mInitializeWasPerformed;
    using MotherType::mMaxIterationNumber;
    using MotherType::mpA; // Tangent matrix
    using MotherType::mpb; // Residual vector of iteration i
    using MotherType::mpBuilderAndSolver;
    using MotherType::mpDx; // Delta x of iteration i
    using MotherType::mpScheme;

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    GeoMechanicsNewtonRaphsonStrategy(ModelPart&                    model_part,
                                      typename TSchemeType::Pointer pScheme,
                                      typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
                                      typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
                                      Parameters&                             rParameters,
                                      int                                     MaxIterations = 30,
                                      bool CalculateReactions                               = false,
                                      bool ReformDofSetAtEachStep                           = false,
                                      bool MoveMeshFlag                                     = false,
                                      bool ReadForce = false)
        : ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver>(
              model_part, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag, ReadForce)
    {
        // only include validation with c++11 since raw_literals do not exist in c++03 
        Parameters default_parameters(R"(
                {
                    "min_iteration":    2,
                    "number_cycles":    5,
                    "increase_factor":  2.0,
                    "reduction_factor": 0.5,
                    "max_piping_iterations": 50,
                    "desired_iterations": 4,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "rebuild_level": 2,
                    "iteration_method": "broyden",
                    "arc_length_settings": {
                        "initial_arc_length"        : 0.0,
                        "max_arc_length_factor"     : 10.0,
                        "min_arc_length_factor"     : 0.1,
                        "desired_iterations"        : 10,
                        "beta"                      : 0.0,
                        "max_load_factor"           : 1.0,
                        "max_arc_length_reductions" : 5
                    },
                    "broyden_settings": {
                       "absolute_global_error": 1e-12,
                       "relative_global_error": 1e-2,
                       "relative_local_error": 1e-1
                    }   
                }  )");

        // Validate against defaults -- this also ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        // Forward arc-length settings (if requested) to the base strategy.
        if (rParameters["iteration_method"].GetString() == "arc_length") {
            this->SetUseArcLength(true);
            const Parameters al = rParameters["arc_length_settings"];
            this->SetArcLengthSettings(
                al["initial_arc_length"].GetDouble(),
                al["max_arc_length_factor"].GetDouble(),
                al["min_arc_length_factor"].GetDouble(),
                static_cast<unsigned int>(al["desired_iterations"].GetInt()),
                al["beta"].GetDouble(),
                al["max_load_factor"].GetDouble(),
                static_cast<unsigned int>(al["max_arc_length_reductions"].GetInt())
            );
        }
		else if (rParameters["iteration_method"].GetString() == "broyden") {
			this->SetUseBroyden(true);
			const Parameters broyden_settings = rParameters["broyden_settings"];
			this->SetBroydenSettings(
				broyden_settings["absolute_global_error"].GetDouble(),
				broyden_settings["relative_global_error"].GetDouble(),
				broyden_settings["relative_local_error"].GetDouble()
			);

		}
    }
}; // Class GeoMechanicsNewtonRaphsonStrategy

} // namespace Kratos
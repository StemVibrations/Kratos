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

// System includes
#include <algorithm>

// Project includes
#include "factories/linear_solver_factory.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

    template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
    class ResidualBasedNewtonRaphsonStrategyTwo
        : public ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
    {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedNewtonRaphsonStrategyTwo);

        /// The definition of the linear solver factory type
        using LinearSolverFactoryType = LinearSolverFactory<TSparseSpace, TDenseSpace>;
        using BaseType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
        using MotherType = ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;
        using TConvergenceCriteriaType = ConvergenceCriteria<TSparseSpace, TDenseSpace>;
        using TBuilderAndSolverType = typename BaseType::TBuilderAndSolverType;
        using TSchemeType = typename BaseType::TSchemeType;
        using DofsArrayType = typename BaseType::DofsArrayType;
        using TSystemMatrixType = typename BaseType::TSystemMatrixType;
        using TSystemVectorType = typename BaseType::TSystemVectorType;
        using MotherType::mCalculateReactionsFlag;
        using MotherType::mInitializeWasPerformed;
        using MotherType::mMaxIterationNumber;
        using MotherType::mpA; // Tangent matrix
        using MotherType::mpb; // Residual vector of iteration i
        using MotherType::mpBuilderAndSolver;
        using MotherType::mpConvergenceCriteria;
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
        ResidualBasedNewtonRaphsonStrategyTwo(ModelPart& rModelPart,
            typename TSchemeType::Pointer pScheme,
            typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
            typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
            Parameters& rParameters,
            int                                     MaxIterations = 30,
            bool CalculateReactions = false,
            bool ReformDofSetAtEachStep = false,
            bool MoveMeshFlag = false)
            : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(
                rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            // only include validation with c++11 since raw_literals do not exist in c++03
            Parameters default_parameters(R"(
                {
                    "quasi_newton_type": "broyden",
                    "quasi_newton_restart_interval": 50,
                    "quasi_newton_max_rank" : 10,
                        "extrapolate_previous_increment": false,
                        "relaxation_factor": 0.8333333333333334

                }  )");

            // Validate against defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);

            // Select the quasi-Newton update scheme (low-rank secant approximation of the Jacobian/its inverse).
            if (const std::string quasi_newton_type = rParameters["quasi_newton_type"].GetString();
                quasi_newton_type == "lbfgs") {
                mQuasiNewtonType = QuasiNewtonType::LBFGS;
            }
            else if (quasi_newton_type == "broyden") {
                mQuasiNewtonType = QuasiNewtonType::Broyden;
            }
			else if (quasi_newton_type == "none") {
				mQuasiNewtonType = QuasiNewtonType::None; // default to Broyden if none is specified
			}
            else {
                KRATOS_ERROR << "Unknown 'quasi_newton_type': '" << quasi_newton_type
                    << "'. Valid options are 'broyden' and 'lbfgs'." << std::endl;
            }

            mRestartInterval = rParameters["quasi_newton_restart_interval"].GetInt();
            mMaxRank = rParameters["quasi_newton_max_rank"].GetInt();
            mRelaxationFactor = rParameters["relaxation_factor"].GetDouble();
            KRATOS_ERROR_IF(mRelaxationFactor <= 0.0)
                << "'relaxation_factor' must be positive, got " << mRelaxationFactor << std::endl;

            // it is required to use a direct solver without scaling for the quasi-newton strategy, as the strategy depends on reusing the factorization of the initial stiffness matrix.
            auto linear_solver_settings = Parameters(R"(
            {
                "solver_type": "LinearSolversApplication.sparse_lu",
                "scaling": false
            }  )");

            mpLinearSolver = LinearSolverFactoryType().Create(linear_solver_settings);
        }

        bool SolveSolutionStep() override
        {
            KRATOS_TRY

                ModelPart& r_model_part = BaseType::GetModelPart();
            typename TSchemeType::Pointer p_scheme = MotherType::GetScheme();
            typename TBuilderAndSolverType::Pointer p_builder_and_solver = MotherType::GetBuilderAndSolver();
            auto& r_dof_set = p_builder_and_solver->GetDofSet();

            TSystemMatrixType& rA = *mpA;
            TSystemVectorType& rDx = *mpDx;
            TSystemVectorType& rb = *mpb;

            // Are there master-slave constraints active? If so, we must operate in the
            // reduced (constrained) space: A_hat = T^T A T, b_hat = T^T b, and recover dx = T dx_hat.
            const auto has_constraints = !r_model_part.MasterSlaveConstraints().empty();

            unsigned int iteration_number = 1;
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            // must be evaluated before the first residual of the step is built
            this->InitializeStageUnbalanceIfNeeded(rb.size());

            // L-BFGS pairs describe the curvature along the iterates of a single step; pairs of a previous
            // (converged or cut back) step are not valid for the current one
            mLbfgsRankStorage.Clear();

            bool rebuild_lhs = false;
            if (BaseType::mRebuildLevel > 0 || !BaseType::mStiffnessMatrixIsBuilt) {
                rebuild_lhs = true;
            }
            else {
                BuildReducedResidual(rb);
            }

            TSystemVectorType rb_new(rb.size());

            for (; iteration_number <= mMaxIterationNumber; ++iteration_number) {
                r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
                if (iteration_number > 1) {
                    p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
                    mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
                }

                // Quasi-Newton solve: rDx = H_k * rb   (H_k approximates rA_k^{-1})
                this->QuasiNewtonSolve(rA, rDx, rb, iteration_number, rebuild_lhs);
                rebuild_lhs = false;

                UpdateDatabaseReduced(rA, rDx, rb, has_constraints);
                BuildReducedResidual(rb_new);

                TSystemVectorType delta_b = rb - rb_new;
                this->UpdateRankStorage(rA, rDx, delta_b);

                TSparseSpace::Copy(rb_new, rb);
                MotherType::EchoInfo(iteration_number);
                p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
                mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

                is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
                if (is_converged) {
                    break;
                }
            }

            if (!is_converged) {
                MotherType::MaxIterationsExceeded();
            }
            else {
                KRATOS_INFO_IF("GeoMechanicsQuasiNewtonStrategy", this->GetEchoLevel() > 0)
                    << "Convergence achieved after " << iteration_number << " / " << mMaxIterationNumber
                    << std::endl;
            }

            if (mCalculateReactionsFlag) {
                p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
            }

            return is_converged;

            KRATOS_CATCH("")
        }

    private:
        /// The available quasi-Newton update schemes.
        enum class QuasiNewtonType { Broyden, LBFGS, None };

        struct RankStorage {
            virtual void Clear() = 0;
            virtual ~RankStorage() = default;
        };

        struct BroydenRankStorage : RankStorage {
            std::vector<TSystemVectorType> u_list; // u_i
            std::vector<TSystemVectorType> v_list; // v_i
            std::vector<TSystemVectorType> z_list; // z_i = rA0^{-1} u_i (cached)

            using RankStorage::Clear;

            void Clear() override
            {
                u_list.clear();
                v_list.clear();
                z_list.clear();
            }
        };

        struct LBFGSRankStorage : RankStorage {
            std::vector<TSystemVectorType> dx_list;  // step differences
            std::vector<TSystemVectorType> db_list;  // gradient/residual differences
            std::vector<double>            rho_list; // curvature scalars

            using RankStorage::Clear;

            void Clear() override
            {
                dx_list.clear();
                db_list.clear();
                rho_list.clear();
            }
        };

        QuasiNewtonType    mQuasiNewtonType = QuasiNewtonType::Broyden;
        LBFGSRankStorage   mLbfgsRankStorage;
        BroydenRankStorage mBroydenRankStorage;
        unsigned int       mRestartInterval = 100; // rebuild K0 every N iterations to refresh curvature
        unsigned int mMaxRank = 10; // maximum number of low-rank updates to store (Broyden or BFGS)
        double       mRelaxationFactor = 1.0; // scales the quasi-Newton increment before it is applied
        typename TLinearSolver::Pointer mpLinearSolver;

        // The out-of-balance force vector at the start of the stage (in the constrained space), i.e. the
        // unbalance that results from (de)activating elements and changing materials or loads. Only the load fraction
        // of this unbalance is applied in a step, such that each
        // step converges to an equilibrium state:
        //     r = f_ext - f_int - (1 - load_fraction) * r_stage
        TSystemVectorType mStageUnbalance;
        bool              mIsStageUnbalanceInitialized = false;
        double            mStageStartTime              = 0.0;

        double GetCurrentLoadFraction(const ProcessInfo& rProcessInfo) const
        {
            const double t0 = rProcessInfo[START_TIME];
            const double t1 = rProcessInfo[END_TIME];
            const double t = rProcessInfo[TIME];
            const double dt = t1 - t0;
            return (std::abs(dt) > 0.0) ? std::clamp((t - t0) / dt, 0.0, 1.0) : 1.0;
        }

        void InitializeStageUnbalanceIfNeeded(std::size_t SystemSize)
        {
            const double stage_start_time = BaseType::GetModelPart().GetProcessInfo()[START_TIME];
            if (mIsStageUnbalanceInitialized && mStageStartTime == stage_start_time &&
                mStageUnbalance.size() == SystemSize) {
                return;
            }

            // Evaluated for the (converged) state at the start of the stage, before any iteration. After a
            // cut back of the first step, the state is reset, so the stored unbalance remains valid.
            mStageUnbalance.resize(SystemSize, false);
            BuildUnbalancedReducedResidual(mStageUnbalance);
            mStageStartTime              = stage_start_time;
            mIsStageUnbalanceInitialized = true;
        }

        /// Removes the part of the stage unbalance that is not yet applied in the current step
        void ApplyLoadFraction(TSystemVectorType& rb) const
        {
            const double load_fraction =
                this->GetCurrentLoadFraction(BaseType::GetModelPart().GetProcessInfo());

            std::cout << "load fraction: " << load_fraction << std::endl;

            if (mIsStageUnbalanceInitialized && mStageUnbalance.size() == rb.size()) {
                TSparseSpace::UnaliasedAdd(rb, load_fraction - 1.0, mStageUnbalance);
            }
        }

        void UpdateBroydenRank(TSystemMatrixType& rA_0, const TSystemVectorType& rDx, const TSystemVectorType& rDb)
        {
            // ---- Build the new Broyden rank-1 column  ----
            //   rA_k rDx = rA_0 rDx + U (V^T rDx)
            //   db_hat = rA_k rDx
            //   u_col = rDb - y_hat
            //   v_col = rDx / (rDx . rDx)
            //   z_col = rA_0^{-1} u_col

            TSystemVectorType db_hat(rDx.size());

            const double d_x_squared = TSparseSpace::Dot(rDx, rDx);
            if (d_x_squared > std::numeric_limits<double>::epsilon()) {
                const std::size_t k = mBroydenRankStorage.u_list.size();
                TSparseSpace::Mult(rA_0, rDx, db_hat);
                if (k > 0) {
                    Vector v_list_dot_dx(k);
                    for (std::size_t i = 0; i < k; ++i) {
                        v_list_dot_dx[i] = TSparseSpace::Dot(mBroydenRankStorage.v_list[i], rDx);
                    }
                    for (std::size_t i = 0; i < k; ++i) {
                        noalias(db_hat) += v_list_dot_dx[i] * mBroydenRankStorage.u_list[i];
                    }
                }
                TSystemVectorType u_col = rDb - db_hat;
                TSystemVectorType v_col = rDx / d_x_squared;

                TSystemVectorType z_col(rDx.size());
                mpLinearSolver->PerformSolutionStep(rA_0, z_col, u_col);

                if (mBroydenRankStorage.u_list.size() >= mMaxRank) {
                    mBroydenRankStorage.u_list.erase(mBroydenRankStorage.u_list.begin());
                    mBroydenRankStorage.v_list.erase(mBroydenRankStorage.v_list.begin());
                    mBroydenRankStorage.z_list.erase(mBroydenRankStorage.z_list.begin());
                }
                mBroydenRankStorage.u_list.push_back(u_col);
                mBroydenRankStorage.v_list.push_back(v_col);
                mBroydenRankStorage.z_list.push_back(z_col);
            }
        }

        void UpdateLBFGSRank(const TSystemVectorType& rDx, const TSystemVectorType& rDb)
        {
            const double db_dot_dx = TSparseSpace::Dot(rDb, rDx);
            if (db_dot_dx > std::numeric_limits<double>::epsilon()) {
                if (mLbfgsRankStorage.dx_list.size() >= mMaxRank) {
                    mLbfgsRankStorage.dx_list.erase(mLbfgsRankStorage.dx_list.begin());
                    mLbfgsRankStorage.db_list.erase(mLbfgsRankStorage.db_list.begin());
                    mLbfgsRankStorage.rho_list.erase(mLbfgsRankStorage.rho_list.begin());
                }
                mLbfgsRankStorage.dx_list.push_back(rDx);
                mLbfgsRankStorage.db_list.push_back(rDb);
                mLbfgsRankStorage.rho_list.push_back(1.0 / db_dot_dx);
            }
        }

        /// Builds A and b, applies master-slave constraints
        /// Dirichlet conditions, then factorizes the resulting (reduced) rA_0.
        void BuildAndConstrainA0(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rb, RankStorage& rRankStorage)
        {
            rRankStorage.Clear();

            auto       p_builder_and_solver = MotherType::GetBuilderAndSolver();
            auto       p_scheme = MotherType::GetScheme();
            ModelPart& r_model_part = BaseType::GetModelPart();

            TSparseSpace::SetToZero(rA_0);
            TSparseSpace::SetToZero(rb);

            // Build raw A and b together (ApplyConstraints needs the raw b to form T^T b consistently).
            p_builder_and_solver->Build(p_scheme, r_model_part, rA_0, rb);
            if (!r_model_part.MasterSlaveConstraints().empty()) {
                p_builder_and_solver->ApplyConstraints(p_scheme, r_model_part, rA_0, rb); // A <- T^T A T ; b <- T^T b ; builds mT
            }
            p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA_0, rDx, rb);
            this->ApplyLoadFraction(rb);

            mpLinearSolver->InitializeSolutionStep(rA_0, rDx, rb); // factorize rA_0

            BaseType::mStiffnessMatrixIsBuilt = true;
        }

        /// Builds the residual and reduces it to the constrained space: b_hat = T^T b (slaves zeroed).
        /// Only the load fraction of the stage unbalance is included.
        void BuildReducedResidual(TSystemVectorType& rb)
        {
            BuildUnbalancedReducedResidual(rb);
            this->ApplyLoadFraction(rb);
        }

        /// Builds the full residual f_ext - f_int and reduces it to the constrained space
        void BuildUnbalancedReducedResidual(TSystemVectorType& rb)
        {
            auto       p_builder_and_solver = MotherType::GetBuilderAndSolver();
            auto       p_scheme = MotherType::GetScheme();
            ModelPart& r_model_part = BaseType::GetModelPart();

            TSparseSpace::SetToZero(rb);

            // note that BuildRHS also applies the Dirichlet conditions on the RHS
            p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);

            if (!r_model_part.MasterSlaveConstraints().empty()) {
                p_builder_and_solver->ApplyRHSConstraints(p_scheme, r_model_part, rb);
            }
        }

        void UpdateDatabaseReduced(TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb, bool HasConstraints)
        {
            auto p_builder_and_solver = MotherType::GetBuilderAndSolver();

            // Relax the increment in place, such that rDx holds the increment that is actually applied.
            // The secant pair (rDx, delta_b) and the convergence criteria then use the applied increment.
            if (mRelaxationFactor != 1.0) {
                TSparseSpace::InplaceMult(rDx, mRelaxationFactor);
            }

            if (HasConstraints) {
                auto& rT = p_builder_and_solver->GetConstraintRelationMatrix();
                TSystemVectorType dx_full(rDx.size());
                TSparseSpace::Mult(rT, rDx, dx_full); // dx_full = T * dx_reduced
                MotherType::UpdateDatabase(rA, dx_full, rb, BaseType::MoveMeshFlag());
            }
            else {
                MotherType::UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());
            }
        }

        void LBfgsSolve(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rb)
        {
            const std::size_t k = mLbfgsRankStorage.dx_list.size();
            Vector            alpha(k);

            TSystemVectorType q(rb.size());
            TSparseSpace::Copy(rb, q); // q = rb

            // First loop: newest to oldest
            for (auto i = k; i-- > 0;) {
                alpha[i] = mLbfgsRankStorage.rho_list[i] * TSparseSpace::Dot(mLbfgsRankStorage.dx_list[i], q);
                noalias(q) -= alpha[i] * mLbfgsRankStorage.db_list[i];
            }

            // Seed: rDx = H_0 * q = rA_0^{-1} q (reuses the cached factorization)
            TSparseSpace::SetToZero(rDx);
            mpLinearSolver->PerformSolutionStep(rA_0, rDx, q);

            // Second loop: oldest to newest
            for (std::size_t i = 0; i < k; ++i) {
                const double beta =
                    mLbfgsRankStorage.rho_list[i] * TSparseSpace::Dot(mLbfgsRankStorage.db_list[i], rDx);
                noalias(rDx) += (alpha[i] - beta) * mLbfgsRankStorage.dx_list[i];
            }
        }

        /// Solves rDx = (rA_0 + U V^T)^{-1} rb (Sherman-Morrison-Woodbury).
        void ShermanMorrisSolve(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rb)
        {
            mpLinearSolver->PerformSolutionStep(rA_0, rDx, rb); // rDx = rA_0^{-1} rb

            // add the low rank correction, rDx = rDx - Z * (I + V^T Z)^{-1} (V^T rDx), where Z = rA_0^{-1} U
            const std::size_t k = mBroydenRankStorage.u_list.size();
            if (k > 0) {
                Vector v_list_dot_dx(k);
                for (std::size_t i = 0; i < k; ++i) {
                    v_list_dot_dx[i] = TSparseSpace::Dot(mBroydenRankStorage.v_list[i], rDx);
                }
                // Calculates: aux_matrix = I + V^T Z, where Z = rA_0^{-1} U
                Matrix aux_matrix(k, k);
                for (std::size_t i = 0; i < k; ++i) {
                    for (std::size_t j = 0; j < k; ++j) {
                        aux_matrix(i, j) =
                            (i == j ? 1.0 : 0.0) + TSparseSpace::Dot(mBroydenRankStorage.v_list[i],
                                mBroydenRankStorage.z_list[j]);
                    }
                }

                // Calculates: alpha = aux_matrix^{-1} (V^T rDx)
                Vector alpha;
                SolveSmallDense(aux_matrix, alpha, v_list_dot_dx);
                for (std::size_t j = 0; j < k; ++j) {
                    noalias(rDx) -= alpha[j] * mBroydenRankStorage.z_list[j];
                }
            }
        }

        /**
         * @brief Solve a small dense k x k linear system M * x = b using
         *        Gaussian elimination with partial pivoting.
         *        Note: rM and rRhs are taken by value (mutated locally).
         */
        void SolveSmallDense(Matrix rM, Vector& rSolution, Vector rRhs)
        {
            const std::size_t k = rRhs.size();
            rSolution.resize(k, false);
            if (k == 0) return;

            for (std::size_t i = 0; i < k; ++i) {
                std::size_t pivot = i;
                double      pivot_val = std::abs(rM(i, i));
                for (std::size_t r = i + 1; r < k; ++r) {
                    const double v = std::abs(rM(r, i));
                    if (v > pivot_val) {
                        pivot_val = v;
                        pivot = r;
                    }
                }
                if (pivot_val < std::numeric_limits<double>::epsilon()) {
                    KRATOS_ERROR
                        << "Singular matrix in GeoMechanicsQuasiNewtonStrategy SolveSmallDense." << std::endl;
                }
                if (pivot != i) {
                    for (std::size_t c = i; c < k; ++c)
                        std::swap(rM(i, c), rM(pivot, c));
                    std::swap(rRhs[i], rRhs[pivot]);
                }
                const double diag = rM(i, i);
                for (std::size_t r = i + 1; r < k; ++r) {
                    const double factor = rM(r, i) / diag;
                    if (factor == 0.0) continue;
                    for (std::size_t c = i; c < k; ++c) {
                        rM(r, c) -= factor * rM(i, c);
                    }
                    rRhs[r] -= factor * rRhs[i];
                }
            }
            for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(k) - 1; i >= 0; --i) {
                double s = rRhs[i];
                for (std::size_t c = i + 1; c < k; ++c) {
                    s -= rM(i, c) * rSolution[c];
                }
                rSolution[i] = s / rM(i, i);
            }
        }

        void QuasiNewtonSolve(TSystemMatrixType& rA_0,
            TSystemVectorType& rDx,
            TSystemVectorType& rb,
            const unsigned int IterationNumber,
            bool               RebuildLhs)
        {
            // rebuild the LHS if the restart interval is reached (to refresh curvature)
            if (IterationNumber > 1 && (IterationNumber - 1) % mRestartInterval == 0) {
                RebuildLhs = true;
            }

            if (mQuasiNewtonType == QuasiNewtonType::LBFGS) {
                if (RebuildLhs) {
                    BuildAndConstrainA0(rA_0, rDx, rb, mLbfgsRankStorage);
                }
                LBfgsSolve(rA_0, rDx, rb);
            }
            else if (mQuasiNewtonType == QuasiNewtonType::Broyden) {
                if (RebuildLhs) {
                    BuildAndConstrainA0(rA_0, rDx, rb, mBroydenRankStorage);
                }
                ShermanMorrisSolve(rA_0, rDx, rb);
            }
			else if (mQuasiNewtonType == QuasiNewtonType::None) {
				if (RebuildLhs) {
					BuildAndConstrainA0(rA_0, rDx, rb, mBroydenRankStorage); // use Broyden storage to clear rank
				}
				mpLinearSolver->PerformSolutionStep(rA_0, rDx, rb);
			}
            else {
                KRATOS_ERROR << "Unknown 'quasi_newton_type'. Valid options are 'broyden' and 'lbfgs'."
                    << std::endl;
            }
        }

        void UpdateRankStorage(TSystemMatrixType& rA_0, TSystemVectorType& rDx, TSystemVectorType& rDeltaB)
        {
            if (mQuasiNewtonType == QuasiNewtonType::LBFGS) {
                this->UpdateLBFGSRank(rDx, rDeltaB);

            }
            else if (mQuasiNewtonType == QuasiNewtonType::Broyden) {
                this->UpdateBroydenRank(rA_0, rDx, rDeltaB);
            }
			else if (mQuasiNewtonType == QuasiNewtonType::None) {
				// do nothing
			}
            else {
                KRATOS_ERROR << "Unknown 'quasi_newton_type'. Valid options are 'broyden' and 'lbfgs'."
                    << std::endl;
            }
        }

    }; // Class ResidualBasedNewtonRaphsonStrategyTwo

} // namespace Kratos

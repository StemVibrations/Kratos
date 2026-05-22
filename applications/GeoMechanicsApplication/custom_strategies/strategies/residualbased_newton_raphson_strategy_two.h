//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#pragma once

// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "utilities/builtin_timer.h"
#include "includes/serializer.h"
#include "includes/file_serializer.h"
#include "geo_mechanics_application_variables.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}

///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ResidualBasedNewtonRaphsonStrategy
 * @ingroup KratosCore
 * @brief This is the base Newton Raphson strategy
 * @details This strategy iterates until the convergence is achieved (or the maximum number of iterations is surpassed) using a Newton Raphson algorithm
 * @author Riccardo Rossi
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class ResidualBasedNewtonRaphsonStrategyTwo
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>
{
  public:
    ///@name Type Definitions
    ///@{
    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace> TConvergenceCriteriaType;

    // Counted pointer of ClassName
    KRATOS_CLASS_POINTER_DEFINITION(ResidualBasedNewtonRaphsonStrategyTwo);

    typedef SolvingStrategy<TSparseSpace, TDenseSpace> SolvingStrategyType;

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver> ClassType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;

    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /// The definition of the element container type
    typedef PointerVectorSet<Element, IndexedObject> ElementsContainerType;

    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ResidualBasedNewtonRaphsonStrategyTwo() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ResidualBasedNewtonRaphsonStrategyTwo(ModelPart& rModelPart)
        : ResidualBasedNewtonRaphsonStrategyTwo(rModelPart, ResidualBasedNewtonRaphsonStrategyTwo::GetDefaultParameters())
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit ResidualBasedNewtonRaphsonStrategyTwo(ModelPart& rModelPart, Parameters ThisParameters)
        : BaseType(rModelPart),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            // Tells to the builder and solver if the reactions have to be Calculated or not
            p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

            // Tells to the Builder And Solver if the system matrix and vectors need to
            // be reshaped at each step or not
            p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
        } else {
            KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategy") << "BuilderAndSolver is not initialized. Please assign one before settings flags" << std::endl;
        }

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();
    }

    /**
     * Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        bool ReadForce=false)
        : BaseType(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mMaxIterationNumber(MaxIterations),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false),
		  mReadForce(ReadForce)
    {
        KRATOS_TRY; 

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));

        // Tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    explicit ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        bool ReadForce = false)
        : BaseType(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mpBuilderAndSolver(pNewBuilderAndSolver),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions),
          mMaxIterationNumber(MaxIterations),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false),
		  mReadForce(ReadForce)
    {
        KRATOS_TRY

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // Tells to the builder and solver if the reactions have to be Calculated or not
        p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /**
     * @brief Constructor specifying the builder and solver
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param MaxIterations The maximum number of non-linear iterations to be considered when solving the problem
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
    explicit ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        int MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false)
        : ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
    {
        KRATOS_TRY

        KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategyTwo") << "This constructor is deprecated, please use the constructor without linear solver" << std::endl;

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // We check if the linear solver considered for the builder and solver is consistent
        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
        KRATOS_ERROR_IF(p_linear_solver != pNewLinearSolver) << "Inconsistent linear solver in strategy and builder and solver. Considering the linear solver assigned to builder and solver :\n" << p_linear_solver->Info() << "\n instead of:\n" << pNewLinearSolver->Info() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * Constructor with Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param Settings Settings used in the strategy
     */
    ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        Parameters Settings)
        : BaseType(rModelPart),
          mpScheme(pScheme),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY;

        // Setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer(
            new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>(pNewLinearSolver));

        // Tells to the builder and solver if the reactions have to be Calculated or not
        mpBuilderAndSolver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        // be reshaped at each step or not
        mpBuilderAndSolver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("");
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Settings Settings used in the strategy
     */
    ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : BaseType(rModelPart),
          mpScheme(pScheme),
          mpBuilderAndSolver(pNewBuilderAndSolver),
          mpConvergenceCriteria(pNewConvergenceCriteria),
          mInitializeWasPerformed(false),
          mKeepSystemConstantDuringIterations(false)
    {
        KRATOS_TRY
        // Validate and assign defaults
        Settings = this->ValidateAndAssignParameters(Settings, this->GetDefaultParameters());
        this->AssignSettings(Settings);
        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // Tells to the builder and solver if the reactions have to be Calculated or not
        p_builder_and_solver->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        // Tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        p_builder_and_solver->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        // Set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        // By default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        mpA = TSparseSpace::CreateEmptyMatrixPointer();
        mpDx = TSparseSpace::CreateEmptyVectorPointer();
        mpb = TSparseSpace::CreateEmptyVectorPointer();

        KRATOS_CATCH("")
    }

    /**
     * @brief Constructor specifying the builder and solver and using Parameters
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param pNewBuilderAndSolver The builder and solver employed
     * @param Parameters Settings used in the strategy
     */
    KRATOS_DEPRECATED_MESSAGE("Constructor deprecated, please use the constructor without linear solver")
        ResidualBasedNewtonRaphsonStrategyTwo(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        Parameters Settings)
        : ResidualBasedNewtonRaphsonStrategyTwo<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewConvergenceCriteria, pNewBuilderAndSolver, Settings)
    {
        KRATOS_TRY

        KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategyTwo") << "This constructor is deprecated, please use the constructor without linear solver" << std::endl;

        // Getting builder and solver
        auto p_builder_and_solver = GetBuilderAndSolver();

        // We check if the linear solver considered for the builder and solver is consistent
        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
        KRATOS_ERROR_IF(p_linear_solver != pNewLinearSolver) << "Inconsistent linear solver in strategy and builder and solver. Considering the linear solver assigned to builder and solver :\n" << p_linear_solver->Info() << "\n instead of:\n" << pNewLinearSolver->Info() << std::endl;

        KRATOS_CATCH("")
    }

    /**
     * @brief Destructor.
     * @details In trilinos third party library, the linear solver's preconditioner should be freed before the system matrix. We control the deallocation order with Clear().
     */
    ~ResidualBasedNewtonRaphsonStrategyTwo() override
    {
        // If the linear solver has not been deallocated, clean it before
        // deallocating mpA. This prevents a memory error with the the ML
        // solver (which holds a reference to it).
        // NOTE: The linear solver is hold by the B&S
        auto p_builder_and_solver = this->GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->Clear();
        }

        // Deallocating system vectors to avoid errors in MPI. Clear calls
        // TrilinosSpace::Clear for the vectors, which preserves the Map of
        // current vectors, performing MPI calls in the process. Due to the
        // way Python garbage collection works, this may happen after
        // MPI_Finalize has already been called and is an error. Resetting
        // the pointers here prevents Clear from operating with the
        // (now deallocated) vectors.
        mpA.reset();
        mpDx.reset();
        mpb.reset();

        Clear();
    }

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    /**
     * @brief Set method for the builder and solver
     * @param pNewBuilderAndSolver The pointer to the builder and solver considered
     */
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver)
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };

    /**
     * @brief Get method for the builder and solver
     * @return mpBuilderAndSolver: The pointer to the builder and solver considered
     */
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    /**
     * @brief This method sets the flag mFullUpdateFlag
     * @param UseOldStiffnessInFirstIterationFlag The flag that tells if
     */
    void SetUseOldStiffnessInFirstIterationFlag(bool UseOldStiffnessInFirstIterationFlag)
    {
        mUseOldStiffnessInFirstIteration = UseOldStiffnessInFirstIterationFlag;
    }

    /**
     * @brief This method returns the flag mFullUpdateFlag
     * @return The flag that tells if
     */
    bool GetUseOldStiffnessInFirstIterationFlag()
    {
        return mUseOldStiffnessInFirstIteration;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    /**
     * @brief This method sets the flag mMaxIterationNumber
     * @param MaxIterationNumber This is the maximum number of on linear iterations
     */
    void SetMaxIterationNumber(unsigned int MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }

    /**
     * @brief This method gets the flag mMaxIterationNumber
     * @return mMaxIterationNumber: This is the maximum number of on linear iterations
     */
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    /**
     * @brief It sets the level of echo for the solving strategy
     * @param Level The level to set
     * @details The different levels of echo are:
     * - 0: Mute... no echo at all
     * - 1: Printing time and basic information
     * - 2: Printing linear solver data
     * - 3: Print of debug information: Echo of stiffness matrix, Dx, b...
     */

    void SetEchoLevel(int Level) override
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
        mpConvergenceCriteria->SetEchoLevel(Level);
    }

    //*********************************************************************************
    /**OPERATIONS ACCESSIBLE FROM THE INPUT: **/

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    typename SolvingStrategyType::Pointer Create(
        ModelPart& rModelPart,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(rModelPart, ThisParameters);
    }

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY
        const DataCommunicator &r_comm = BaseType::GetModelPart().GetCommunicator().GetDataCommunicator();
        //OPERATIONS THAT SHOULD BE DONE ONCE - internal check to avoid repetitions
        //if the operations needed were already performed this does nothing
        if (mInitializeWasPerformed == false)
            Initialize();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        DofsArrayType& r_dof_set = GetBuilderAndSolver()->GetDofSet();

        GetScheme()->Predict(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);

        // Applying constraints if needed
        auto& r_constraints_array = BaseType::GetModelPart().MasterSlaveConstraints();
        const int local_number_of_constraints = r_constraints_array.size();
        const int global_number_of_constraints = r_comm.SumAll(local_number_of_constraints);
        if(global_number_of_constraints != 0) {
            const auto& r_process_info = BaseType::GetModelPart().GetProcessInfo();

            block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.ResetSlaveDofs(r_process_info);
            });
            block_for_each(r_constraints_array, [&r_process_info](MasterSlaveConstraint& rConstraint){
                rConstraint.Apply(r_process_info);
            });

            // The following is needed since we need to eventually compute time derivatives after applying
            // Master slave relations
            TSparseSpace::SetToZero(rDx);
            this->GetScheme()->Update(BaseType::GetModelPart(), r_dof_set, rA, rDx, rb);
        }

        // Move the mesh if needed
        if (this->MoveMeshFlag() == true)
            BaseType::MoveMesh();

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        if (mInitializeWasPerformed == false)
        {
            //pointers needed in the solution
            typename TSchemeType::Pointer p_scheme = GetScheme();
            typename TConvergenceCriteriaType::Pointer p_convergence_criteria = mpConvergenceCriteria;

            //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (p_scheme->SchemeIsInitialized() == false)
                p_scheme->Initialize(BaseType::GetModelPart());

            //Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ElementsAreInitialized() == false)
                p_scheme->InitializeElements(BaseType::GetModelPart());

            //Initialize The Conditions - OPERATIONS TO BE DONE ONCE
            if (p_scheme->ConditionsAreInitialized() == false)
                p_scheme->InitializeConditions(BaseType::GetModelPart());

            //initialisation of the convergence criteria
            if (p_convergence_criteria->IsInitialized() == false)
                p_convergence_criteria->Initialize(BaseType::GetModelPart());

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY;

        // save external force vector to be able to restart the incremental solver
        const auto file_name = "incremental_solver_serialization";
        //Serializer::TraceType const& rTrace = Serializer::TraceType ::SERIALIZER_NO_TRACE;
        auto serializer = FileSerializer(file_name, Serializer::TraceType::SERIALIZER_NO_TRACE, false);
        save(serializer);

		std::cout << "Clearing strategy..." << std::endl;

        // Setting to zero the internal flag to ensure that the dof sets are recalculated. Also clear the linear solver stored in the B&S
        auto p_builder_and_solver = GetBuilderAndSolver();
        if (p_builder_and_solver != nullptr) {
            p_builder_and_solver->SetDofSetIsInitializedFlag(false);
            p_builder_and_solver->Clear();
        }
		std::cout << "Clearing builder and solver..." << std::endl;
        // Clearing the system of equations
        if (mpA != nullptr)
            SparseSpaceType::Clear(mpA);
        if (mpDx != nullptr)
            SparseSpaceType::Clear(mpDx);
        if (mpb != nullptr)
            SparseSpaceType::Clear(mpb);

		std::cout << "Clearing system of equations..." << std::endl;
        // Clearing scheme
        auto p_scheme = GetScheme();
        if (p_scheme != nullptr) {
            GetScheme()->Clear();
        }

		std::cout << "Clearing scheme..." << std::endl;
        mInitializeWasPerformed = false;

        KRATOS_CATCH("");
    }

    /**
     * @brief This should be considered as a "post solution" convergence check which is useful for coupled analysis - the convergence criteria used is the one used inside the "solve" step
     */
    bool IsConverged() override
    {
        KRATOS_TRY;

        TSystemMatrixType& rA = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb = *mpb;

        if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
        {
            ModelPart& r_model_part = BaseType::GetModelPart();
            TSparseSpace::SetToZero(rb);
            const double load_fraction = GetCurrentLoadFraction(r_model_part.GetProcessInfo());
            TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
            TSystemVectorType int_force = ZeroVector(rb.size());
            BuildInternalForceVector(r_model_part, int_force);
            rb = load_fraction * d_ext_force + mPreviousExternalForceVector - int_force;

            //GetBuilderAndSolver()->BuildRHS(GetScheme(), BaseType::GetModelPart(), rb);
        }

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(), GetBuilderAndSolver()->GetDofSet(), rA, rDx, rb);

        KRATOS_CATCH("");
    }

    /**
     * @brief This operations should be called before printing the results when non trivial results
     * (e.g. stresses)
     * Need to be calculated given the solution of the step
     * @details This operations should be called only when needed, before printing as it can involve a non
     * negligible cost
     */
    void CalculateOutputData() override
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        GetScheme()->CalculateOutputData(BaseType::GetModelPart(),
                                         GetBuilderAndSolver()->GetDofSet(),
                                         rA, rDx, rb);
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        KRATOS_TRY;

        // Pointers needed in the solution
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        ModelPart& r_model_part = BaseType::GetModelPart();

		//r_model_part.GetProcessInfo().SetValue(INACCURATE_PLASTIC_POINTS, 0);

        // Set up the system, operation performed just once unless it is required
        // to reform the dof set at each iteration
        BuiltinTimer system_construction_time;
        if (!p_builder_and_solver->GetDofSetIsInitializedFlag() || mReformDofSetAtEachStep)
        {
            // Setting up the list of the DOFs to be solved
            BuiltinTimer setup_dofs_time;
            p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", BaseType::GetEchoLevel() > 0) << "Setup Dofs Time: "
                << setup_dofs_time << std::endl;

            // Shaping correctly the system
            BuiltinTimer setup_system_time;
            p_builder_and_solver->SetUpSystem(r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", BaseType::GetEchoLevel() > 0) << "Setup System Time: "
                << setup_system_time << std::endl;

            // Setting up the Vectors involved to the correct size
            BuiltinTimer system_matrix_resize_time;
            p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, mpA, mpDx, mpb,
                                                                r_model_part);
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", BaseType::GetEchoLevel() > 0) << "System Matrix Resize Time: "
                << system_matrix_resize_time << std::endl;

            if (mReadForce) {
                const auto file_name = "incremental_solver_serialization";
                //Serializer::TraceType const& rTrace = Serializer::TraceType ::SERIALIZER_NO_TRACE;
                auto serializer = FileSerializer(file_name, Serializer::TraceType::SERIALIZER_NO_TRACE, false);
                load(serializer);
            }
            else
            {
                mPreviousExternalForceVector.resize(mpb->size(), false);
                TSparseSpace::SetToZero(mPreviousExternalForceVector);

            }


        }

        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", BaseType::GetEchoLevel() > 0) << "System Construction Time: "
            << system_construction_time << std::endl;

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        // Initial operations ... things that are constant over the Solution Step
        p_builder_and_solver->InitializeSolutionStep(r_model_part, rA, rDx, rb);

        // Initial operations ... things that are constant over the Solution Step
        p_scheme->InitializeSolutionStep(r_model_part, rA, rDx, rb);

        // Initialisation of the convergence criteria
        if (mpConvergenceCriteria->GetActualizeRHSflag())
        {
            mIniExternalForceVector.resize(rb.size(), false);
            TSparseSpace::SetToZero(mIniExternalForceVector);
            BuildExternalForceVector(r_model_part, mIniExternalForceVector);

            TSparseSpace::SetToZero(rb);
            double load_fraction = GetCurrentLoadFraction(r_model_part.GetProcessInfo());
            TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
            TSystemVectorType int_force = ZeroVector(rb.size());
            BuildInternalForceVector(r_model_part, int_force);
            rb = load_fraction * d_ext_force + mPreviousExternalForceVector - int_force;
			//rb = mIniExternalForceVector - int_force;
            //p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
        }

        mpConvergenceCriteria->InitializeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        if (mpConvergenceCriteria->GetActualizeRHSflag()) {
            TSparseSpace::SetToZero(rb);
        }

        KRATOS_CATCH("");
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        ModelPart& r_model_part = BaseType::GetModelPart();

        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        //Finalisation of the solution step,
        //operations to be done after achieving convergence, for example the
        //Final Residual Vector (mb) has to be saved in there
        //to avoid error accumulation

        p_scheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        p_builder_and_solver->FinalizeSolutionStep(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeSolutionStep(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        // Arc-length: adapt arc-length and promote step state on success.
        if (mUseArcLength && mLastStepConverged) {
            AdaptArcLength();
            OnArcLengthConverged();
            mLastStepConverged = false;
        }

        //Cleaning memory after the solution
        p_scheme->Clean();

        if (mReformDofSetAtEachStep == true) //deallocate the systemvectors
        {
            this->Clear();
        }
        KRATOS_CATCH("");
    }

    /**
    * @brief Gets the current solution vector
    * @param rDofSet The set of degrees of freedom (DOFs)
    * @param this_solution The vector that will be filled with the current solution values
    * @details This method retrieves the current solution values for the provided DOF set.
    * The provided solution vector will be resized to match the size of the DOF set if necessary,
    * and will be filled with the solution values corresponding to each DOF. Each value is accessed
    * using the equation ID associated with each DOF.
    */
    void GetCurrentSolution(DofsArrayType& rDofSet, Vector& this_solution) {
        this_solution.resize(rDofSet.size());
        block_for_each(rDofSet, [&](const auto& r_dof) {
            this_solution[r_dof.EquationId()] = r_dof.GetSolutionStepValue();
        });
    }

    /**
     * @brief Gets the matrix of non-converged solutions and the DOF set
     * @return A tuple containing the matrix of non-converged solutions and the DOF set
     * @details This method returns a tuple where the first element is a matrix of non-converged solutions and the second element is the corresponding DOF set.
     */
    std::tuple<Matrix, DofsArrayType> GetNonconvergedSolutions(){
        return std::make_tuple(mNonconvergedSolutionsMatrix, GetBuilderAndSolver()->GetDofSet());
    }

    /**
    * @brief Sets the state for storing non-converged solutions.
    * @param state The state to set for storing non-converged solutions (true to enable, false to disable).
    * @details This method enables or disables the storage of non-converged solutions at each iteration
    * by setting the appropriate flag. When the flag is set to true, non-converged solutions will be stored.
    */
    void SetUpNonconvergedSolutionsFlag(bool state) {
        mStoreNonconvergedSolutionsFlag = state;
    }

    void AssembleRHS(
        TSystemVectorType& b,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId
    )
    {
        unsigned int local_size = RHS_Contribution.size();

        for (unsigned int i_local = 0; i_local < local_size; i_local++) {
            unsigned int i_global = EquationId[i_local];

            // ASSEMBLING THE SYSTEM VECTOR
            double& b_value = b[i_global];
            const double& rhs_value = RHS_Contribution[i_local];

            AtomicAdd(b_value, rhs_value);
        }
    }

    void BuildInternalForceVector(ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();


        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements

        const int nelements = static_cast<int>(pElements.size());
#pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
#pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < nelements; i++) {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                // If the element is active
                if (it->IsActive()) {
                    //calculate elemental Right Hand Side Contribution
                    it->Calculate(INTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")
    }


    void BuildExternalForceVector(
        ModelPart& rModelPart,
        TSystemVectorType& b)
    {
        KRATOS_TRY

        //Getting the Elements
        ElementsArrayType& pElements = rModelPart.Elements();

        //getting the array of the conditions
        ConditionsArrayType& ConditionsArray = rModelPart.Conditions();

        const ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

        //contributions to the system
        LocalSystemVectorType RHS_Contribution = LocalSystemVectorType(0);

        //vector containing the localization in the system of the different
        //terms
        Element::EquationIdVectorType EquationId;

        // assemble all elements
        //for (typename ElementsArrayType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)

        const int nelements = static_cast<int>(pElements.size());
#pragma omp parallel firstprivate(nelements, RHS_Contribution, EquationId)
        {
#pragma omp for schedule(guided, 512) nowait
            for (int i = 0; i < nelements; i++) {
                typename ElementsArrayType::iterator it = pElements.begin() + i;
                // If the element is active
                if (it->IsActive()) {
                    //calculate elemental Right Hand Side Contribution
                    it->Calculate(EXTERNAL_FORCES_VECTOR, RHS_Contribution, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }

            RHS_Contribution.resize(0, false);

            // assemble all conditions
            const int nconditions = static_cast<int>(ConditionsArray.size());
#pragma omp for schedule(guided, 512)
            for (int i = 0; i < nconditions; i++) {
                auto it = ConditionsArray.begin() + i;
                // If the condition is active
                if (it->IsActive()) {
                    //calculate elemental contribution
                    it->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
                    it->EquationIdVector(EquationId, CurrentProcessInfo);

                    //pScheme->CalculateRHSContribution(*it, RHS_Contribution, EquationId, CurrentProcessInfo);

                    //assemble the elemental contribution
                    AssembleRHS(b, RHS_Contribution, EquationId);
                }
            }
        }

        KRATOS_CATCH("")

    }


    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        // Arc-length dispatch: when enabled, run Crisfield arc-length driver instead of plain Newton-Raphson.
        if (mUseArcLength) {
            return SolveArcLengthSolutionStep();
        }
        // Broyden / Woodbury quasi-Newton dispatch: K0 is factored once per step and reused
        // with rank-1 Broyden updates applied via the Sherman-Morrison-Woodbury identity.
        if (mUseBroyden) {
            return SolveBroydenSolutionStep();
        }

        double damping_factor = 1.0;
        // Pointers needed in the solution
        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();
        std::vector<Vector> NonconvergedSolutions;

        if (mStoreNonconvergedSolutionsFlag) {
            Vector initial;
            GetCurrentSolution(r_dof_set,initial);
            NonconvergedSolutions.push_back(initial);
        }

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        //initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool residual_is_updated = false;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

		mIniExternalForceVector.resize(rb.size(), false);
        TSparseSpace::SetToZero(mIniExternalForceVector);

		BuildExternalForceVector(r_model_part, mIniExternalForceVector);

		mIniFNorm = TSparseSpace::TwoNorm(mIniExternalForceVector);

		double rb_norm = TSparseSpace::TwoNorm(rb);
        double abs_tol = 1e-9;
        double rel_tol = 1e-2;
        bool is_converged = true;
		int max_inaccurate_plastic_points = 0;
        double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
        //mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        //bool is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

        // Function to perform the building and the solving phase.
        if (BaseType::mRebuildLevel > 0 || BaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            if (mUseOldStiffnessInFirstIteration){
                p_builder_and_solver->BuildAndSolveLinearizedOnPreviousIteration(p_scheme, r_model_part, rA, rDx, rb,BaseType::MoveMeshFlag());
            } else {

                p_builder_and_solver->BuildLHS(p_scheme, r_model_part, rA);
                if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol)) {
                    return false;
                }
                r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
                r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;
				//p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                //rb = load_fraction * rb;
                this->BuildRHS(r_model_part, r_dof_set, rb);
                
                p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);
                std::cout << "rb_norm: " << TSparseSpace::TwoNorm(rb) << " prev_force_vector_norm: " << TSparseSpace::TwoNorm(mPreviousExternalForceVector) << " ini_force_vector_norm: " << TSparseSpace::TwoNorm(mIniExternalForceVector) << std::endl;

                auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
                p_linear_solver->InitializeSolutionStep(rA, rDx, rb);
                p_linear_solver->PerformSolutionStep(rA, rDx, rb);
                //p_linear_solver->FinalizeSolutionStep(rA, rDx, rb);
                //p_builder_and_solver->SystemSolve(rA, rDx, rb);
            }
        } else {
            TSparseSpace::SetToZero(rDx);  // Dx = 0.00;
            TSparseSpace::SetToZero(rb);

            //p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            //rb = load_fraction * rb;
            this->BuildRHS(r_model_part, r_dof_set, rb);
            if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol)) {
                return false;
            }

            auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
            p_linear_solver->PerformSolutionStep(rA, rDx, rb);
            //p_linear_solver->FinalizeSolutionStep(rA, rDx, rb);

            //p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }
        BaseType::mStiffnessMatrixIsBuilt = true;

        // Debugging info
        EchoInfo(iteration_number);

        // Updating the results stored in the database
        UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        if (mStoreNonconvergedSolutionsFlag) {
            Vector first;
            GetCurrentSolution(r_dof_set,first);
            NonconvergedSolutions.push_back(first);
        }

        if (is_converged) {
            //if (mpConvergenceCriteria->GetActualizeRHSflag()) {
            //    TSparseSpace::SetToZero(rb);


            //    double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
            //    TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
            //    TSystemVectorType int_force = ZeroVector(rb.size());
            //    BuildInternalForceVector(r_model_part, int_force);
            //    rb = load_fraction * d_ext_force + mPreviousExternalForceVector - int_force;

            //    //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
            //    block_for_each(p_builder_and_solver->GetDofSet(), [&](Dof<double>& rDof) {
            //        const std::size_t i = rDof.EquationId();

            //        if (rDof.IsFixed())
            //            rb[i] = 0.0;
            //        });
            //    //p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            //}

			rb_norm = TSparseSpace::TwoNorm(rb);

            is_converged = rb_norm / mIniFNorm < rel_tol || rb_norm < abs_tol;
            std::cout << "rb_norm: " << rb_norm << "; IniFNorm: " << mIniFNorm << "; rel_tol: " << rel_tol << "; abs_tol: " << abs_tol << "; is_converged: " << is_converged << std::endl;

            //is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
            if (is_converged)
            {
                damping_factor = 1.0;
            }

        }

        //Iteration Cycle... performed only for NonLinearProblems
        while (is_converged == false &&
               iteration_number++ < mMaxIterationNumber)
        {
            //setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
			r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

            mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            is_converged = true;
            //rb_norm = TSparseSpace::TwoNorm(rb);
            //is_converged = rb_norm / mIniFNorm < rel_tol || rb_norm < abs_tol;
            //std::cout << "rb_norm: " << rb_norm << "; IniFNorm: " << mIniFNorm << "; rel_tol: " << rel_tol << "; abs_tol: " << abs_tol << "; is_converged: " << is_converged << std::endl;

            //is_converged = mpConvergenceCriteria->PreCriteria(r_model_part, r_dof_set, rA, rDx, rb);

            //call the linear system solver to find the correction mDx for the
            //it is not called if there is no system to solve
			//std::cout << "rebuild level; " << "stiffness matrix is built; " << "keep system constant during iterations: " << std::endl;
			//std::cout << BaseType::mRebuildLevel << " " << BaseType::mStiffnessMatrixIsBuilt << " " << GetKeepSystemConstantDuringIterations() << std::endl;
            if (SparseSpaceType::Size(rDx) != 0) 
            {
                if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
                {
                    if (GetKeepSystemConstantDuringIterations() == false)
                    {
                        const auto timer = BuiltinTimer();

                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

						std::cout << "Building LHS at iteration " << iteration_number << std::endl;
                        p_builder_and_solver->BuildLHS(p_scheme, r_model_part, rA);
                        if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol)) {
                            return false;
                        }
                        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
                        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;
						//p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                        //rb = load_fraction * rb;
                        this->BuildRHS(r_model_part, r_dof_set, rb);

                        p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

						std::cout << "rb_norm: " << TSparseSpace::TwoNorm(rb) << " prev_force_vector_norm: " << TSparseSpace::TwoNorm(mPreviousExternalForceVector) << " ini_force_vector_norm: " << TSparseSpace::TwoNorm(mIniExternalForceVector) << std::endl;

                        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
						p_linear_solver->InitializeSolutionStep(rA, rDx, rb);
                        p_linear_solver->PerformSolutionStep(rA, rDx, rb);



                    }
                    else
                    {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        this->BuildRHS(r_model_part, r_dof_set, rb);
						//p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                        //rb = load_fraction * rb;
                        if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol)) {
                            return false;
                        }

                        std::cout << "rb_norm: " << TSparseSpace::TwoNorm(rb) << " prev_force_vector_norm: " << TSparseSpace::TwoNorm(mPreviousExternalForceVector) << " ini_force_vector_norm: " << TSparseSpace::TwoNorm(mIniExternalForceVector) << std::endl;

                        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
                        //p_linear_solver->InitializeSolutionStep(rA, rDx, rb);
                        p_linear_solver->PerformSolutionStep(rA, rDx, rb);
                        //p_builder_and_solver->SystemSolve(rA, rDx, rb);

                        //p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                }
                else
                {
					//std::cout << "ATTENTION: the system is not being rebuilt, but the flag to keep it constant is set to false. This can lead to unexpected results. Please check the settings of the strategy." << std::endl;
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    //p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    //rb = load_fraction * rb;
                    this->BuildRHS(r_model_part, r_dof_set, rb);
                    if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol)) {
                        return false;
                    }
                   
                    std::cout << "rb_norm: " << TSparseSpace::TwoNorm(rb) << " prev_force_vector_norm: " << TSparseSpace::TwoNorm(mPreviousExternalForceVector) << " ini_force_vector_norm: " << TSparseSpace::TwoNorm(mIniExternalForceVector) << std::endl;

                    auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();
                    //p_linear_solver->InitializeSolutionStep(rA, rDx, rb);
                    p_linear_solver->PerformSolutionStep(rA, rDx, rb);
                    //p_builder_and_solver->SystemSolve(rA, rDx, rb);

                }
            }
            else
            {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            EchoInfo(iteration_number);

            TSystemVectorType damped_rDx(rDx.size());
			TSparseSpace::Copy(rDx, damped_rDx);
            TSparseSpace::InplaceMult(damped_rDx, damping_factor);
            UpdateDatabase(rA, damped_rDx, rb, BaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            if (mStoreNonconvergedSolutionsFlag == true){
                Vector ith;
                GetCurrentSolution(r_dof_set,ith);
                NonconvergedSolutions.push_back(ith);
            }

            residual_is_updated = false;

            if (is_converged == true)
            {
				//std::cout << "GetActualizedRHS: " << mpConvergenceCriteria->GetActualizeRHSflag() << std::endl;
    //            if (mpConvergenceCriteria->GetActualizeRHSflag() == true)
    //            {
    //                TSparseSpace::SetToZero(rb);


    //                double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
    //                TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
    //                TSystemVectorType int_force = ZeroVector(rb.size());
    //                BuildInternalForceVector(r_model_part, int_force);
    //                rb = load_fraction * d_ext_force + mPreviousExternalForceVector - int_force;

    //                //NOTE: dofs are assumed to be numbered consecutively in the BlockBuilderAndSolver
    //                block_for_each(p_builder_and_solver->GetDofSet(), [&](Dof<double>& rDof) {
    //                    const std::size_t i = rDof.EquationId();

    //                    if (rDof.IsFixed())
    //                        rb[i] = 0.0;
    //                    });

    //                //p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
    //                residual_is_updated = true;
    //            }
                rb_norm = TSparseSpace::TwoNorm(rb);
                is_converged = rb_norm / mIniFNorm < rel_tol || rb_norm < abs_tol;
                std::cout << "rb_norm: " << rb_norm << "; IniFNorm: " << mIniFNorm << "; rel_tol: " << rel_tol << "; abs_tol: " << abs_tol << "; is_converged: " << is_converged << std::endl;

                //is_converged = mpConvergenceCriteria->PostCriteria(r_model_part, r_dof_set, rA, rDx, rb);
				/*if (is_converged)
				{
                    current_damping_step += 1;
                    factor_calculated = damping_factor * ;
				}*/
            }
        }
        //plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= mMaxIterationNumber) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number << " / "
                << mMaxIterationNumber << " iterations" << std::endl;
        }
        //recalculate residual if needed
        //(note that some convergence criteria need it to be recalculated)
        if (residual_is_updated == false)
        {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

            //    TSparseSpace::SetToZero(mb);
            //    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        //calculate reactions if required
        if (mCalculateReactionsFlag == true)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
        if (mStoreNonconvergedSolutionsFlag) {
            mNonconvergedSolutionsMatrix = Matrix( r_dof_set.size(), NonconvergedSolutions.size() );
            for (std::size_t i = 0; i < NonconvergedSolutions.size(); ++i) {
                block_for_each(r_dof_set, [&](const auto& r_dof) {
                    mNonconvergedSolutionsMatrix(r_dof.EquationId(), i) = NonconvergedSolutions[i](r_dof.EquationId());
                });
            }
        }
        return is_converged;
    }

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetBuilderAndSolver()->Check(BaseType::GetModelPart());

        GetScheme()->Check(BaseType::GetModelPart());

        mpConvergenceCriteria->Check(BaseType::GetModelPart());

        return 0;

        KRATOS_CATCH("")
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name"                                : "newton_raphson_strategy",
            "use_old_stiffness_in_first_iteration": false,
            "max_iteration"                       : 10,
            "reform_dofs_at_each_step"            : false,
            "compute_reactions"                   : false,
            "use_arc_length"                      : false,
            "use_broyden"                         : false,
            "arc_length_settings"                 : {
                "initial_arc_length"        : 0.0,
                "max_arc_length_factor"     : 10.0,
                "min_arc_length_factor"     : 0.1,
                "desired_iterations"        : 4,
                "beta"                      : 0.0,
                "max_load_factor"           : 1.0,
                "max_arc_length_reductions" : 5
            },
            "builder_and_solver_settings"         : {},
            "convergence_criteria_settings"       : {},
            "linear_solver_settings"              : {},
            "scheme_settings"                     : {}
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "newton_raphson_strategy";
    }

    ///@}
    ///@name Operators

    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief This method returns the LHS matrix
     * @return The LHS matrix
     */
    TSystemMatrixType &GetSystemMatrix() override
    {
        TSystemMatrixType &mA = *mpA;

        return mA;
    }

    /**
     * @brief This method returns the RHS vector
     * @return The RHS vector
     */
    TSystemVectorType& GetSystemVector() override
    {
        TSystemVectorType& mb = *mpb;

        return mb;
    }

    /**
     * @brief This method returns the solution vector
     * @return The Dx vector
     */
    TSystemVectorType& GetSolutionVector() override
    {
        TSystemVectorType& mDx = *mpDx;

        return mDx;
    }

    /**
     * @brief Set method for the flag mKeepSystemConstantDuringIterations
     * @param Value If we consider constant the system of equations during the iterations
     */
    void SetKeepSystemConstantDuringIterations(bool Value)
    {
        mKeepSystemConstantDuringIterations = Value;
    }

    /**
     * @brief Enable / disable arc-length control (Crisfield's spherical/cylindrical method).
     *        When enabled, SolveSolutionStep() dispatches to the arc-length driver and
     *        the externally-supplied load fraction (derived from TIME) is replaced by the
     *        internal load factor lambda controlled by the arc-length constraint.
     */
    void SetUseArcLength(bool Value) { mUseArcLength = Value; }
    bool GetUseArcLength() const     { return mUseArcLength; }

    /**
     * @brief Enable / disable the Broyden + Woodbury quasi-Newton solver.
     *        When enabled, the tangent K0 is built and factored only ONCE per solution step
     *        (via p_linear_solver->InitializeSolutionStep), and subsequent iterations reuse
     *        that factorization (p_linear_solver->PerformSolutionStep) together with rank-1
     *        Broyden updates resolved through the Sherman-Morrison-Woodbury identity.
     */
    void SetUseBroyden(bool Value) { mUseBroyden = Value; }
    bool GetUseBroyden() const     { return mUseBroyden; }

    void SetBroydenSettings(
        double AbsoluteGlobalError,
        double RelativeGlobalError,
        double RelativeLocalError)
    {
		mMaxGlobalAbsoluteError = AbsoluteGlobalError;
		mMaxGlobalRelativeError = RelativeGlobalError;
		mMaxLocalRelativeError = RelativeLocalError;
    }

    /**
     * @brief Configure arc-length parameters at runtime.
     * @param InitialArcLength       Reference arc-length. If <= 0, it is auto-computed from the first predictor.
     * @param MaxArcLengthFactor     Upper bound multiplier: max radius = factor * InitialArcLength.
     * @param MinArcLengthFactor     Lower bound multiplier: min radius = factor * InitialArcLength.
     * @param DesiredIterations      Target iteration count used for arc-length adaptation.
     * @param Beta                   0 = cylindrical, 1 = spherical, in-between for mixed.
     * @param MaxLoadFactor          Termination cap for the cumulative load factor lambda.
     * @param MaxArcLengthReductions Max retries (halvings) when the constraint quadratic has no real root.
     */
    void SetArcLengthSettings(double InitialArcLength,
                              double MaxArcLengthFactor   = 10.0,
                              double MinArcLengthFactor   = 0.1,
                              unsigned int DesiredIterations = 4,
                              double Beta                 = 0.0,
                              double MaxLoadFactor        = 1.0,
                              unsigned int MaxArcLengthReductions = 5)
    {
        mArcLength0             = InitialArcLength;
        mArcLength              = InitialArcLength;
        mMaxArcLengthFactor     = MaxArcLengthFactor;
        mMinArcLengthFactor     = MinArcLengthFactor;
        mDesiredIterations      = DesiredIterations;
        mBeta                   = Beta;
        mLambdaMax              = MaxLoadFactor;
        mMaxArcLengthReductions = MaxArcLengthReductions;
    }

    /// Returns the current converged load factor (only meaningful when arc-length is enabled).
    double GetArcLengthLambda() const { return mLambda; }

    /**
     * @brief Get method for the flag mKeepSystemConstantDuringIterations
     * @return True if we consider constant the system of equations during the iterations, false otherwise
     */
    bool GetKeepSystemConstantDuringIterations()
    {
        return mKeepSystemConstantDuringIterations;
    }

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ResidualBasedNewtonRaphsonStrategyTwo";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
  protected:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    typename TSchemeType::Pointer mpScheme = nullptr; /// The pointer to the time scheme employed
    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver = nullptr; /// The pointer to the builder and solver employed
    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria = nullptr; /// The pointer to the convergence criteria employed

    TSystemVectorPointerType mpDx; /// The increment in the solution
    TSystemVectorPointerType mpb; /// The RHS vector of the system of equations
    TSystemVectorType mIniExternalForceVector;
    TSystemVectorType mPreviousExternalForceVector;
    TSystemMatrixPointerType mpA; /// The LHS matrix of the system of equations
   
    bool mReadForce = false; /// A flag to set if the force is read from file or not

	double mIniFNorm; /// The initial residual norm, used for convergence criteria that require it

    // ===== Arc-length (Crisfield) state =====
    bool   mUseArcLength            = false;  ///< If true, SolveSolutionStep uses arc-length control
    double mBeta                    = 0.0;    ///< 0 = cylindrical, 1 = spherical
    double mArcLength               = 0.0;    ///< Current arc-length radius
    double mArcLength0              = 0.0;    ///< Reference arc-length (computed from first predictor if <= 0)
    double mMaxArcLengthFactor      = 10.0;   ///< Upper bound: mArcLength <= mMaxArcLengthFactor * mArcLength0
    double mMinArcLengthFactor      = 0.1;    ///< Lower bound: mArcLength >= mMinArcLengthFactor * mArcLength0
    unsigned int mDesiredIterations = 4;      ///< Target iterations per step for arc-length adaptation
    double mLambda                  = 0.0;    ///< Converged load factor at the start of the current step
    double mDeltaLambda             = 0.0;    ///< Load-factor increment for the current step
    double mLambdaMax               = 1.0;    ///< Cap on (mLambda + mDeltaLambda)
    unsigned int mMaxArcLengthReductions = 5; ///< Max number of arc-length halvings on snap-back failure
    int    mPredictorSignHint       = 1;      ///< Sign hint for predictor (carried between steps)
    TSystemVectorType mDeltaU;                ///< Accumulated displacement increment over the current step
    TSystemVectorType mQ;                     ///< Reference load vector for the step (= mIniExternalForceVector - mPreviousExternalForceVector)
    bool   mLastStepConverged       = false;  ///< Set in SolveSolutionStep / used in FinalizeSolutionStep for adaptation
    unsigned int mLastStepIterations = 1;     ///< Iteration count of the last converged arc-length step

    // ===== Broyden / Woodbury quasi-Newton state =====
    bool mUseBroyden = true; ///< If true, SolveSolutionStep uses the Broyden + Woodbury solver
    double mMaxGlobalAbsoluteError = 1e-12;
    double mMaxGlobalRelativeError = 1e-3;
    double mMaxLocalRelativeError = 1e-1;

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    * @details Default = false
        - true  : Reforme at each time step
        - false : Form just one (more efficient)
     */
    bool mReformDofSetAtEachStep;

    /**
     * @brief Flag telling if it is needed or not to compute the reactions
     * @details default = true
     */
    bool mCalculateReactionsFlag;

    /**
     * @brief Flag telling if a full update of the database will be performed at the first iteration
     * @details default = false
     */
    bool mUseOldStiffnessInFirstIteration = false;

    unsigned int mMaxIterationNumber; /// The maximum number of iterations, 30 by default

    bool mInitializeWasPerformed; /// Flag to set as initialized the strategy

    bool mKeepSystemConstantDuringIterations; // Flag to allow keeping system matrix constant during iterations

    /**
     * @brief This matrix stores the non-converged solutions
     * @details The matrix is structured such that each column represents the solution vector at a specific non-converged iteration.
     */
    Matrix mNonconvergedSolutionsMatrix;

    /**
     * @brief Flag indicating whether to store non-converged solutions
     * @details Only when set to true (by calling the SetUpNonconvergedSolutionsGathering method) will the non-converged solutions at each iteration be stored.
     */
    bool mStoreNonconvergedSolutionsFlag = false;

    bool LocalConvergenceAchieved(ProcessInfo& rProcessInfo, const double relative_tolerance)
    {
        double max_inaccurate_plastic_points = relative_tolerance * rProcessInfo[N_PLASTIC_POINTS] + 3;

        if (rProcessInfo[INACCURATE_PLASTIC_POINTS] > max_inaccurate_plastic_points)
        {
            std::cout << "Warning: number of inaccurate plastic points: " << rProcessInfo[INACCURATE_PLASTIC_POINTS] << ", max allowed: " << max_inaccurate_plastic_points << std::endl;
            return false;
        }
		return true;
    }

	void BuildRHS(
		ModelPart& r_model_part,
        DofsArrayType& r_dof_set,
		TSystemVectorType& rb)
	{
        const double load_fraction = GetCurrentLoadFraction(r_model_part.GetProcessInfo());

        TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
        
        TSystemVectorType int_force = ZeroVector(rb.size());

        // reset plastic points
        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

        BuildInternalForceVector(r_model_part, int_force);
        
        TSystemVectorType ini_int_force = ZeroVector(rb.size());

		// initial internal force is equal to external force from previous step, as the system is in equilibrium at the beginning of the step
        TSparseSpace::Copy(mPreviousExternalForceVector, ini_int_force);

		TSystemVectorType d_int_force = ZeroVector(rb.size());
		TSparseSpace::ScaleAndAdd(1.0, int_force, -1.0, ini_int_force, d_int_force); // ini_int_force = mPreviousExternalForceVector - int_force;
        //TSystemVectorType d_int_force = int_force - ini_int_force;

        // norm int force
        //double int_force_norm = TSparseSpace::TwoNorm(int_force);
        //std::cout << "int_force_norm: " << int_force_norm << std::endl;

		// load fraction * (d_Fext - d_Fint) + Fext_prev - Fint_prev
		TSparseSpace::ScaleAndAdd(load_fraction, d_ext_force, -load_fraction, d_int_force, rb); //
        //rb = load_fraction * (d_ext_force - d_int_force);
		//rb = load_fraction * d_ext_force + mPreviousExternalForceVector - int_force;
        //rb = (d_ext_force - d_int_force);

		// apply Dirichlet conditions RHS
        //NOTE: dofs are assumed to be numbered consecutively
        block_for_each(r_dof_set, [&](Dof<double>& rDof) {
            const std::size_t i = rDof.EquationId();

            if (rDof.IsFixed())
                rb[i] = 0.0;
            });
	}

    /**
     * @brief Returns the active load fraction used to scale the external-load increment in BuildRHS.
     *        - Standard mode: linear in time, (TIME - START_TIME) / (END_TIME - START_TIME).
     *        - Arc-length mode: (mLambda + mDeltaLambda).
     */
    double GetCurrentLoadFraction(const ProcessInfo& rProcessInfo) const
    {
        if (mUseArcLength) {
            return mLambda + mDeltaLambda;
        }
        const double t0 = rProcessInfo[START_TIME];
        const double t1 = rProcessInfo[END_TIME];
        const double t  = rProcessInfo[TIME];
        const double dt = t1 - t0;
        return (std::abs(dt) > 0.0) ? (t - t0) / dt : 1.0;
    }

    /**
     * @brief Apply Dirichlet zeroing on a vector b in-place (used for predictor / corrector RHS).
     */
    void ZeroOutFixedDofs(DofsArrayType& r_dof_set, TSystemVectorType& rb)
    {
        block_for_each(r_dof_set, [&](Dof<double>& rDof) {
            const std::size_t i = rDof.EquationId();
            if (rDof.IsFixed())
                rb[i] = 0.0;
        });
    }

    /**
     * @brief Top-level driver for one solution step using Crisfield arc-length control.
     *        Mirrors the non-linear loop structure of SolveSolutionStep() but uses
     *        the predictor / corrector defined below and adapts the load factor lambda.
     */
    bool SolveArcLengthSolutionStep()
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        auto p_scheme           = GetScheme();
        auto p_bs               = GetBuilderAndSolver();
        auto& r_dof_set         = p_bs->GetDofSet();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        const double abs_tol = 1e-9;
        const double rel_tol = 1e-2;
		double local_converged = true;
        double max_d_lambda = 1.0;

        // Build the reference load vector q for this step (target Fext step increment).
        mIniExternalForceVector.resize(rb.size(), false);
        TSparseSpace::SetToZero(mIniExternalForceVector);
        BuildExternalForceVector(r_model_part, mIniExternalForceVector);
        mIniFNorm = TSparseSpace::TwoNorm(mIniExternalForceVector);

        mQ.resize(rb.size(), false);
        TSparseSpace::SetToZero(mQ);
        noalias(mQ) = mIniExternalForceVector - mPreviousExternalForceVector;
        const double q_norm = TSparseSpace::TwoNorm(mQ);

        // Edge case: no load increment this step -> degrade to plain Newton-Raphson.
        if (q_norm < std::numeric_limits<double>::epsilon()) {
            KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategyTwo")
                << "Arc-length: ||q||=0 for this step, falling back to load-control Newton-Raphson." << std::endl;
            const bool saved = mUseArcLength;
            mUseArcLength = false;
            const bool ok = this->SolveSolutionStep();
            mUseArcLength = saved;
            return ok;
        }

        // Reset step-local state and allocate Δu.
        mDeltaLambda = 0.0;
        mDeltaU.resize(rb.size(), false);
        TSparseSpace::SetToZero(mDeltaU);

        // Retry loop: if the constraint quadratic has no real root, halve the arc-length and retry.
        bool   step_converged = false;
        unsigned int iterations_used = 0;

        for (unsigned int retry = 0; retry <= mMaxArcLengthReductions && !step_converged; ++retry) {

            // Reset step-local state at each retry.
            mDeltaLambda = 0.0;
            TSparseSpace::SetToZero(mDeltaU);

            // ---- Predictor ----
            TSystemVectorType delta_ut(rb.size());
            TSparseSpace::SetToZero(delta_ut);

            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = 0;
            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
            r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

            // Build tangent K (predictor): honor mRebuildLevel / mStiffnessMatrixIsBuilt.
            // - mRebuildLevel == 0 : never rebuild after first build (only if not yet built)
            // - mRebuildLevel >= 1 : always rebuild at the predictor of every step
            const bool build_predictor_lhs =
                (BaseType::mRebuildLevel >= 1) || (BaseType::mStiffnessMatrixIsBuilt == false);
            if (build_predictor_lhs) {
                TSparseSpace::SetToZero(rA);
                p_bs->BuildLHS(p_scheme, r_model_part, rA);
                local_converged = LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol);
				if (!local_converged) {
					return false;
				}

                p_bs->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);
                BaseType::mStiffnessMatrixIsBuilt = true;
            }

            // Solve K * delta_ut = q
            TSystemVectorType q_dir = mQ;
            ZeroOutFixedDofs(r_dof_set, q_dir);
            auto p_linear_solver = p_bs->GetLinearSystemSolver();
            if (build_predictor_lhs) {
                p_linear_solver->InitializeSolutionStep(rA, delta_ut, q_dir);
            }
            p_linear_solver->PerformSolutionStep(rA, delta_ut, q_dir);

            const double ut_dot_ut = TSparseSpace::Dot(delta_ut, delta_ut);
            const double denom = std::sqrt(ut_dot_ut + mBeta * mBeta * q_norm * q_norm);

            // Auto-initialise arc-length on very first arc-length step if user provided <= 0.
            if (mArcLength0 <= 0.0) {
                mArcLength0 = denom;            // size of the very first predictor
                mArcLength  = mArcLength0;
            }

            // Predictor magnitude and sign.
            double d_lambda = (denom > 0.0) ? (mArcLength / denom) : 0.0;
            d_lambda *= static_cast<double>(mPredictorSignHint);
			if (d_lambda > max_d_lambda) {
				d_lambda = max_d_lambda;
			}
			else if (d_lambda < -max_d_lambda) {
				d_lambda = -max_d_lambda;
			}

			std::cout << "d_lambda: " << d_lambda << "; ut_dot_ut: " << ut_dot_ut << "; q_norm: " << q_norm << "; denom: " << denom << std::endl;

            // Cap by mLambdaMax: do not overshoot the final target load factor.
            if (mLambda + d_lambda > mLambdaMax) {
                d_lambda = mLambdaMax - mLambda;
            }

            mDeltaLambda = d_lambda;
            TSystemVectorType delta_u_predictor(rb.size());
            TSparseSpace::Assign(delta_u_predictor, d_lambda, delta_ut);
            noalias(mDeltaU) = delta_u_predictor;

            // Apply predictor displacement to the DOF database.
            UpdateDatabase(rA, delta_u_predictor, rb, BaseType::MoveMeshFlag());

            // ---- Corrector iterations ----
            unsigned int iteration_number = 1;
            bool   converged = false;
            bool   discriminant_failure = false;

            while (!converged && iteration_number <= mMaxIterationNumber) {
                r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
                p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
                r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
                r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;
                mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

                // Rebuild tangent (corrector): honor mRebuildLevel and mKeepSystemConstantDuringIterations.
                // - mRebuildLevel  > 1 AND !KeepSystemConstant : rebuild each corrector iteration (full Newton)
                // - mRebuildLevel == 1                         : keep K from predictor (modified Newton within the step)
                // - mRebuildLevel == 0                         : never rebuild (initial-stiffness method)
                const bool rebuild_corrector_lhs =
                    (BaseType::mRebuildLevel > 1) && (mKeepSystemConstantDuringIterations == false);
                if (rebuild_corrector_lhs) {
                    TSparseSpace::SetToZero(rA);
                    p_bs->BuildLHS(p_scheme, r_model_part, rA);
                    local_converged = LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol);
					if (!local_converged) {
						return false;
					}

                    p_bs->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);
                    BaseType::mStiffnessMatrixIsBuilt = true;
                }

                // Residual at current (u, lambda): rb = (mLambda+mDeltaLambda) * d_ext_force - d_int_force
                //                                    + previous_step_residual_terms
                // The existing BuildRHS already uses GetCurrentLoadFraction() so this is consistent.
                TSparseSpace::SetToZero(rb);
                this->BuildRHS(r_model_part, r_dof_set, rb);
                //local_converged = LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol);

                // Solve K * delta_ubar = rb (residual correction)
                TSystemVectorType delta_ubar(rb.size());
                TSparseSpace::SetToZero(delta_ubar);
                if (rebuild_corrector_lhs) {
                    p_linear_solver->InitializeSolutionStep(rA, delta_ubar, rb);
                }
                p_linear_solver->PerformSolutionStep(rA, delta_ubar, rb);

                // Solve K * delta_ut = q (tangent direction)
                TSparseSpace::SetToZero(delta_ut);
                TSystemVectorType q_dir2 = mQ;
                ZeroOutFixedDofs(r_dof_set, q_dir2);
                p_linear_solver->PerformSolutionStep(rA, delta_ut, q_dir2);

                // Crisfield quadratic: a δλ² + b δλ + c = 0 with
                //   x = Δu + δubar
                //   a = δut·δut + β² ||q||²
                //   b = 2 [ x·δut + β² ΔΛ ||q||² ]
                //   c = x·x   + β² ΔΛ² ||q||² - L²
                TSystemVectorType x(rb.size());
                noalias(x) = mDeltaU + delta_ubar;

                const double a = TSparseSpace::Dot(delta_ut, delta_ut) + mBeta*mBeta*q_norm*q_norm;
                const double b = 2.0 * (TSparseSpace::Dot(x, delta_ut) + mBeta*mBeta*mDeltaLambda*q_norm*q_norm);
                const double c = TSparseSpace::Dot(x, x) + mBeta*mBeta*mDeltaLambda*mDeltaLambda*q_norm*q_norm
                                 - mArcLength*mArcLength;

                const double discriminant = b*b - 4.0*a*c;
                double d_lam = 0.0;
				std::cout << "a: " << a << "; b: " << b << "; c: " << c << "; discriminant: " << discriminant << std::endl;
                if (discriminant < 0.0 || a <= 0.0) {
                    // Snap-back / ill-conditioned step: bail out to the retry loop with smaller arc-length.
                    discriminant_failure = true;
                    break;
                } else {
                    const double sq = std::sqrt(discriminant);
                    const double r1 = (-b + sq) / (2.0 * a);
                    const double r2 = (-b - sq) / (2.0 * a);
                    // Crisfield angle criterion: pick the root maximizing (Δu + δu)·Δu
                    TSystemVectorType du1(rb.size());
                    TSystemVectorType du2(rb.size());
                    noalias(du1) = delta_ubar + r1 * delta_ut;
                    noalias(du2) = delta_ubar + r2 * delta_ut;
                    const double dot1 = TSparseSpace::Dot(mDeltaU + du1, mDeltaU);
                    const double dot2 = TSparseSpace::Dot(mDeltaU + du2, mDeltaU);
                    d_lam = (dot1 >= dot2) ? r1 : r2;

                    if (d_lam > max_d_lambda) {
                        d_lam = max_d_lambda;
                    }
                    else if (d_lam < -max_d_lambda) {
                        d_lam = -max_d_lambda;
                    }
                }

                // Cap by mLambdaMax
                if (mLambda + mDeltaLambda + d_lam > mLambdaMax) {
                    d_lam = mLambdaMax - (mLambda + mDeltaLambda);
                }

                // Total iterative correction: δu = δubar + δλ * δut
                TSystemVectorType delta_u(rb.size());
                noalias(delta_u) = delta_ubar + d_lam * delta_ut;

                // Update accumulated step state
                mDeltaLambda += d_lam;
                noalias(mDeltaU) += delta_u;

                // Apply to DOF database
                UpdateDatabase(rA, delta_u, rb, BaseType::MoveMeshFlag());

                p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
                mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
                EchoInfo(iteration_number);

                // Convergence: residual norm vs reference scaled by current load increment.
                TSparseSpace::SetToZero(rb);
                this->BuildRHS(r_model_part, r_dof_set, rb);
                local_converged = LocalConvergenceAchieved(r_model_part.GetProcessInfo(), rel_tol);

                const double rb_norm  = TSparseSpace::TwoNorm(rb);
                const double ref_norm = std::max(mIniFNorm, q_norm * std::abs(mDeltaLambda));
				converged = ((rb_norm < abs_tol) || (ref_norm > 0.0 && (rb_norm / ref_norm) < rel_tol)) && local_converged;

                KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo[ArcLength]", BaseType::GetEchoLevel() > 0)
                    << "iter " << iteration_number
                    << "  lambda=" << (mLambda + mDeltaLambda)
                    << "  ||rb||=" << rb_norm
                    << "  ref=" << ref_norm
                    << "  converged=" << converged << std::endl;

                if (!converged) ++iteration_number;
            }

            iterations_used = iteration_number;

            if (discriminant_failure || !converged) {
                // Roll back DOF database by reapplying -mDeltaU (so next retry starts from converged state).
                TSystemVectorType rollback(rb.size());
                TSparseSpace::Assign(rollback, -1.0, mDeltaU);
                UpdateDatabase(rA, rollback, rb, BaseType::MoveMeshFlag());

                // Halve the arc-length AND the per-step load-fraction cap in lockstep.
                // Without halving max_d_lambda, the predictor would remain clamped to the same
                // value across retries and trial N would reproduce trial 1 byte-for-byte.
                const double new_arc = std::max(mArcLength * 0.5, mMinArcLengthFactor * mArcLength0);
                const double scale   = (mArcLength > 0.0) ? (new_arc / mArcLength) : 0.5;
                mArcLength    = new_arc;
                max_d_lambda *= scale;

                // Force the predictor of the next retry to rebuild K at the rolled-back state,
                // regardless of mRebuildLevel. The DOF state has changed, so the cached K is stale.
                //BaseType::mStiffnessMatrixIsBuilt = false;

                KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategyTwo[ArcLength]")
                    << "Step failed (discriminant<0 or max-iter): halving arc-length to "
                    << mArcLength << ", max_d_lambda=" << max_d_lambda
                    << " (retry " << retry+1 << "/" << mMaxArcLengthReductions << ")" << std::endl;
                continue;
            }

            step_converged = true;
        }

        if (!step_converged) {
            MaxIterationsExceeded();
            return false;
        }

        // Bookkeeping for FinalizeSolutionStep + next-step predictor sign.
        mLastStepIterations = std::max<unsigned int>(iterations_used, 1u);
        mLastStepConverged  = true;
        mPredictorSignHint  = (mDeltaLambda >= 0.0) ? 1 : -1;

        // Calculate reactions if required.
        if (mCalculateReactionsFlag) {
            p_bs->CalculateReactions(p_scheme, r_model_part, rA, *mpDx, rb);
        }

        return true;

        KRATOS_CATCH("")
    }

    /// Arc-length adaptation: enlarges/shrinks the arc-length based on iterations used vs desired.
    void AdaptArcLength()
    {
        if (mArcLength0 <= 0.0) return;
        const double ratio = std::sqrt(static_cast<double>(mDesiredIterations) /
                                       static_cast<double>(std::max<unsigned int>(mLastStepIterations, 1u)));
        mArcLength = std::min(mMaxArcLengthFactor * mArcLength0,
                              std::max(mMinArcLengthFactor * mArcLength0, mArcLength * ratio));
    }

    /// Update the converged load factor / previous external force vector at the end of a converged step.
    void OnArcLengthConverged()
    {
        // Promote step into the running totals.
        mLambda += mDeltaLambda;
        if (mLambda > mLambdaMax) mLambda = mLambdaMax;
        // Update mPreviousExternalForceVector so that the *next* step's q = (new mIniExt - prev) is correct.
        // Since mIniExternalForceVector is rebuilt each step from element/condition contributions, we update
        // the "previously applied" external force to reflect the converged scaled load.
        if (mPreviousExternalForceVector.size() == mIniExternalForceVector.size()) {
            noalias(mPreviousExternalForceVector) = mPreviousExternalForceVector + mDeltaLambda * mQ;
        }
        mDeltaLambda = 0.0;
    }

    /**
     * @brief Solve a small dense k x k linear system M * x = b using
     *        Gaussian elimination with partial pivoting. Returns false if M is singular.
     *        Used for the Woodbury (I + V^T Z) capacitance system.
     *        Note: rM and rRhs are taken by value (mutated locally).
     */
    bool SolveSmallDense(Matrix rM, Vector rRhs, Vector& rSolution)
    {
        const std::size_t k = rRhs.size();
        rSolution.resize(k, false);
        if (k == 0) return true;

        // Forward elimination with partial pivoting
        for (std::size_t i = 0; i < k; ++i) {
            std::size_t pivot = i;
            double pivot_val = std::abs(rM(i, i));
            for (std::size_t r = i + 1; r < k; ++r) {
                const double v = std::abs(rM(r, i));
                if (v > pivot_val) { pivot_val = v; pivot = r; }
            }
            if (pivot_val < std::numeric_limits<double>::epsilon()) {
                return false; // singular
            }
            if (pivot != i) {
                for (std::size_t c = i; c < k; ++c) std::swap(rM(i, c), rM(pivot, c));
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
        // Back substitution
        for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(k) - 1; i >= 0; --i) {
            double s = rRhs[i];
            for (std::size_t c = i + 1; c < k; ++c) {
                s -= rM(i, c) * rSolution[c];
            }
            rSolution[i] = s / rM(i, i);
        }
        return true;
    }

    /**
     * @brief Broyden quasi-Newton solver using the Sherman-Morrison-Woodbury identity.
     *
     * The tangent K0 is built and LU-factored ONCE at the start of the step
     * (p_linear_solver->InitializeSolutionStep). All subsequent iterations reuse that
     * factorization (p_linear_solver->PerformSolutionStep) and apply accumulated rank-1
     * Broyden updates K_k = K0 + U V^T via:
     *
     *     (K0 + U V^T)^{-1} r = w - Z (I + V^T Z)^{-1} (V^T w),
     *     where w = K0^{-1} r and Z = K0^{-1} U.
     *
     * Each new iteration appends one column:
     *     s  = du,
     *     y  = R_new - R              (here R = -rb, so y_kratos = rb - rb_new),
     *     u  = y - K_k s,
     *     v  = s / (s . s),
     *     z  = K0^{-1} u              (sparse triangular solve, reused factorization).
     */
    bool SolveBroydenSolutionStep()
    {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        const double abs_tol = mMaxGlobalAbsoluteError;
		const double rel_tol = mMaxGlobalRelativeError;
        const double local_rel_tol = mMaxLocalRelativeError;

        // ---- Reference external force for relative convergence ----
        mIniExternalForceVector.resize(rb.size(), false);
        TSparseSpace::SetToZero(mIniExternalForceVector);
        BuildExternalForceVector(r_model_part, mIniExternalForceVector);
        mIniFNorm = TSparseSpace::TwoNorm(mIniExternalForceVector);

        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

        // ---- Build K0 ONCE and LU-factor it ONCE ----
        //TSparseSpace::SetToZero(rA);
        TSparseSpace::SetToZero(rDx);
        TSparseSpace::SetToZero(rb);


        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(rA);
            p_builder_and_solver->BuildLHS(p_scheme, r_model_part, rA);
        }

        if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), local_rel_tol)) {
            return false;
        }
        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

        // Initial residual rb (= -R in Python notation: K * du = rb)
        this->BuildRHS(r_model_part, r_dof_set, rb);
        p_builder_and_solver->ApplyDirichletConditions(p_scheme, r_model_part, rA, rDx, rb);

        // One-time LU decomposition of K0 
        auto p_linear_solver = p_builder_and_solver->GetLinearSystemSolver();

        if (BaseType::mRebuildLevel > 1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            p_linear_solver->InitializeSolutionStep(rA, rDx, rb);
            BaseType::mStiffnessMatrixIsBuilt = true;
        }

        

        // ---- Broyden low-rank storage: K_k = K0 + sum_i U_i V_i^T ----
        std::vector<TSystemVectorType> U_list;   // numerators
        std::vector<TSystemVectorType> V_list;   // denominators
        std::vector<TSystemVectorType> Z_list;   // K0^{-1} U_i (cached)

        const std::size_t n = rb.size();
        TSystemVectorType w(n);
        TSystemVectorType rb_new(n);
        TSystemVectorType K0s(n);
        TSystemVectorType u_col(n);
        TSystemVectorType z_col(n);
        TSystemVectorType rb_solve(n);

        bool is_converged = false;

        for (; iteration_number <= mMaxIterationNumber; ++iteration_number) {

            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
            if (iteration_number > 1) {
                p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
                mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
            }
            r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
            r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

            // ---- Convergence test on current residual rb ----
            const double rb_norm = TSparseSpace::TwoNorm(rb);
            is_converged = (rb_norm < abs_tol) ||
                           (mIniFNorm > 0.0 && (rb_norm / mIniFNorm) < rel_tol);
            std::cout << "[Broyden] iter " << iteration_number
                      << "  ||rb||=" << rb_norm
                      << "  ||F0||=" << mIniFNorm
                      << "  k="      << U_list.size()
                      << "  converged=" << is_converged << std::endl;
            if (is_converged) break;

            // ---- Woodbury solve:  du = (K0 + U V^T)^{-1} rb ----
            // w = K0^{-1} rb  (sparse triangular solve; reuses cached LU factorization)
            TSparseSpace::SetToZero(w);
            TSparseSpace::Copy(rb, rb_solve); // PerformSolutionStep may consume the RHS
            p_linear_solver->PerformSolutionStep(rA, w, rb_solve);

            TSystemVectorType du(n);
            const std::size_t k = U_list.size();
            if (k > 0) {
                // V^T w  (k-vector)
                Vector VtW(k);
                for (std::size_t i = 0; i < k; ++i) {
                    VtW[i] = TSparseSpace::Dot(V_list[i], w);
                }
                // M = I + V^T Z   (k x k)
                Matrix M(k, k);
                for (std::size_t i = 0; i < k; ++i) {
                    for (std::size_t j = 0; j < k; ++j) {
                        M(i, j) = (i == j ? 1.0 : 0.0) + TSparseSpace::Dot(V_list[i], Z_list[j]);
                    }
                }
                Vector alpha;
                if (!SolveSmallDense(M, VtW, alpha)) {
                    KRATOS_WARNING("ResidualBasedNewtonRaphsonStrategyTwo[Broyden]")
                        << "Woodbury capacitance matrix singular at iter " << iteration_number
                        << " (k=" << k << "); falling back to pure K0 step." << std::endl;
                    noalias(du) = w;
                } else {
                    // du = w - Z * alpha
                    noalias(du) = w;
                    for (std::size_t j = 0; j < k; ++j) {
                        noalias(du) -= alpha[j] * Z_list[j];
                    }
                }
            } else {
                noalias(du) = w;
            }

            // ---- Apply increment to the DOF database ----
            TSparseSpace::Copy(du, rDx);

            double relaxation_factor = 1.5;
			rDx /= relaxation_factor; // under-relaxation to improve robustness of Broyden updates
            UpdateDatabase(rA, rDx, rb, BaseType::MoveMeshFlag());

            // ---- Compute the new residual at the updated state ----
            TSparseSpace::SetToZero(rb_new);
            this->BuildRHS(r_model_part, r_dof_set, rb_new);
            if (!LocalConvergenceAchieved(r_model_part.GetProcessInfo(), local_rel_tol)) {
                return false;
            }

            // ---- Build the new Broyden rank-1 column ----
            //   s     = du
            //   K_k s = K0 s + U (V^T s)
            //   y     = rb - rb_new  
            //   u_col = y - K_k s
            //   v_col = s / (s . s)
            //   z_col = K0^{-1} u_col   (reuses factorization)
            const double s_dot_s = TSparseSpace::Dot(du, du);
            if (s_dot_s > std::numeric_limits<double>::epsilon()) {
                TSparseSpace::Mult(rA, du, K0s); // K0 * s
                if (k > 0) {
                    Vector VtS(k);
                    for (std::size_t i = 0; i < k; ++i) {
                        VtS[i] = TSparseSpace::Dot(V_list[i], du);
                    }
                    for (std::size_t i = 0; i < k; ++i) {
                        noalias(K0s) += VtS[i] * U_list[i];
                    }
                }

                TSystemVectorType y_col(n);
                noalias(y_col) = rb - rb_new;

                noalias(u_col) = y_col - K0s;

                TSystemVectorType v_col(n);
                noalias(v_col) = du / s_dot_s;

                TSparseSpace::SetToZero(z_col);
                TSystemVectorType u_solve(n);
                TSparseSpace::Copy(u_col, u_solve);
                p_linear_solver->PerformSolutionStep(rA, z_col, u_solve);

                U_list.push_back(u_col);
                V_list.push_back(v_col);
                Z_list.push_back(z_col);
            }

            // Advance: R := R_new   <=>   rb := rb_new
            TSparseSpace::Copy(rb_new, rb);

            // Per-iteration housekeeping
            EchoInfo(iteration_number);
            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        }

        if (!is_converged) {
            MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo[Broyden]", this->GetEchoLevel() > 0)
                << "Convergence achieved after " << iteration_number
                << " / " << mMaxIterationNumber << " Broyden iterations (rank-"
                << U_list.size() << " update)" << std::endl;
        }

        if (mCalculateReactionsFlag) {
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);
        }

        return is_converged;

        KRATOS_CATCH("")
    }


    ///@}
    ///@name Private Operators
    ///@{

    /**
     * @brief Here the database is updated
     * @param A The LHS matrix of the system of equations
     * @param Dx The incremement in the solution
     * @param b The RHS vector of the system of equations
     * @param MoveMesh The flag that allows to move the mesh
     */
    virtual void UpdateDatabase(
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb,
        const bool MoveMesh)
    {
        typename TSchemeType::Pointer p_scheme = GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = GetBuilderAndSolver();

        p_scheme->Update(BaseType::GetModelPart(), p_builder_and_solver->GetDofSet(), rA, rDx, rb);

        // Move the mesh if needed
        if (MoveMesh == true)
            BaseType::MoveMesh();
    }

    /**
     * @brief This method returns the components of the system of equations depending of the echo level
     * @param IterationNumber The non linear iteration in the solution loop
     */
    virtual void EchoInfo(const unsigned int IterationNumber)
    {
        TSystemMatrixType& rA  = *mpA;
        TSystemVectorType& rDx = *mpDx;
        TSystemVectorType& rb  = *mpb;

        if (this->GetEchoLevel() == 2) //if it is needed to print the debug info
        {
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 3) //if it is needed to print the debug info
        {
            KRATOS_INFO("LHS") << "SystemMatrix = " << rA << std::endl;
            KRATOS_INFO("Dx")  << "Solution obtained = " << rDx << std::endl;
            KRATOS_INFO("RHS") << "RHS  = " << rb << std::endl;
        }
        else if (this->GetEchoLevel() == 4) //print to matrix market file
        {
            std::stringstream matrix_market_name;
            matrix_market_name << "A_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm";
            TSparseSpace::WriteMatrixMarketMatrix((char *)(matrix_market_name.str()).c_str(), rA, false);

            std::stringstream matrix_market_vectname;
            matrix_market_vectname << "b_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_vectname.str()).c_str(), rb);

            std::stringstream matrix_market_dxname;
            matrix_market_dxname << "dx_" << BaseType::GetModelPart().GetProcessInfo()[TIME] << "_" << IterationNumber << ".mm.rhs";
            TSparseSpace::WriteMatrixMarketVector((char *)(matrix_market_dxname.str()).c_str(), rDx);

            std::stringstream dof_data_name;
            unsigned int rank=BaseType::GetModelPart().GetCommunicator().MyPID();
            dof_data_name << "dofdata_" << BaseType::GetModelPart().GetProcessInfo()[TIME]
                << "_" << IterationNumber << "_rank_"<< rank << ".csv";
            WriteDofInfo(dof_data_name.str(), rDx);
        }
    }

    /**
     * @brief This method prints information after reach the max number of iterations
     */
    virtual void MaxIterationsExceeded()
    {
        KRATOS_INFO_IF("ResidualBasedNewtonRaphsonStrategyTwo", this->GetEchoLevel() > 0)
            << "ATTENTION: max iterations ( " << mMaxIterationNumber
            << " ) exceeded!" << std::endl;
    }

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
        mMaxIterationNumber = ThisParameters["max_iteration"].GetInt();
        mReformDofSetAtEachStep = ThisParameters["reform_dofs_at_each_step"].GetBool();
        mCalculateReactionsFlag = ThisParameters["compute_reactions"].GetBool();
        mUseOldStiffnessInFirstIteration = ThisParameters["use_old_stiffness_in_first_iteration"].GetBool();

        // Arc-length settings
        std::cout << "use_arc_length: " << ThisParameters["use_arc_length"].GetBool();
        if (ThisParameters.Has("use_arc_length")) {
            mUseArcLength = ThisParameters["use_arc_length"].GetBool();
        }
        if (ThisParameters.Has("use_broyden")) {
            mUseBroyden = ThisParameters["use_broyden"].GetBool();
        }
        if (ThisParameters.Has("arc_length_settings")) {
            const Parameters al = ThisParameters["arc_length_settings"];
            if (al.Has("initial_arc_length"))        mArcLength0            = al["initial_arc_length"].GetDouble();
            if (al.Has("max_arc_length_factor"))     mMaxArcLengthFactor    = al["max_arc_length_factor"].GetDouble();
            if (al.Has("min_arc_length_factor"))     mMinArcLengthFactor    = al["min_arc_length_factor"].GetDouble();
            if (al.Has("desired_iterations"))        mDesiredIterations     = al["desired_iterations"].GetInt();
            if (al.Has("beta"))                      mBeta                  = al["beta"].GetDouble();
            if (al.Has("max_load_factor"))           mLambdaMax             = al["max_load_factor"].GetDouble();
            if (al.Has("max_arc_length_reductions")) mMaxArcLengthReductions= al["max_arc_length_reductions"].GetInt();
            mArcLength = mArcLength0;
        }

        // Saving the convergence criteria to be used
        if (ThisParameters["convergence_criteria_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

        // Saving the scheme
        if (ThisParameters["scheme_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }

        // Setting up the default builder and solver
        if (ThisParameters["builder_and_solver_settings"].Has("name")) {
            KRATOS_ERROR << "IMPLEMENTATION PENDING IN CONSTRUCTOR WITH PARAMETERS" << std::endl;
        }
    }

    void WriteDofInfo(std::string FileName, const TSystemVectorType& rDX)
    {
        std::ofstream out(FileName);

        out.precision(15);
        out << "EquationId,NodeId,VariableName,IsFixed,Value,coordx,coordy,coordz" << std::endl;
        for(const auto& rdof : GetBuilderAndSolver()->GetDofSet()) {
            const auto& coords = BaseType::GetModelPart().Nodes()[rdof.Id()].Coordinates();
            out << rdof.EquationId() << "," << rdof.Id() << "," << rdof.GetVariable().Name() << "," << rdof.IsFixed() << ","
                        << rdof.GetSolutionStepValue() << "," <<  "," << coords[0]  << "," << coords[1]  << "," << coords[2]<< "\n";
        }
        out.close();
    }

    void save(FileSerializer& rSerializer) const
    {
        rSerializer.save("PreviousExternalForceVector", mIniExternalForceVector);
        rSerializer.save("ArcLengthLambda", mLambda);
        rSerializer.save("ArcLength", mArcLength);
        rSerializer.save("ArcLength0", mArcLength0);
        rSerializer.save("ArcLengthSignHint", mPredictorSignHint);
    }

    void load(FileSerializer& rSerializer)
    {
        rSerializer.load("PreviousExternalForceVector", mPreviousExternalForceVector);
        rSerializer.load("ArcLengthLambda", mLambda);
        rSerializer.load("ArcLength", mArcLength);
        rSerializer.load("ArcLength0", mArcLength0);
        rSerializer.load("ArcLengthSignHint", mPredictorSignHint);
    }

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /**
     * Copy constructor.
     */

    ResidualBasedNewtonRaphsonStrategyTwo(const ResidualBasedNewtonRaphsonStrategyTwo&Other){};

    ///@}

}; /* Class ResidualBasedNewtonRaphsonStrategyTwo */

///@}

///@name Type Definitions
///@{

///@}

} /* namespace Kratos. */



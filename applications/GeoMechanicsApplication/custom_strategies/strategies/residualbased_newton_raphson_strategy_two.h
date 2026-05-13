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
            double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
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
            double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
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
        double load_fraction = (r_model_part.GetProcessInfo()[TIME] - r_model_part.GetProcessInfo()[START_TIME]) / (r_model_part.GetProcessInfo()[END_TIME] - r_model_part.GetProcessInfo()[START_TIME]);
        
        TSystemVectorType d_ext_force = mIniExternalForceVector - mPreviousExternalForceVector;
        
        TSystemVectorType int_force = ZeroVector(rb.size());

        // reset plastic points
        r_model_part.GetProcessInfo()[INACCURATE_PLASTIC_POINTS] = 0;
        r_model_part.GetProcessInfo()[N_PLASTIC_POINTS] = 0;

        BuildInternalForceVector(r_model_part, int_force);
        
        TSystemVectorType ini_int_force = ZeroVector(rb.size());

		// initial internal force is equal to external force from previous step, as the system is in equilibrium at the beginning of the step
        TSparseSpace::Copy(mPreviousExternalForceVector, ini_int_force);
        TSystemVectorType d_int_force = int_force - ini_int_force;

        // norm int force
        //double int_force_norm = TSparseSpace::TwoNorm(int_force);
        //std::cout << "int_force_norm: " << int_force_norm << std::endl;

		// load fraction * (d_Fext - d_Fint) + Fext_prev - Fint_prev
        rb = load_fraction * (d_ext_force - d_int_force);

		// apply Dirichlet conditions RHS
        //NOTE: dofs are assumed to be numbered consecutively
        block_for_each(r_dof_set, [&](Dof<double>& rDof) {
            const std::size_t i = rDof.EquationId();

            if (rDof.IsFixed())
                rb[i] = 0.0;
            });
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
    }

    void load(FileSerializer& rSerializer)
    {
        rSerializer.load("PreviousExternalForceVector", mPreviousExternalForceVector);
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
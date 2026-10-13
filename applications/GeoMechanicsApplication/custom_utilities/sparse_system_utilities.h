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

#include <algorithm>
#include <vector>

// External includes
#include "containers/system_vector.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "spaces/ublas_space.h"

namespace Kratos::Geo
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) SparseSystemUtilities
{
public:
    /// Vector type definition
    using SystemVectorType = UblasSpace<double, CompressedMatrix, Vector>::VectorType;
    using SystemMatrixType = UblasSpace<double, CompressedMatrix, Vector>::MatrixType;
    /// DoF array type definition
    using DofsArrayType = ModelPart::DofsArrayType;

    static void GetTotalSolutionStepValueVector(SystemVectorType&    rTotalSolutionStepValues,
                                                const DofsArrayType& rDofSet,
                                                const ModelPart&     rModelPart,
                                                const IndexType      BufferIndex);

    static void GetUFirstAndSecondDerivativeVector(SystemVectorType&    rFirstDerivativeVector,
                                                   SystemVectorType&    rSecondDerivativeVector,
                                                   const DofsArrayType& rDofSet,
                                                   const ModelPart&     rModelPart,
                                                   const IndexType      BufferIndex);

    static void SetUFirstAndSecondDerivativeVector(const SystemVectorType& rFirstDerivativeVector,
                                                   const SystemVectorType& rSecondDerivativeVector,
                                                   const DofsArrayType&    rDofSet,
                                                   ModelPart&              rModelPart);

    /// <summary>
    /// In this function, the rows and columns of the secondary matrix corresponding to the fixed DOFs are zeroed,
    /// the diagonal is kept as it is, as for secondary matrices, zeros are allowed on the diagonal.
    /// </summary>
    /// <param name="rDofSet"></param>
    /// <param name="rSecondaryMatrix"></param>
    static void ApplyDirichletConditionsSecondaryMatrix(const DofsArrayType& rDofSet,
                                                        SystemMatrixType&    rSecondaryMatrix);

    /// <summary>
    /// Checks if the two matrices have the same signature on the diagonal, meaning that they have non-zeros in the same positions on the diagonal.
    /// </summary>
    /// <param name="rPrimaryMatrix"></param>
    /// <param name="rSecondaryMatrix"></param>
    /// <param name="rDofSet"></param>
    /// <returns></returns>
    static bool MatricesHaveSameDiagonalSignature(const SystemMatrixType& rPrimaryMatrix,
                                                  const SystemMatrixType& rSecondaryMatrix,
                                                  const DofsArrayType&    rDofSet);

    /// Returns true when the dof is part of the (consecutively numbered) dof set, i.e. when it is
    /// part of the system of equations. Dofs that have only been added to a node (e.g. the
    /// translations of an embedded beam node or the rotations of a soil node) are not.
    static bool IsPartOfDofSet(const Dof<double>& rDof, const DofsArrayType& rDofSet);

private:
    static void GetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  const Node&             rNode,
                                                  SystemVectorType&       rFirstDerivativeVector,
                                                  SystemVectorType&       rSecondDerivativeVector,
                                                  const DofsArrayType&    rDofSet,
                                                  const IndexType         BufferIndex);

    static void SetDerivativesForOptionalVariable(const Variable<double>& rVariable,
                                                  Node&                   rNode,
                                                  const SystemVectorType& rFirstDerivativeVector,
                                                  const SystemVectorType& rSecondDerivativeVector,
                                                  const DofsArrayType&    rDofSet);

    static void GetDerivativesForVariable(const Variable<double>& rVariable,
                                          const Node&             rNode,
                                          SystemVectorType&       rFirstDerivativeVector,
                                          SystemVectorType&       rSecondDerivativeVector,
                                          const DofsArrayType&    rDofSet,
                                          const IndexType         BufferIndex);

    static void SetDerivativesForVariable(const Variable<double>& rVariable,
                                          Node&                   rNode,
                                          const SystemVectorType& rFirstDerivativeVector,
                                          const SystemVectorType& rSecondDerivativeVector,
                                          const DofsArrayType&    rDofSet);

    /// <summary>
    /// Checks if current row has a non-zero entry on the diagonal,
    /// </summary>
    /// <param name="RowIndex"></param>
    /// <param name="rCsrIndices1"></param>
    /// <param name="rCsrIndices2"></param>
    /// <param name="rCsrValues"></param>
    /// <returns></returns>
    static bool HasNonZeroDiagonalEntryOnCurrentRow(const std::size_t RowIndex,
                                                    const unbounded_array<std::size_t>& rCsrIndices1,
                                                    const unbounded_array<std::size_t>& rCsrIndices2,
                                                    const unbounded_array<double>& rCsrValues);
};
} // namespace Kratos::Geo

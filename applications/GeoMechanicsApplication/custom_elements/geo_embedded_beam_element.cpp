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

#include "custom_elements/geo_embedded_beam_element.h"
#include "custom_utilities/dof_utilities.hpp"
#include "custom_utilities/equation_of_motion_utilities.hpp"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/variables.h"
#include "structural_mechanics_application_variables.h"

#include <array>

namespace
{

using namespace Kratos;

std::vector<const Variable<double>*> GetRotationVariables(unsigned int Dimension)
{
    return Dimension == 2 ? std::vector<const Variable<double>*>{&ROTATION_Z}
                          : std::vector<const Variable<double>*>{&ROTATION_X, &ROTATION_Y, &ROTATION_Z};
}

const std::array<const Variable<double>*, 3>& GetDisplacementVariables()
{
    static const std::array<const Variable<double>*, 3> result = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
    return result;
}

} // namespace

namespace Kratos
{

template <unsigned int TDim>
GeoEmbeddedBeamElement<TDim>::GeoEmbeddedBeamElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
{
}

template <unsigned int TDim>
GeoEmbeddedBeamElement<TDim>::GeoEmbeddedBeamElement(IndexType                        NewId,
                                                     typename GeometryType::Pointer   pGeometry,
                                                     typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim>
Element::Pointer GeoEmbeddedBeamElement<TDim>::Create(IndexType             NewId,
                                                      const NodesArrayType& rNodes,
                                                      typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoEmbeddedBeamElement>(NewId, this->GetGeometry().Create(rNodes), pProperties);
}

template <unsigned int TDim>
Element::Pointer GeoEmbeddedBeamElement<TDim>::Create(IndexType                      NewId,
                                                      typename GeometryType::Pointer pGeometry,
                                                      typename PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoEmbeddedBeamElement>(NewId, pGeometry, pProperties);
}

template <unsigned int TDim>
Element::Pointer GeoEmbeddedBeamElement<TDim>::Clone(IndexType NewId, const NodesArrayType& rNodes) const
{
    KRATOS_TRY

    auto p_new_element = Kratos::make_intrusive<GeoEmbeddedBeamElement>(
        NewId, this->GetGeometry().Create(rNodes), this->pGetProperties());
    p_new_element->SetData(this->GetData());
    p_new_element->Set(Flags(*this));
    p_new_element->SetIntegrationMethod(this->mThisIntegrationMethod);
    p_new_element->SetConstitutiveLawVector(this->mConstitutiveLawVector);
    p_new_element->mHostInterpolation = mHostInterpolation;
    return p_new_element;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    BaseType::Initialize(rCurrentProcessInfo);

    if (this->Has(NEIGHBOUR_ELEMENTS)) {
        mHostInterpolation.Initialize(this->GetGeometry(), this->GetValue(NEIGHBOUR_ELEMENTS),
                                      msLocalCoordinatesTolerance);
    }

    KRATOS_ERROR_IF_NOT(mHostInterpolation.IsInitialized())
        << "No candidate host elements have been assigned to embedded beam element " << this->Id()
        << ". Please assign them to NEIGHBOUR_ELEMENTS, e.g. by means of the "
           "FindEmbeddedBeamHostElementsProcess."
        << std::endl;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::FinalizeNonLinearIteration(rCurrentProcessInfo);
    UpdateEmbeddedBeamNodalKinematics();
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    BaseType::FinalizeSolutionStep(rCurrentProcessInfo);
    UpdateEmbeddedBeamNodalKinematics();
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::UpdateEmbeddedBeamNodalKinematics()
{
    // The translations of embedded beam nodes are not degrees of freedom of the system. For
    // postprocessing purposes, they are interpolated from the host nodes.
    const std::array<const Variable<array_1d<double, 3>>*, 3> variables = {&DISPLACEMENT, &VELOCITY, &ACCELERATION};

    auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < msNumberOfBeamNodes; ++i) {
        auto& r_node = r_geometry[i];
        if (mHostInterpolation.IsInterpolatedFromItself(r_node, i)) continue;

        for (const auto* p_variable : variables) {
            if (!r_node.SolutionStepsDataHas(*p_variable)) continue;

            const auto value = mHostInterpolation.Interpolate(*p_variable, i);
            r_node.SetLock();
            r_node.FastGetSolutionStepValue(*p_variable) = value;
            r_node.UnSetLock();
        }
    }
}

template <unsigned int TDim>
std::vector<Dof<double>*> GeoEmbeddedBeamElement<TDim>::GetDofs() const
{
    const auto& r_interpolation_nodes = mHostInterpolation.GetInterpolationNodes();
    KRATOS_ERROR_IF(r_interpolation_nodes.empty())
        << "Embedded beam element " << this->Id() << " has not been initialized" << std::endl;

    std::vector<Dof<double>*> result;
    result.reserve(r_interpolation_nodes.size() * TDim + msNumberOfBeamNodes * msNumberOfRotations);

    const auto& r_displacement_variables = GetDisplacementVariables();
    for (const auto& rp_node : r_interpolation_nodes) {
        for (IndexType d = 0; d < TDim; ++d) {
            result.push_back(rp_node->pGetDof(*r_displacement_variables[d]));
        }
    }

    const auto rotation_variables = GetRotationVariables(TDim);
    for (const auto& r_node : this->GetGeometry()) {
        for (const auto* p_variable : rotation_variables) {
            result.push_back(r_node.pGetDof(*p_variable));
        }
    }

    return result;
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::GetValuesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSolutionStepValues(GetDofs(), Step);
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivatives(GetDofs(), Step);
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivatives(GetDofs(), Step);
}

template <unsigned int TDim>
Matrix GeoEmbeddedBeamElement<TDim>::CalculateTransformationMatrix() const
{
    const auto number_of_interpolation_nodes = mHostInterpolation.GetInterpolationNodes().size();
    KRATOS_ERROR_IF(number_of_interpolation_nodes == 0)
        << "Embedded beam element " << this->Id() << " has not been initialized" << std::endl;

    const auto index_of_first_rotation = number_of_interpolation_nodes * TDim;
    Matrix     result                  = ZeroMatrix(msNumberOfBeamDofs,
                                                    index_of_first_rotation + msNumberOfBeamNodes * msNumberOfRotations);

    for (IndexType i = 0; i < msNumberOfBeamNodes; ++i) {
        const auto  first_beam_dof = i * msDofsPerBeamNode;
        const auto& r_indices      = mHostInterpolation.GetInterpolationNodeIndices(i);
        const auto& r_weights      = mHostInterpolation.GetInterpolationWeights(i);

        // U^b_i = N^s(X_i) U^s
        for (IndexType k = 0; k < r_indices.size(); ++k) {
            for (IndexType d = 0; d < TDim; ++d) {
                result(first_beam_dof + d, r_indices[k] * TDim + d) += r_weights[k];
            }
        }

        // Theta^b_i = Theta_i
        for (IndexType r = 0; r < msNumberOfRotations; ++r) {
            result(first_beam_dof + TDim + r, index_of_first_rotation + i * msNumberOfRotations + r) = 1.0;
        }
    }

    return result;
}

template <unsigned int TDim>
const EmbeddedHostInterpolation& GeoEmbeddedBeamElement<TDim>::GetHostInterpolation() const
{
    return mHostInterpolation;
}

template <unsigned int TDim>
auto GeoEmbeddedBeamElement<TDim>::CalculateGlobalToLocalRotationMatrix() const
    -> BoundedMatrix<double, msNumberOfBeamDofs, msNumberOfBeamDofs>
{
    // The local degrees of freedom of a beam node are [u, v, theta] in 2D and
    // [u, v, w, theta_x, -theta_y, theta_z] in 3D, where the sign of the rotation about the local
    // y axis is inverted, such that both bending planes follow the same convention (theta = dv/dx)
    const auto local_axes = EmbeddedBeamUtilities::CalculateLocalAxes(this->GetGeometry(), TDim);

    BoundedMatrix<double, msNumberOfBeamDofs, msNumberOfBeamDofs> result =
        ZeroMatrix(msNumberOfBeamDofs, msNumberOfBeamDofs);
    for (IndexType node = 0; node < msNumberOfBeamNodes; ++node) {
        const auto offset = node * msDofsPerBeamNode;
        for (IndexType i = 0; i < TDim; ++i) {
            for (IndexType j = 0; j < TDim; ++j) {
                result(offset + i, offset + j) = local_axes(i, j);
            }
        }

        if constexpr (TDim == 2) {
            result(offset + 2, offset + 2) = 1.0;
        } else {
            for (IndexType i = 0; i < 3; ++i) {
                const auto sign = (i == 1) ? -1.0 : 1.0;
                for (IndexType j = 0; j < 3; ++j) {
                    result(offset + 3 + i, offset + 3 + j) = sign * local_axes(i, j);
                }
            }
        }
    }

    return result;
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::GetNodalValuesVector(VectorType& rNodalValues) const
{
    Vector dof_values;
    this->GetValuesVector(dof_values, 0);
    const Vector global_beam_values = prod(CalculateTransformationMatrix(), dof_values);

    if (rNodalValues.size() != msNumberOfBeamDofs) rNodalValues.resize(msNumberOfBeamDofs, false);
    noalias(rNodalValues) = prod(CalculateGlobalToLocalRotationMatrix(), global_beam_values);
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::RotateLHS(MatrixType& rLHS, const GeometryType&)
{
    const auto   rotation_matrix = CalculateGlobalToLocalRotationMatrix();
    const Matrix auxiliary       = prod(rLHS, rotation_matrix);
    noalias(rLHS)                = prod(trans(rotation_matrix), auxiliary);
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::RotateRHS(VectorType& rRHS, const GeometryType&)
{
    const auto   rotation_matrix = CalculateGlobalToLocalRotationMatrix();
    const Vector local_rhs       = rRHS;
    noalias(rRHS)                = prod(trans(rotation_matrix), local_rhs);
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::RotateAll(MatrixType& rLHS, VectorType& rRHS, const GeometryType& rGeometry)
{
    RotateLHS(rLHS, rGeometry);
    RotateRHS(rRHS, rGeometry);
}

template <unsigned int TDim>
array_1d<double, 3> GeoEmbeddedBeamElement<TDim>::GetLocalAxesBodyForce(
    const Element&, const typename GeometryType::IntegrationPointsArrayType&, IndexType) const
{
    return ZeroVector(3);
}

template <unsigned int TDim>
double GeoEmbeddedBeamElement<TDim>::GetBeamCrossArea()
{
    if constexpr (TDim == 2) {
        return this->GetCrossArea();
    } else {
        return this->GetProperties()[CROSS_AREA];
    }
}

template <unsigned int TDim>
Vector GeoEmbeddedBeamElement<TDim>::CalculateBodyForceVector()
{
    KRATOS_TRY

    // Consistent body forces in local axes, using the shape functions of the Timoshenko beam
    const auto& r_properties         = this->GetProperties();
    const auto& r_geometry           = this->GetGeometry();
    const auto  length               = EmbeddedBeamUtilities::CalculateReferenceLength(r_geometry);
    const auto  area                 = GetBeamCrossArea();
    const auto  local_axes           = EmbeddedBeamUtilities::CalculateLocalAxes(r_geometry, TDim);
    const auto& r_integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    VectorType local_body_forces = ZeroVector(msNumberOfBeamDofs);
    VectorType global_size_N(msNumberOfBeamDofs);
    VectorType N_u(msNumberOfBeamNodes);
    VectorType N_shape(4);
    for (IndexType integration_point = 0; integration_point < r_integration_points.size(); ++integration_point) {
        const array_1d<double, 3> local_body_force =
            prod(local_axes, StructuralMechanicsElementUtilities::GetBodyForce(
                                 *this, r_integration_points, integration_point));
        const auto xi = r_integration_points[integration_point].X();
        const auto integration_coefficient =
            r_integration_points[integration_point].Weight() * 0.5 * length * area;

        this->GetNu0ShapeFunctionsValues(N_u, length, 0.0, xi);
        this->GlobalSizeAxialVector(global_size_N, N_u);
        noalias(local_body_forces) += global_size_N * local_body_force[0] * integration_coefficient;

        if constexpr (TDim == 2) {
            this->GetShapeFunctionsValues(
                N_shape, length, StructuralMechanicsElementUtilities::CalculatePhi(r_properties, length), xi);
            this->GlobalSizeVector(global_size_N, N_shape);
            noalias(local_body_forces) += global_size_N * local_body_force[1] * integration_coefficient;
        } else {
            this->GetShapeFunctionsValues(
                N_shape, length, StructuralMechanicsElementUtilities::CalculatePhi(r_properties, length), xi);
            this->GlobalSizeVectorTransversalY(global_size_N, N_shape);
            noalias(local_body_forces) += global_size_N * local_body_force[1] * integration_coefficient;

            this->GetShapeFunctionsValues(
                N_shape, length,
                StructuralMechanicsElementUtilities::CalculatePhi(r_properties, length, 1), xi);
            this->GlobalSizeVectorTransversalZ(global_size_N, N_shape);
            noalias(local_body_forces) += global_size_N * local_body_force[2] * integration_coefficient;
        }
    }

    return prod(trans(CalculateGlobalToLocalRotationMatrix()), local_body_forces);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
Matrix GeoEmbeddedBeamElement<TDim>::CalculateLumpedBeamMassMatrix()
{
    KRATOS_TRY

    const auto& r_properties = this->GetProperties();
    const auto  length       = EmbeddedBeamUtilities::CalculateReferenceLength(this->GetGeometry());
    const auto  density      = r_properties.Has(DENSITY) ? r_properties[DENSITY] : 0.0;
    const auto  half_mass_per_unit_inertia = 0.5 * density * length;
    const auto  translational_mass         = half_mass_per_unit_inertia * GetBeamCrossArea();

    const auto second_moment_of_area_y = r_properties.Has(I22) ? r_properties[I22] : 0.0;
    const auto second_moment_of_area_z = r_properties.Has(I33) ? r_properties[I33] : 0.0;

    // Lumped rotary inertia about the local axes
    BoundedMatrix<double, 3, 3> rotary_inertia = ZeroMatrix(3, 3);
    if constexpr (TDim == 2) {
        rotary_inertia(2, 2) = half_mass_per_unit_inertia * second_moment_of_area_z;
    } else {
        const auto local_axes = EmbeddedBeamUtilities::CalculateLocalAxes(this->GetGeometry(), TDim);
        BoundedMatrix<double, 3, 3> local_rotary_inertia = ZeroMatrix(3, 3);
        local_rotary_inertia(0, 0) =
            half_mass_per_unit_inertia * (second_moment_of_area_y + second_moment_of_area_z);
        local_rotary_inertia(1, 1) = half_mass_per_unit_inertia * second_moment_of_area_y;
        local_rotary_inertia(2, 2) = half_mass_per_unit_inertia * second_moment_of_area_z;
        const BoundedMatrix<double, 3, 3> auxiliary = prod(local_rotary_inertia, local_axes);
        noalias(rotary_inertia)                     = prod(trans(local_axes), auxiliary);
    }

    Matrix result = ZeroMatrix(msNumberOfBeamDofs, msNumberOfBeamDofs);
    for (IndexType node = 0; node < msNumberOfBeamNodes; ++node) {
        const auto offset = node * msDofsPerBeamNode;
        for (IndexType d = 0; d < TDim; ++d) {
            result(offset + d, offset + d) = translational_mass;
        }

        if constexpr (TDim == 2) {
            result(offset + 2, offset + 2) = rotary_inertia(2, 2);
        } else {
            for (IndexType i = 0; i < 3; ++i) {
                for (IndexType j = 0; j < 3; ++j) {
                    result(offset + 3 + i, offset + 3 + j) = rotary_inertia(i, j);
                }
            }
        }
    }

    return result;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                                                        VectorType&        rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MatrixType beam_left_hand_side;
    VectorType beam_right_hand_side;
    BaseType::CalculateLocalSystem(beam_left_hand_side, beam_right_hand_side, rCurrentProcessInfo);
    noalias(beam_right_hand_side) += CalculateBodyForceVector();

    // Eq. (11): R^bs = N^bs^T R^b N^bs
    const auto   transformation_matrix = CalculateTransformationMatrix();
    const Matrix auxiliary             = prod(beam_left_hand_side, transformation_matrix);
    rLeftHandSideMatrix                = prod(trans(transformation_matrix), auxiliary);
    rRightHandSideVector               = prod(trans(transformation_matrix), beam_right_hand_side);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::CalculateLeftHandSide(MatrixType&        rLeftHandSideMatrix,
                                                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MatrixType beam_left_hand_side;
    BaseType::CalculateLeftHandSide(beam_left_hand_side, rCurrentProcessInfo);

    const auto   transformation_matrix = CalculateTransformationMatrix();
    const Matrix auxiliary             = prod(beam_left_hand_side, transformation_matrix);
    rLeftHandSideMatrix                = prod(trans(transformation_matrix), auxiliary);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::CalculateRightHandSide(VectorType&        rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    VectorType beam_right_hand_side;
    BaseType::CalculateRightHandSide(beam_right_hand_side, rCurrentProcessInfo);
    noalias(beam_right_hand_side) += CalculateBodyForceVector();

    rRightHandSideVector = prod(trans(CalculateTransformationMatrix()), beam_right_hand_side);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo&)
{
    KRATOS_TRY

    const auto   transformation_matrix = CalculateTransformationMatrix();
    const Matrix auxiliary = prod(CalculateLumpedBeamMassMatrix(), transformation_matrix);
    rMassMatrix            = prod(trans(transformation_matrix), auxiliary);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    MatrixType mass_matrix;
    CalculateMassMatrix(mass_matrix, rCurrentProcessInfo);
    MatrixType stiffness_matrix;
    CalculateLeftHandSide(stiffness_matrix, rCurrentProcessInfo);

    const auto& r_properties = this->GetProperties();
    rDampingMatrix           = GeoEquationOfMotionUtilities::CalculateDampingMatrix(
        r_properties.Has(RAYLEIGH_ALPHA) ? r_properties[RAYLEIGH_ALPHA] : rCurrentProcessInfo[RAYLEIGH_ALPHA],
        r_properties.Has(RAYLEIGH_BETA) ? r_properties[RAYLEIGH_BETA] : rCurrentProcessInfo[RAYLEIGH_BETA],
        mass_matrix, stiffness_matrix);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
int GeoEmbeddedBeamElement<TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    const auto result = BaseType::Check(rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(mHostInterpolation.IsInitialized())
        << "Embedded beam element " << this->Id() << " has not been initialized" << std::endl;

    const auto& r_displacement_variables = GetDisplacementVariables();
    for (const auto& rp_node : mHostInterpolation.GetInterpolationNodes()) {
        for (IndexType d = 0; d < TDim; ++d) {
            KRATOS_ERROR_IF_NOT(rp_node->HasDofFor(*r_displacement_variables[d]))
                << "Node " << rp_node->Id() << " of (host of) embedded beam element " << this->Id()
                << " does not have degree of freedom " << r_displacement_variables[d]->Name() << std::endl;
        }
    }

    for (const auto& r_node : this->GetGeometry()) {
        for (const auto* p_variable : GetRotationVariables(TDim)) {
            KRATOS_ERROR_IF_NOT(r_node.HasDofFor(*p_variable))
                << "Node " << r_node.Id() << " of embedded beam element " << this->Id()
                << " does not have degree of freedom " << p_variable->Name()
                << ". Please make sure the rotation degrees of freedom are added (e.g. by setting "
                   "'rotation_dofs' to true in the solver settings)."
                << std::endl;
        }
    }

    return result;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
std::string GeoEmbeddedBeamElement<TDim>::Info() const
{
    return "GeoEmbeddedBeamElement" + std::to_string(TDim) + "D2N #" + std::to_string(this->Id());
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
}

template <unsigned int TDim>
void GeoEmbeddedBeamElement<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
}

template class GeoEmbeddedBeamElement<2>;
template class GeoEmbeddedBeamElement<3>;

} // namespace Kratos

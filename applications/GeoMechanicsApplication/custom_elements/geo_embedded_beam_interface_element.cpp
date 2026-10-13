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

#include "custom_elements/geo_embedded_beam_interface_element.h"
#include "custom_utilities/dof_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace
{

using namespace Kratos;

const std::array<const Variable<double>*, 3>& GetDisplacementVariables()
{
    static const std::array<const Variable<double>*, 3> result = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};
    return result;
}

} // namespace

namespace Kratos
{

template <unsigned int TDim>
GeoEmbeddedBeamInterfaceElement<TDim>::GeoEmbeddedBeamInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

template <unsigned int TDim>
GeoEmbeddedBeamInterfaceElement<TDim>::GeoEmbeddedBeamInterfaceElement(IndexType NewId,
                                                                       GeometryType::Pointer pGeometry,
                                                                       PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim>
Element::Pointer GeoEmbeddedBeamInterfaceElement<TDim>::Create(IndexType             NewId,
                                                               const NodesArrayType& rNodes,
                                                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoEmbeddedBeamInterfaceElement>(NewId, GetGeometry().Create(rNodes), pProperties);
}

template <unsigned int TDim>
Element::Pointer GeoEmbeddedBeamInterfaceElement<TDim>::Create(IndexType             NewId,
                                                               GeometryType::Pointer pGeometry,
                                                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoEmbeddedBeamInterfaceElement>(NewId, pGeometry, pProperties);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(GetGeometry().PointsNumber() == msNumberOfBeamNodes)
        << "Embedded beam interface element " << Id() << " requires a geometry with "
        << msNumberOfBeamNodes << " nodes" << std::endl;

    if (Has(NEIGHBOUR_ELEMENTS)) {
        mHostInterpolation.Initialize(GetGeometry(), GetValue(NEIGHBOUR_ELEMENTS), msLocalCoordinatesTolerance);
    }
    KRATOS_ERROR_IF_NOT(mHostInterpolation.IsInitialized())
        << "No candidate host elements have been assigned to embedded beam interface element " << Id()
        << ". Please assign them to NEIGHBOUR_ELEMENTS, e.g. by means of the "
           "FindEmbeddedBeamHostElementsProcess."
        << std::endl;

    // The beam nodes come first, followed by the host nodes that are not beam nodes
    mNodes.clear();
    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        mNodes.push_back(GetGeometry()(i));
    }

    const auto& r_interpolation_nodes = mHostInterpolation.GetInterpolationNodes();
    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        mHostNodeIndices[i].clear();
        for (const auto interpolation_node_index : mHostInterpolation.GetInterpolationNodeIndices(i)) {
            const auto& rp_node = r_interpolation_nodes[interpolation_node_index];
            const auto it = std::find_if(mNodes.begin(), mNodes.end(), [&rp_node](const auto& rpOther) {
                return rpOther->Id() == rp_node->Id();
            });
            if (it == mNodes.end()) {
                mNodes.push_back(rp_node);
                mHostNodeIndices[i].push_back(mNodes.size() - 1);
            } else {
                mHostNodeIndices[i].push_back(static_cast<std::size_t>(std::distance(mNodes.begin(), it)));
            }
        }
    }

    // The element is initialized at the start of every stage. Since DISPLACEMENT is reset at the
    // start of a stage, the converged bonding forces of the previous stage become the initial
    // bonding forces and the plastic slip is measured from the start of this stage. An inactive
    // element does not transfer any forces, such that it starts unloaded when it is reactivated.
    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        if (mConvergedTraction[i].size() != TDim || !IsActive()) {
            mConvergedTraction[i] = ZeroVector(TDim);
        }
        mInitialTraction[i]      = mConvergedTraction[i];
        mConvergedPlasticSlip[i] = 0.0;
    }

    // The host elements may not have been initialized yet, hence their stresses are not used here
    UpdateMaximumShearForces(rCurrentProcessInfo, false);

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    UpdateMaximumShearForces(rCurrentProcessInfo, true);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::FinalizeSolutionStep(const ProcessInfo&)
{
    KRATOS_TRY

    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        if (!mHostInterpolation.IsEmbedded(i)) continue;

        Vector traction;
        Matrix tangent;
        double plastic_slip = 0.0;
        CalculateTractionAndTangent(i, CalculateLocalRelativeDisplacement(i), traction, tangent, plastic_slip);
        mConvergedPlasticSlip[i] = plastic_slip;
        mConvergedTraction[i]    = traction;
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
bool GeoEmbeddedBeamInterfaceElement<TDim>::HasShearStrength() const
{
    return GetProperties().Has(GEO_COHESION) || GetProperties().Has(GEO_FRICTION_ANGLE);
}

template <unsigned int TDim>
double GeoEmbeddedBeamInterfaceElement<TDim>::GetPerimeter() const
{
    // In 2D (per unit width), the beam is in contact with the soil on both of its faces
    if constexpr (TDim == 2) {
        return 2.0;
    } else {
        return Globals::Pi * GetProperties()[DIAMETER];
    }
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::UpdateMaximumShearForces(const ProcessInfo& rCurrentProcessInfo,
                                                                     bool UseHostStresses)
{
    KRATOS_TRY

    if (!HasShearStrength()) {
        mMaximumShearForce.fill(std::numeric_limits<double>::max());
        return;
    }

    const auto& r_properties = GetProperties();
    const auto  cohesion     = r_properties.Has(GEO_COHESION) ? r_properties[GEO_COHESION] : 0.0;
    const auto  tan_friction_angle =
        r_properties.Has(GEO_FRICTION_ANGLE)
             ? std::tan(MathUtils<>::DegreesToRadians(r_properties[GEO_FRICTION_ANGLE]))
             : 0.0;

    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        const auto compressive_normal_stress =
            (UseHostStresses && tan_friction_angle > 0.0)
                ? std::max(0.0, -CalculateNormalStressOnShaft(i, rCurrentProcessInfo))
                : 0.0;
        mMaximumShearForce[i] = GetPerimeter() * (cohesion + compressive_normal_stress * tan_friction_angle);
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
double GeoEmbeddedBeamInterfaceElement<TDim>::CalculateNormalStressOnShaft(std::size_t BeamNodeIndex,
                                                                           const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    auto* p_host_element = mHostInterpolation.GetHostElement(BeamNodeIndex);
    if (!p_host_element) return 0.0;

    std::vector<Vector> stress_vectors;
    p_host_element->CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, stress_vectors, rCurrentProcessInfo);
    if (stress_vectors.empty() || stress_vectors.front().size() == 0) return 0.0;

    Vector average_stress_vector = ZeroVector(stress_vectors.front().size());
    for (const auto& r_stress_vector : stress_vectors) {
        noalias(average_stress_vector) += r_stress_vector;
    }
    average_stress_vector /= static_cast<double>(stress_vectors.size());

    const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor(average_stress_vector);
    const auto   local_axes    = EmbeddedBeamUtilities::CalculateLocalAxes(GetGeometry(), TDim);

    // Normal stress on the plane perpendicular to local axis `AxisIndex`
    auto normal_stress = [&stress_tensor, &local_axes](std::size_t AxisIndex) {
        auto result = 0.0;
        for (std::size_t i = 0; i < stress_tensor.size1(); ++i) {
            for (std::size_t j = 0; j < stress_tensor.size2(); ++j) {
                result += local_axes(AxisIndex, i) * stress_tensor(i, j) * local_axes(AxisIndex, j);
            }
        }
        return result;
    };

    if constexpr (TDim == 2) {
        return normal_stress(1);
    } else {
        return 0.5 * (normal_stress(1) + normal_stress(2));
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
std::vector<Dof<double>*> GeoEmbeddedBeamInterfaceElement<TDim>::GetDofs() const
{
    KRATOS_ERROR_IF(mNodes.empty())
        << "Embedded beam interface element " << Id() << " has not been initialized" << std::endl;

    const auto&               r_displacement_variables = GetDisplacementVariables();
    std::vector<Dof<double>*> result;
    result.reserve(mNodes.size() * TDim);
    for (const auto& rp_node : mNodes) {
        for (std::size_t d = 0; d < TDim; ++d) {
            result.push_back(rp_node->pGetDof(*r_displacement_variables[d]));
        }
    }
    return result;
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::GetValuesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSolutionStepValues(GetDofs(), Step);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::GetFirstDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractFirstTimeDerivatives(GetDofs(), Step);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    rValues = Geo::DofUtilities::ExtractSecondTimeDerivatives(GetDofs(), Step);
}

template <unsigned int TDim>
Matrix GeoEmbeddedBeamInterfaceElement<TDim>::CalculateLocalRelativeDisplacementMatrix(std::size_t BeamNodeIndex) const
{
    // Eq. (12) and (19): du = u^b_i - N^s(X_i) U^s, expressed in the local axes of the beam
    const auto number_of_dofs                      = mNodes.size() * TDim;
    Matrix     global_relative_displacement_matrix = ZeroMatrix(TDim, number_of_dofs);
    if (mHostInterpolation.IsEmbedded(BeamNodeIndex)) {
        const auto& r_weights = mHostInterpolation.GetInterpolationWeights(BeamNodeIndex);
        for (std::size_t d = 0; d < TDim; ++d) {
            global_relative_displacement_matrix(d, BeamNodeIndex * TDim + d) += 1.0;
            for (std::size_t k = 0; k < r_weights.size(); ++k) {
                global_relative_displacement_matrix(d, mHostNodeIndices[BeamNodeIndex][k] * TDim + d) -=
                    r_weights[k];
            }
        }
    }

    const auto local_axes = EmbeddedBeamUtilities::CalculateLocalAxes(GetGeometry(), TDim);
    Matrix     rotation_matrix(TDim, TDim);
    for (std::size_t i = 0; i < TDim; ++i) {
        for (std::size_t j = 0; j < TDim; ++j) {
            rotation_matrix(i, j) = local_axes(i, j);
        }
    }

    return prod(rotation_matrix, global_relative_displacement_matrix);
}

template <unsigned int TDim>
Vector GeoEmbeddedBeamInterfaceElement<TDim>::CalculateLocalRelativeDisplacement(std::size_t BeamNodeIndex) const
{
    Vector dof_values;
    GetValuesVector(dof_values, 0);
    return prod(CalculateLocalRelativeDisplacementMatrix(BeamNodeIndex), dof_values);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateTractionAndTangent(std::size_t BeamNodeIndex,
                                                                        const Vector& rLocalRelativeDisplacement,
                                                                        Vector& rTraction,
                                                                        Matrix& rTangent,
                                                                        double& rPlasticSlip) const
{
    // Eq. (13)-(14) with a Mohr-Coulomb limit on the tangential component (elastic perfectly plastic),
    // starting from the bonding forces at the start of the stage
    const auto  shear_stiffness    = GetProperties()[INTERFACE_SHEAR_STIFFNESS];
    const auto  normal_stiffness   = GetProperties()[INTERFACE_NORMAL_STIFFNESS];
    const auto& r_initial_traction = mInitialTraction[BeamNodeIndex];

    rTangent = ZeroMatrix(TDim, TDim);
    rTraction.resize(TDim, false);
    for (std::size_t i = 1; i < TDim; ++i) {
        rTangent(i, i) = normal_stiffness;
        rTraction[i]   = r_initial_traction[i] + normal_stiffness * rLocalRelativeDisplacement[i];
    }

    rPlasticSlip                 = mConvergedPlasticSlip[BeamNodeIndex];
    const auto trial_shear_force = r_initial_traction[0] + shear_stiffness * (rLocalRelativeDisplacement[0] -
                                                                              mConvergedPlasticSlip[BeamNodeIndex]);
    const auto maximum_shear_force = mMaximumShearForce[BeamNodeIndex];
    if (std::abs(trial_shear_force) <= maximum_shear_force) {
        rTangent(0, 0) = shear_stiffness;
        rTraction[0]   = trial_shear_force;
    } else {
        const auto sign = trial_shear_force > 0.0 ? 1.0 : -1.0;
        rTraction[0]    = sign * maximum_shear_force;
        rPlasticSlip += sign * (std::abs(trial_shear_force) - maximum_shear_force) / shear_stiffness;
    }
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateAll(MatrixType* pLeftHandSideMatrix,
                                                         VectorType* pRightHandSideVector) const
{
    KRATOS_TRY

    const auto number_of_dofs = mNodes.size() * TDim;
    if (pLeftHandSideMatrix) *pLeftHandSideMatrix = ZeroMatrix(number_of_dofs, number_of_dofs);
    if (pRightHandSideVector) *pRightHandSideVector = ZeroVector(number_of_dofs);

    Vector dof_values;
    GetValuesVector(dof_values, 0);

    // Nodal integration: each beam node represents half of the length of the beam
    const auto integration_coefficient =
        0.5 * EmbeddedBeamUtilities::CalculateReferenceLength(GetGeometry());
    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        if (!mHostInterpolation.IsEmbedded(i)) continue;

        const auto   relative_displacement_matrix = CalculateLocalRelativeDisplacementMatrix(i);
        const Vector relative_displacement        = prod(relative_displacement_matrix, dof_values);

        Vector traction;
        Matrix tangent;
        double plastic_slip = 0.0;
        CalculateTractionAndTangent(i, relative_displacement, traction, tangent, plastic_slip);

        if (pLeftHandSideMatrix) {
            const Matrix auxiliary = prod(tangent, relative_displacement_matrix);
            noalias(*pLeftHandSideMatrix) +=
                integration_coefficient * prod(trans(relative_displacement_matrix), auxiliary);
        }
        if (pRightHandSideVector) {
            noalias(*pRightHandSideVector) -=
                integration_coefficient * prod(trans(relative_displacement_matrix), traction);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                 VectorType& rRightHandSideVector,
                                                                 const ProcessInfo&)
{
    CalculateAll(&rLeftHandSideMatrix, &rRightHandSideVector);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo&)
{
    CalculateAll(&rLeftHandSideMatrix, nullptr);
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                                   const ProcessInfo&)
{
    CalculateAll(nullptr, &rRightHandSideVector);
}

template <unsigned int TDim>
GeometryData::IntegrationMethod GeoEmbeddedBeamInterfaceElement<TDim>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                                                         std::vector<Vector>& rOutput,
                                                                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Other variables (e.g. stresses requested by output or processes that loop over all elements
    // of a model part) are not available for this element and are left to the base class
    if (rVariable != GEO_RELATIVE_DISPLACEMENT_VECTOR && rVariable != GEO_EFFECTIVE_TRACTION_VECTOR) {
        Element::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
        return;
    }

    rOutput.resize(msNumberOfBeamNodes);
    for (std::size_t i = 0; i < msNumberOfBeamNodes; ++i) {
        const auto relative_displacement = CalculateLocalRelativeDisplacement(i);
        if (rVariable == GEO_RELATIVE_DISPLACEMENT_VECTOR) {
            rOutput[i] = relative_displacement;
        } else if (mHostInterpolation.IsEmbedded(i)) {
            Matrix tangent;
            double plastic_slip = 0.0;
            CalculateTractionAndTangent(i, relative_displacement, rOutput[i], tangent, plastic_slip);
        } else {
            rOutput[i] = ZeroVector(TDim);
        }
    }

    KRATOS_CATCH("")
}

template <unsigned int TDim>
int GeoEmbeddedBeamInterfaceElement<TDim>::Check(const ProcessInfo&) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mHostInterpolation.IsInitialized())
        << Info() << " has not been initialized" << std::endl;

    const auto& r_properties = GetProperties();
    for (const auto* p_variable : {&INTERFACE_SHEAR_STIFFNESS, &INTERFACE_NORMAL_STIFFNESS}) {
        KRATOS_ERROR_IF_NOT(r_properties.Has(*p_variable))
            << p_variable->Name() << " is not defined for " << Info() << " (properties "
            << r_properties.Id() << ")" << std::endl;
        KRATOS_ERROR_IF(r_properties[*p_variable] <= 0.0)
            << p_variable->Name() << " of " << Info() << " (properties " << r_properties.Id()
            << ") must be positive" << std::endl;
    }

    if (HasShearStrength()) {
        KRATOS_ERROR_IF(r_properties.Has(GEO_COHESION) && r_properties[GEO_COHESION] < 0.0)
            << "GEO_COHESION of " << Info() << " must not be negative" << std::endl;
        KRATOS_ERROR_IF(r_properties.Has(GEO_FRICTION_ANGLE) && (r_properties[GEO_FRICTION_ANGLE] < 0.0 || r_properties[GEO_FRICTION_ANGLE] >= 90.0))
            << "GEO_FRICTION_ANGLE of " << Info() << " must be in the range [0, 90) degrees" << std::endl;
        if constexpr (TDim == 3) {
            KRATOS_ERROR_IF(!r_properties.Has(DIAMETER) || r_properties[DIAMETER] <= 0.0)
                << "A positive DIAMETER is required to compute the perimeter of " << Info() << std::endl;
        }
    }

    const auto& r_displacement_variables = GetDisplacementVariables();
    for (const auto& rp_node : mNodes) {
        for (std::size_t d = 0; d < TDim; ++d) {
            KRATOS_ERROR_IF_NOT(rp_node->HasDofFor(*r_displacement_variables[d]))
                << "Node " << rp_node->Id() << " of " << Info() << " does not have degree of freedom "
                << r_displacement_variables[d]->Name() << std::endl;
        }
    }

    return 0;

    KRATOS_CATCH("")
}

template <unsigned int TDim>
const EmbeddedHostInterpolation& GeoEmbeddedBeamInterfaceElement<TDim>::GetHostInterpolation() const
{
    return mHostInterpolation;
}

template <unsigned int TDim>
std::string GeoEmbeddedBeamInterfaceElement<TDim>::Info() const
{
    return "GeoEmbeddedBeamInterfaceElement" + std::to_string(TDim) + "D2N #" + std::to_string(Id());
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
}

template <unsigned int TDim>
void GeoEmbeddedBeamInterfaceElement<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
}

template class GeoEmbeddedBeamInterfaceElement<2>;
template class GeoEmbeddedBeamInterfaceElement<3>;

} // namespace Kratos

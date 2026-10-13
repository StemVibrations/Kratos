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

#include "custom_utilities/embedded_beam_utilities.h"
#include "includes/element.h"
#include "includes/kratos_export_api.h"

#include <array>
#include <string>
#include <vector>

namespace Kratos
{

/**
 * @class GeoEmbeddedBeamInterfaceElement
 * @brief Interface between an embedded beam and the soil in which it is embedded, according to
 * Sadek & Shahrour (2004), "A three dimensional embedded beam element for reinforced
 * geomaterials", Int. J. Numer. Anal. Meth. Geomech. 28:931-946, section 2.3.
 * @details This element is defined on the (2-noded) geometry of a beam whose nodes have their
 * own translational degrees of freedom (e.g. a LinearTimoshenkoBeamElement2D2N/3D2N). It connects
 * the translations of the beam nodes to the displacement of the soil at the location of these
 * nodes, which is interpolated from the nodes of the host elements, Eq. (19)-(21). The relative
 * displacement, Eq. (12), is expressed in the local axes of the beam (one tangential and one
 * (2D) or two (3D) normal components) and is linked to the bonding force through, Eq. (13)-(14):
 *     F^int = C^int du^int,   C^int = diag(K_t, K_n, K_n)
 * where K_t (INTERFACE_SHEAR_STIFFNESS) and K_n (INTERFACE_NORMAL_STIFFNESS) are stiffnesses per
 * unit length of the beam. The relative displacements are evaluated at the beam nodes (nodal
 * integration, each node representing half the length of the beam).
 *
 * The tangential bonding force per unit length follows a Mohr-Coulomb criterion:
 *     |t_t| <= P (c + max(0, -sigma'_n) tan(phi))
 * where c is GEO_COHESION, phi is GEO_FRICTION_ANGLE (in degrees), P is the perimeter of the beam
 * (pi * DIAMETER in 3D, two faces per unit width in 2D) and sigma'_n is the effective normal
 * stress of the host element on the shaft (tension positive), averaged over the integration
 * points of the host element and updated at the start of every solution step. When neither
 * GEO_COHESION nor GEO_FRICTION_ANGLE is given, the interface is linear elastic.
 *
 * In a staged analysis, DISPLACEMENT is the displacement of the current stage (it is reset at the
 * start of every stage). Therefore, the converged bonding forces of the previous stage are kept as
 * the initial bonding forces of the current stage, and the plastic slip is measured from the start
 * of the current stage. An element that is inactive in a stage loses its bonding forces.
 *
 * The candidate host elements must be stored in the NEIGHBOUR_ELEMENTS of this element before it
 * is initialized (see FindEmbeddedBeamHostElementsProcess). At beam nodes that are not located
 * inside any host element, no interface forces are generated.
 */
template <unsigned int TDim>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoEmbeddedBeamInterfaceElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoEmbeddedBeamInterfaceElement);

    using Element::CalculateOnIntegrationPoints;

    GeoEmbeddedBeamInterfaceElement() = default;
    GeoEmbeddedBeamInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry);
    GeoEmbeddedBeamInterfaceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType               NewId,
                            GeometryType::Pointer   pGeometry,
                            PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetValuesVector(Vector& rValues, int Step = 0) const override;
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /// Results are given at the two beam nodes, in the order of the nodes of the geometry
    [[nodiscard]] IntegrationMethod GetIntegrationMethod() const override;

    /// Supports GEO_RELATIVE_DISPLACEMENT_VECTOR and GEO_EFFECTIVE_TRACTION_VECTOR (force per
    /// unit length), both in the local axes of the beam [tangential, normal(s)]
    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    [[nodiscard]] const EmbeddedHostInterpolation& GetHostInterpolation() const;

    [[nodiscard]] std::string Info() const override;
    void                      PrintInfo(std::ostream& rOStream) const override;

private:
    static constexpr std::size_t msNumberOfBeamNodes         = 2;
    static constexpr double      msLocalCoordinatesTolerance = 1.0e-6;

    [[nodiscard]] std::vector<Dof<double>*> GetDofs() const;
    [[nodiscard]] Matrix CalculateLocalRelativeDisplacementMatrix(std::size_t BeamNodeIndex) const;
    [[nodiscard]] Vector CalculateLocalRelativeDisplacement(std::size_t BeamNodeIndex) const;
    void                 CalculateTractionAndTangent(std::size_t   BeamNodeIndex,
                                                     const Vector& rLocalRelativeDisplacement,
                                                     Vector&       rTraction,
                                                     Matrix&       rTangent,
                                                     double&       rPlasticSlip) const;
    void CalculateAll(MatrixType* pLeftHandSideMatrix, VectorType* pRightHandSideVector) const;
    void UpdateMaximumShearForces(const ProcessInfo& rCurrentProcessInfo, bool UseHostStresses);
    [[nodiscard]] double CalculateNormalStressOnShaft(std::size_t        BeamNodeIndex,
                                                      const ProcessInfo& rCurrentProcessInfo) const;
    [[nodiscard]] bool   HasShearStrength() const;
    [[nodiscard]] double GetPerimeter() const;

    EmbeddedHostInterpolation mHostInterpolation;

    /// Nodes providing the degrees of freedom: the beam nodes, followed by the host nodes that are not beam nodes
    std::vector<Node::Pointer>                                mNodes;
    std::array<std::vector<std::size_t>, msNumberOfBeamNodes> mHostNodeIndices;
    std::array<double, msNumberOfBeamNodes>                   mConvergedPlasticSlip = {0.0, 0.0};
    std::array<double, msNumberOfBeamNodes>                   mMaximumShearForce    = {0.0, 0.0};

    /// Bonding forces (per unit length, local axes) at the start of the stage, i.e. the converged
    /// bonding forces of the previous stage
    std::array<Vector, msNumberOfBeamNodes> mInitialTraction;
    /// Bonding forces (per unit length, local axes) of the last converged solution step
    std::array<Vector, msNumberOfBeamNodes> mConvergedTraction;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos

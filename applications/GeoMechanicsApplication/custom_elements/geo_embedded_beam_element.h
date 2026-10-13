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

#include "custom_elements/beam_elements/timoshenko_beam_element_3D2N.h"
#include "custom_utilities/embedded_beam_utilities.h"
#include "includes/kratos_export_api.h"

#include <string>
#include <type_traits>
#include <vector>

namespace Kratos
{

template <unsigned int TDim>
using GeoEmbeddedBeamBaseElement =
    std::conditional_t<TDim == 2, LinearTimoshenkoBeamElement2D2N, LinearTimoshenkoBeamElement3D2N>;

/**
 * @class GeoEmbeddedBeamElement
 * @brief Perfectly bonded embedded beam element according to Sadek & Shahrour (2004), "A three
 * dimensional embedded beam element for reinforced geomaterials", Int. J. Numer. Anal. Meth.
 * Geomech. 28:931-946, section 2.2.
 * @details The nodes of the beam do not need to coincide with the nodes of the soil mesh. The
 * translations of a beam node are interpolated from the displacements of the nodes of the
 * (solid) element in which it is located, Eq. (5), whereas the rotations remain independent
 * degrees of freedom of the beam nodes. With the transformation matrix N^bs of Eq. (7), the beam
 * stiffness matrix is expressed in terms of the host-node displacements and the beam-node
 * rotations, Eq. (11):
 *     R^bs = N^bs^T R^b N^bs
 * The same transformation is applied to the internal and external force vectors, the mass matrix
 * and the damping matrix. The degrees of freedom of the element are ordered as
 *     [U^s (translations of all host nodes), Theta_i, Theta_j]
 * The beam itself is a linear Timoshenko beam (LinearTimoshenkoBeamElement2D2N/3D2N of the
 * StructuralMechanicsApplication).
 *
 * The candidate host elements must be stored in the NEIGHBOUR_ELEMENTS of this element before it
 * is initialized (see FindEmbeddedBeamHostElementsProcess). A beam node that is not located inside
 * any of these elements keeps its own translational degrees of freedom. Since the translations of
 * an embedded beam node are not independent degrees of freedom, point loads and boundary
 * conditions on translations must be applied to the soil (or to a beam node that coincides with a
 * soil node).
 */
template <unsigned int TDim>
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoEmbeddedBeamElement : public GeoEmbeddedBeamBaseElement<TDim>
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GeoEmbeddedBeamElement);

    using BaseType = GeoEmbeddedBeamBaseElement<TDim>;
    using BaseType::CalculateOnIntegrationPoints;
    using typename BaseType::DofsVectorType;
    using typename BaseType::EquationIdVectorType;
    using typename BaseType::GeometryType;
    using typename BaseType::IndexType;
    using typename BaseType::MatrixType;
    using typename BaseType::NodesArrayType;
    using typename BaseType::PropertiesType;
    using typename BaseType::SizeType;
    using typename BaseType::VectorType;

    GeoEmbeddedBeamElement() = default;
    GeoEmbeddedBeamElement(IndexType NewId, typename GeometryType::Pointer pGeometry);
    GeoEmbeddedBeamElement(IndexType                        NewId,
                           typename GeometryType::Pointer   pGeometry,
                           typename PropertiesType::Pointer pProperties);

    Element::Pointer Create(IndexType                        NewId,
                            const NodesArrayType&            rNodes,
                            typename PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType                        NewId,
                            typename GeometryType::Pointer   pGeometry,
                            typename PropertiesType::Pointer pProperties) const override;
    Element::Pointer Clone(IndexType NewId, const NodesArrayType& rNodes) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;
    void GetValuesVector(Vector& rValues, int Step = 0) const override;
    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;
    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    /// Returns the beam degrees of freedom in local axes, where the translations are interpolated from the host nodes
    void GetNodalValuesVector(VectorType& rNodalValues) const override;
    void RotateLHS(MatrixType& rLHS, const GeometryType& rGeometry) override;
    void RotateRHS(VectorType& rRHS, const GeometryType& rGeometry) override;
    void RotateAll(MatrixType& rLHS, VectorType& rRHS, const GeometryType& rGeometry) override;

    /// The body forces are added by this element itself (see CalculateBodyForceVector), hence the
    /// base element should not add them
    array_1d<double, 3> GetLocalAxesBodyForce(const Element&,
                                              const typename GeometryType::IntegrationPointsArrayType&,
                                              IndexType) const override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    /// Returns the transformation matrix N^bs (Eq. 7) that maps the degrees of freedom of this
    /// element to the (global) degrees of freedom of the beam [u_i, theta_i, u_j, theta_j]
    [[nodiscard]] Matrix                           CalculateTransformationMatrix() const;
    [[nodiscard]] const EmbeddedHostInterpolation& GetHostInterpolation() const;

    [[nodiscard]] std::string Info() const override;
    void                      PrintInfo(std::ostream& rOStream) const override;

private:
    static constexpr SizeType msNumberOfBeamNodes         = 2;
    static constexpr SizeType msNumberOfRotations         = TDim == 2 ? 1 : 3;
    static constexpr SizeType msDofsPerBeamNode           = TDim + msNumberOfRotations;
    static constexpr SizeType msNumberOfBeamDofs          = msNumberOfBeamNodes * msDofsPerBeamNode;
    static constexpr double   msLocalCoordinatesTolerance = 1.0e-6;

    [[nodiscard]] std::vector<Dof<double>*> GetDofs() const;
    [[nodiscard]] BoundedMatrix<double, msNumberOfBeamDofs, msNumberOfBeamDofs> CalculateGlobalToLocalRotationMatrix() const;
    [[nodiscard]] Vector CalculateBodyForceVector();
    [[nodiscard]] Matrix CalculateLumpedBeamMassMatrix();
    [[nodiscard]] double GetBeamCrossArea();
    void                 UpdateEmbeddedBeamNodalKinematics();

    EmbeddedHostInterpolation mHostInterpolation;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};

} // namespace Kratos

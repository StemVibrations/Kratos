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

#include "custom_processes/deactivate_unsupported_interfaces_process.h"
#include "custom_elements/U_Pw_interface_element.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

#include <algorithm>
#include <unordered_map>
#include <vector>

namespace
{

using namespace Kratos;

bool IsInterfaceElement(const Element& rElement)
{
    return dynamic_cast<const UPwInterfaceElement*>(&rElement) != nullptr;
}

std::size_t GetNumberOfCornerNodes(const Geometry<Node>& rSideGeometry)
{
    switch (rSideGeometry.GetGeometryFamily()) {
        using enum GeometryData::KratosGeometryFamily;
    case Kratos_Linear:
        return 2;
    case Kratos_Triangle:
        return 3;
    case Kratos_Quadrilateral:
        return 4;
    default:
        return rSideGeometry.PointsNumber();
    }
}

bool Contains(const Element& rElement, std::size_t NodeId)
{
    const auto& r_geometry = rElement.GetGeometry();
    return std::any_of(r_geometry.begin(), r_geometry.end(),
                       [NodeId](const auto& rNode) { return rNode.Id() == NodeId; });
}

} // namespace

namespace Kratos
{
using namespace std::string_literals;

DeactivateUnsupportedInterfacesProcess::DeactivateUnsupportedInterfacesProcess(ModelPart& rModelPart)
    : mrModelPart(rModelPart)
{
}

void DeactivateUnsupportedInterfacesProcess::Execute()
{
    KRATOS_TRY

    // The non-interface elements that contain a node, split into active and inactive elements
    using ElementsByNodeId          = std::unordered_map<std::size_t, std::vector<const Element*>>;
    auto active_elements_by_node_id   = ElementsByNodeId{};
    auto inactive_elements_by_node_id = ElementsByNodeId{};
    for (const auto& r_element : mrModelPart.Elements()) {
        if (IsInterfaceElement(r_element)) continue;

        auto& r_elements_by_node_id =
            r_element.IsActive() ? active_elements_by_node_id : inactive_elements_by_node_id;
        for (const auto& r_node : r_element.GetGeometry()) {
            r_elements_by_node_id[r_node.Id()].push_back(&r_element);
        }
    }

    mNumberOfDeactivatedInterfaces = 0;
    for (auto& r_element : mrModelPart.Elements()) {
        const auto* p_interface_element = dynamic_cast<const UPwInterfaceElement*>(&r_element);
        if (!p_interface_element || !r_element.IsActive()) continue;

        const auto& r_geometry = r_element.GetGeometry();
        const auto  number_of_nodes_per_side = r_geometry.PointsNumber() / 2;
        const auto  number_of_corner_nodes =
            GetNumberOfCornerNodes(p_interface_element->GetDisplacementMidGeometry());

        // Whether one of the given elements contains all corner nodes of a side. Note that the nodes
        // of a side are ordered like the nodes of the mid geometry, i.e. the corner nodes come first.
        const auto has_element_at_side = [&](const ElementsByNodeId& rElementsByNodeId, std::size_t IndexOfFirstNode) {
            const auto it = rElementsByNodeId.find(r_geometry[IndexOfFirstNode].Id());
            if (it == rElementsByNodeId.end()) return false;

            return std::any_of(it->second.begin(), it->second.end(), [&](const Element* pCandidate) {
                for (auto i = std::size_t{1}; i < number_of_corner_nodes; ++i) {
                    if (!Contains(*pCandidate, r_geometry[IndexOfFirstNode + i].Id())) return false;
                }
                return true;
            });
        };

        // A side lost its support when its adjacent elements have all been deactivated. Sides that
        // never had an adjacent element (e.g. which are supported by boundary conditions) are kept.
        const auto has_lost_support = [&](std::size_t IndexOfFirstNode) {
            return !has_element_at_side(active_elements_by_node_id, IndexOfFirstNode) &&
                   has_element_at_side(inactive_elements_by_node_id, IndexOfFirstNode);
        };

        if (has_lost_support(0) || has_lost_support(number_of_nodes_per_side)) {
            r_element.Set(ACTIVE, false);
            ++mNumberOfDeactivatedInterfaces;
        }
    }

    KRATOS_CATCH("")
}

std::size_t DeactivateUnsupportedInterfacesProcess::GetNumberOfDeactivatedInterfaces() const
{
    return mNumberOfDeactivatedInterfaces;
}

std::string DeactivateUnsupportedInterfacesProcess::Info() const
{
    return "DeactivateUnsupportedInterfacesProcess"s;
}

} // namespace Kratos

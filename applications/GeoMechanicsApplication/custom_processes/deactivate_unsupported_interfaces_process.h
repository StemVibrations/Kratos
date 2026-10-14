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

#include "includes/kratos_export_api.h"
#include "processes/process.h"

#include <cstddef>
#include <string>

namespace Kratos
{

class ModelPart;

/**
 * @class DeactivateUnsupportedInterfacesProcess
 * @ingroup GeoMechanicsApplication
 * @brief Deactivates the active interface elements of which a side is no longer supported by an
 * active (non-interface) element, e.g. after the soil at one side of a wall has been excavated.
 * @details A side has lost its support when no active non-interface element contains all corner
 * nodes of that side, whereas an inactive one does. Such an interface can't transfer any traction:
 * the nodes at its unsupported side are only held by the interface itself. Keeping it active makes
 * the problem ill-posed, e.g. when the water pressure in the interface pushes the unsupported side
 * away from the other side. Sides without any adjacent element (e.g. which are only supported by
 * boundary conditions) are left untouched.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) DeactivateUnsupportedInterfacesProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(DeactivateUnsupportedInterfacesProcess);

    explicit DeactivateUnsupportedInterfacesProcess(ModelPart& rModelPart);
    DeactivateUnsupportedInterfacesProcess(const DeactivateUnsupportedInterfacesProcess&) = delete;
    DeactivateUnsupportedInterfacesProcess& operator=(const DeactivateUnsupportedInterfacesProcess&) = delete;
    ~DeactivateUnsupportedInterfacesProcess() override = default;

    void Execute() override;

    [[nodiscard]] std::size_t GetNumberOfDeactivatedInterfaces() const;
    [[nodiscard]] std::string Info() const override;

private:
    ModelPart&  mrModelPart;
    std::size_t mNumberOfDeactivatedInterfaces = 0;
};

} // namespace Kratos

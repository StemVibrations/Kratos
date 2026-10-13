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

#include <functional>
#include <string>
#include <vector>

namespace Kratos
{

class Model;
class ModelPart;
class Parameters;

/**
 * @brief Finds, for every node of the embedded beam (interface) elements, the (solid) host element
 * in which it is located and stores these host elements in the NEIGHBOUR_ELEMENTS of the embedded
 * elements. This needs to be done before the embedded elements are initialized.
 * @details Settings:
 *   - "model_part_name" or "model_part_name_list": model part(s) with the embedded elements
 *   - "host_model_part_name" or "host_model_part_name_list": model part(s) with the candidate host
 *     elements. Only elements whose local dimension equals their working space dimension (i.e.
 *     solid elements) are considered.
 *   - "search_tolerance" (optional, default 1e-6): tolerance on the local coordinates used to
 *     decide whether a point is inside a host element.
 * Embedded nodes that are not located inside any host element are reported. Such nodes keep their
 * own translational degrees of freedom.
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) FindEmbeddedBeamHostElementsProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(FindEmbeddedBeamHostElementsProcess);

    FindEmbeddedBeamHostElementsProcess(Model& rModel, const Parameters& rProcessSettings);

    void ExecuteInitialize() override;

    [[nodiscard]] std::string Info() const override;

private:
    std::vector<std::reference_wrapper<ModelPart>> mrEmbeddedModelParts;
    std::vector<std::reference_wrapper<ModelPart>> mrHostModelParts;
    double                                         mSearchTolerance = 1.0e-6;
};

} // namespace Kratos

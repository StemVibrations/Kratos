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

#include "find_embedded_beam_host_elements_process.h"
#include "containers/model.h"
#include "custom_utilities/process_utilities.h"
#include "includes/element.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <set>
#include <sstream>
#include <unordered_set>

using namespace std::string_literals;

namespace
{

using namespace Kratos;

std::vector<std::string> GetHostModelPartNames(const Parameters& rProcessSettings)
{
    KRATOS_ERROR_IF_NOT(rProcessSettings.Has("host_model_part_name") ||
                        rProcessSettings.Has("host_model_part_name_list"))
        << "Please specify 'host_model_part_name' or 'host_model_part_name_list' for "
           "FindEmbeddedBeamHostElementsProcess"
        << std::endl;

    KRATOS_ERROR_IF(rProcessSettings.Has("host_model_part_name") &&
                    rProcessSettings.Has("host_model_part_name_list"))
        << "The parameters 'host_model_part_name' and 'host_model_part_name_list' are mutually "
           "exclusive for FindEmbeddedBeamHostElementsProcess"
        << std::endl;

    auto result = rProcessSettings.Has("host_model_part_name_list")
                      ? rProcessSettings["host_model_part_name_list"].GetStringArray()
                      : std::vector{rProcessSettings["host_model_part_name"].GetString()};
    KRATOS_ERROR_IF(result.empty())
        << "'host_model_part_name_list' needs to contain at least one model part name for "
           "FindEmbeddedBeamHostElementsProcess"
        << std::endl;
    return result;
}

/// A uniform grid of the (padded) axis-aligned bounding boxes of elements, used to efficiently
/// find the element in which a point is located
class HostElementGrid
{
public:
    HostElementGrid(std::vector<Element*> Elements, double RelativeTolerance)
        : mElements(std::move(Elements))
    {
        mLowerCorners.reserve(mElements.size());
        mUpperCorners.reserve(mElements.size());

        array_1d<double, 3> mean_size = ZeroVector(3);
        for (const auto* p_element : mElements) {
            const auto&         r_geometry = p_element->GetGeometry();
            array_1d<double, 3> lower      = r_geometry[0].Coordinates();
            array_1d<double, 3> upper      = lower;
            for (const auto& r_node : r_geometry) {
                for (std::size_t d = 0; d < 3; ++d) {
                    lower[d] = std::min(lower[d], r_node.Coordinates()[d]);
                    upper[d] = std::max(upper[d], r_node.Coordinates()[d]);
                }
            }

            const auto padding = RelativeTolerance * norm_2(upper - lower) + 1.0e-12;
            for (std::size_t d = 0; d < 3; ++d) {
                lower[d] -= padding;
                upper[d] += padding;
            }
            noalias(mean_size) += upper - lower;
            mLowerCorners.push_back(lower);
            mUpperCorners.push_back(upper);
        }

        if (mElements.empty()) return;

        mean_size /= static_cast<double>(mElements.size());
        mLower = mLowerCorners.front();
        mUpper = mUpperCorners.front();
        for (std::size_t i = 0; i < mElements.size(); ++i) {
            for (std::size_t d = 0; d < 3; ++d) {
                mLower[d] = std::min(mLower[d], mLowerCorners[i][d]);
                mUpper[d] = std::max(mUpper[d], mUpperCorners[i][d]);
            }
        }

        // Use cells with the size of an average element, while limiting the total number of cells
        const auto maximum_number_of_cells = std::max<std::size_t>(1000, 8 * mElements.size());
        for (std::size_t d = 0; d < 3; ++d) {
            mCellSize[d] = std::max(mean_size[d], 1.0e-12);
        }
        while (true) {
            std::size_t number_of_cells = 1;
            for (std::size_t d = 0; d < 3; ++d) {
                mNumberOfCells[d] = std::max<std::size_t>(
                    1, static_cast<std::size_t>(std::ceil((mUpper[d] - mLower[d]) / mCellSize[d])));
                number_of_cells *= mNumberOfCells[d];
            }
            if (number_of_cells <= maximum_number_of_cells) break;
            for (auto& r_cell_size : mCellSize) {
                r_cell_size *= 2.0;
            }
        }

        mCells.resize(mNumberOfCells[0] * mNumberOfCells[1] * mNumberOfCells[2]);
        for (std::size_t element_index = 0; element_index < mElements.size(); ++element_index) {
            std::array<std::size_t, 3> first{};
            std::array<std::size_t, 3> last{};
            for (std::size_t d = 0; d < 3; ++d) {
                first[d] = CellCoordinate(mLowerCorners[element_index][d], d);
                last[d]  = CellCoordinate(mUpperCorners[element_index][d], d);
            }
            for (auto i = first[0]; i <= last[0]; ++i) {
                for (auto j = first[1]; j <= last[1]; ++j) {
                    for (auto k = first[2]; k <= last[2]; ++k) {
                        mCells[CellIndex({i, j, k})].push_back(element_index);
                    }
                }
            }
        }
    }

    [[nodiscard]] Element* FindContainingElement(const array_1d<double, 3>& rPoint, double LocalTolerance) const
    {
        if (mElements.empty()) return nullptr;
        for (std::size_t d = 0; d < 3; ++d) {
            if (rPoint[d] < mLower[d] || rPoint[d] > mUpper[d]) return nullptr;
        }

        const auto cell_index = CellIndex({CellCoordinate(rPoint[0], 0), CellCoordinate(rPoint[1], 1),
                                           CellCoordinate(rPoint[2], 2)});
        array_1d<double, 3> local_coordinates;
        for (const auto element_index : mCells[cell_index]) {
            if (!BoundingBoxContains(element_index, rPoint)) continue;
            if (mElements[element_index]->GetGeometry().IsInside(rPoint, local_coordinates, LocalTolerance)) {
                return mElements[element_index];
            }
        }
        return nullptr;
    }

private:
    [[nodiscard]] std::size_t CellCoordinate(double Coordinate, std::size_t Direction) const
    {
        const auto position = (Coordinate - mLower[Direction]) / mCellSize[Direction];
        if (position <= 0.0) return 0;
        return std::min(static_cast<std::size_t>(position), mNumberOfCells[Direction] - 1);
    }

    [[nodiscard]] std::size_t CellIndex(const std::array<std::size_t, 3>& rCellCoordinates) const
    {
        return (rCellCoordinates[2] * mNumberOfCells[1] + rCellCoordinates[1]) * mNumberOfCells[0] +
               rCellCoordinates[0];
    }

    [[nodiscard]] bool BoundingBoxContains(std::size_t ElementIndex, const array_1d<double, 3>& rPoint) const
    {
        for (std::size_t d = 0; d < 3; ++d) {
            if (rPoint[d] < mLowerCorners[ElementIndex][d] || rPoint[d] > mUpperCorners[ElementIndex][d])
                return false;
        }
        return true;
    }

    std::vector<Element*>                 mElements;
    std::vector<array_1d<double, 3>>      mLowerCorners;
    std::vector<array_1d<double, 3>>      mUpperCorners;
    array_1d<double, 3>                   mLower = ZeroVector(3);
    array_1d<double, 3>                   mUpper = ZeroVector(3);
    std::array<double, 3>                 mCellSize{};
    std::array<std::size_t, 3>            mNumberOfCells{1, 1, 1};
    std::vector<std::vector<std::size_t>> mCells;
};

} // namespace

namespace Kratos
{

FindEmbeddedBeamHostElementsProcess::FindEmbeddedBeamHostElementsProcess(Model& rModel, const Parameters& rProcessSettings)
    : mrEmbeddedModelParts{ProcessUtilities::GetModelPartsFromSettings(
          rModel, rProcessSettings, FindEmbeddedBeamHostElementsProcess::Info())}
{
    for (const auto& r_name : GetHostModelPartNames(rProcessSettings)) {
        mrHostModelParts.emplace_back(rModel.GetModelPart(r_name));
    }

    if (rProcessSettings.Has("search_tolerance")) {
        mSearchTolerance = rProcessSettings["search_tolerance"].GetDouble();
    }
    KRATOS_ERROR_IF(mSearchTolerance < 0.0)
        << "The search tolerance of " << FindEmbeddedBeamHostElementsProcess::Info()
        << " must not be negative" << std::endl;
}

void FindEmbeddedBeamHostElementsProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Only solid elements can host an embedded beam
    std::vector<Element*>         candidate_host_elements;
    std::unordered_set<IndexType> candidate_ids;
    for (const auto& r_host_model_part : mrHostModelParts) {
        for (auto& r_element : r_host_model_part.get().Elements()) {
            const auto& r_geometry = r_element.GetGeometry();
            if (r_geometry.LocalSpaceDimension() != r_geometry.WorkingSpaceDimension()) continue;
            if (candidate_ids.insert(r_element.Id()).second)
                candidate_host_elements.push_back(&r_element);
        }
    }
    KRATOS_ERROR_IF(candidate_host_elements.empty())
        << "No solid elements were found in the host model part(s) of " << Info() << std::endl;

    const HostElementGrid grid(std::move(candidate_host_elements), mSearchTolerance);

    std::set<IndexType> free_node_ids;
    std::size_t         number_of_embedded_elements = 0;
    for (const auto& r_embedded_model_part : mrEmbeddedModelParts) {
        for (auto& r_element : r_embedded_model_part.get().Elements()) {
            GlobalPointersVector<Element> host_elements;
            for (const auto& r_node : r_element.GetGeometry()) {
                auto* p_host_element = grid.FindContainingElement(r_node.Coordinates(), mSearchTolerance);
                if (!p_host_element) {
                    free_node_ids.insert(r_node.Id());
                    continue;
                }

                const auto is_known_host = std::any_of(
                    host_elements.ptr_begin(), host_elements.ptr_end(),
                    [p_host_element](const auto& rpHost) { return rpHost.get() == p_host_element; });
                if (!is_known_host) host_elements.push_back(GlobalPointer<Element>(p_host_element));
            }
            r_element.SetValue(NEIGHBOUR_ELEMENTS, host_elements);
            ++number_of_embedded_elements;
        }
    }

    KRATOS_INFO(Info()) << "Assigned host elements to " << number_of_embedded_elements
                        << " embedded element(s)" << std::endl;

    if (!free_node_ids.empty()) {
        std::ostringstream node_ids;
        for (const auto id : free_node_ids) {
            node_ids << " " << id;
        }
        KRATOS_WARNING(Info()) << free_node_ids.size()
                               << " node(s) of embedded elements are not located inside any host "
                                  "element. These nodes keep their own translational degrees of "
                                  "freedom. Node Ids:"
                               << node_ids.str() << std::endl;
    }

    KRATOS_CATCH("")
}

std::string FindEmbeddedBeamHostElementsProcess::Info() const
{
    return "FindEmbeddedBeamHostElementsProcess"s;
}

} // namespace Kratos

// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#include "apply_constant_boundary_hydrostatic_pressure_process.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/model_part.h"

namespace Kratos
{
using namespace std::string_literals;

ApplyConstantBoundaryHydrostaticPressureProcess::ApplyConstantBoundaryHydrostaticPressureProcess(ModelPart& model_part,
                                                                                                 Parameters rParameters)
    : Process(Flags()), mrModelPart(model_part)
{
    KRATOS_TRY

    // only include validation with c++11 since raw_literals do not exist in c++03
    Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "gravity_direction" : 2,
                "reference_coordinate" : [0.0, 0.0, 0.0],
                "specific_weight" : 10000.0,
                "table" : 1
            }  )");

    // Some values have to be input by the user since no meaningful default value exist. For
    // this reason, we try to access them, so that an error is thrown if they don't exist.
    rParameters["reference_coordinate"];
    rParameters["variable_name"];
    rParameters["model_part_name"];

    mIsFixedProvided = rParameters.Has("is_fixed");

    // Now validate against defaults -- this also ensures no type mismatch
    rParameters.ValidateAndAssignDefaults(default_parameters);

    mVariableName        = rParameters["variable_name"].GetString();
    mIsFixed             = rParameters["is_fixed"].GetBool();
    mGravityDirection    = static_cast<unsigned int>(rParameters["gravity_direction"].GetInt());
    mReferenceCoordinate = rParameters["reference_coordinate"].GetVector();
    mSpecificWeight      = rParameters["specific_weight"].GetDouble();

    KRATOS_CATCH("")
}

void ApplyConstantBoundaryHydrostaticPressureProcess::ExecuteInitialize()
{
    KRATOS_TRY

    const auto& r_variable = KratosComponents<Variable<double>>::Get(GetVariableName());

    block_for_each(GetModelPart().Nodes(), [&r_variable, this](Node& rNode) {
        if (mIsFixed) rNode.Fix(r_variable);
        else if (mIsFixedProvided) rNode.Free(r_variable);

        const double pressure = GetSpecificWeight() *
                              (GetReferenceCoordinate()[GetGravityDirection()] - rNode.Coordinates()[GetGravityDirection()]);
        rNode.FastGetSolutionStepValue(r_variable) = std::max(pressure, 0.);
    });

    block_for_each(GetModelPart().Conditions(), [&r_variable, this](Condition& rCondition) {

        const auto& r_geometry = rCondition.GetGeometry();

        // Edge vector
        Vector edge = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();

		bool is_edge_vertical = false;
        if (abs(edge[0]) < 1e-12)
        {
			is_edge_vertical = true;
        }

        // Out-of-plane
        Vector out_of_plane_vector = ZeroVector(3);
        out_of_plane_vector[2] = 1.0;

        // Normal
        Vector normal_vector = ZeroVector(3);
        MathUtils<double>::CrossProduct(normal_vector, edge, out_of_plane_vector);

        const double norm = norm_2(normal_vector);
        if (norm > 1e-12) {
            normal_vector /= norm;
        }

        // Center of condition
        Vector center = ZeroVector(3);
        for (SizeType i = 0; i < r_geometry.size(); ++i) {
            center += r_geometry[i].Coordinates();
        }
        center /= r_geometry.size();

        // Water reference point (must be defined!)
        auto& reference_coordinate = GetReferenceCoordinate();

        Vector water_to_element = center - reference_coordinate;

		double sign = 0.0;
        if (!is_edge_vertical)
        {
            if (normal_vector[1] < 0.0)
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }
        }
        else
        {
			if (water_to_element[0] * normal_vector[0])
			{
				sign = -1.0;
			}
			else
			{
				sign = 1.0;
			}
        }


        rCondition.SetValue(CONTACT_PRESSURE, sign);

       
        });


    KRATOS_CATCH("")
}

std::string ApplyConstantBoundaryHydrostaticPressureProcess::Info() const
{
    return "ApplyConstantBoundaryHydrostaticPressureProcess"s;
}

ModelPart& ApplyConstantBoundaryHydrostaticPressureProcess::GetModelPart() { return mrModelPart; }

const std::string& ApplyConstantBoundaryHydrostaticPressureProcess::GetVariableName() const
{
    return mVariableName;
}

unsigned int ApplyConstantBoundaryHydrostaticPressureProcess::GetGravityDirection() const
{
    return mGravityDirection;
}

const Vector& ApplyConstantBoundaryHydrostaticPressureProcess::GetReferenceCoordinate() const
{
    return mReferenceCoordinate;
}

double ApplyConstantBoundaryHydrostaticPressureProcess::GetSpecificWeight() const
{
    return mSpecificWeight;
}

} // namespace Kratos
import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication import geomechanics_analysis as analysis
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers

import test_helper


def get_column_under_gravity_directory():
    return test_helper.get_file_path(
        os.path.join('test_mohr_coulomb_with_tension_cutoff', 'test_column_under_gravity'))


def read_parameters(test_directory):
    """
    Reads the project parameters of a test case. The output processes are removed, such that the
    (reference) output of the test case is left untouched.
    """
    with open(os.path.join(test_directory, 'ProjectParameters.json'), 'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    if parameters.Has("output_processes"):
        parameters.RemoveValue("output_processes")

    return parameters


def use_global_force_error_criterion(parameters, criterion_settings='{"echo_level": 1}', use_local_error_criteria=False):
    solver_settings = parameters["solver_settings"]
    solver_settings["convergence_criterion"].SetString("global_force_error_criterion")
    for name in ("global_force_error_criterion_settings", "use_local_error_criteria"):
        if solver_settings.Has(name):
            solver_settings.RemoveValue(name)
    solver_settings.AddValue("global_force_error_criterion_settings", Kratos.Parameters(criterion_settings))
    solver_settings.AddBool("use_local_error_criteria", use_local_error_criteria)
    return parameters


def make_next_stage(parameters):
    """Returns the parameters of a stage that follows the given one, without any changes of the model."""
    result = parameters.Clone()
    problem_data = result["problem_data"]
    duration = problem_data["end_time"].GetDouble() - problem_data["start_time"].GetDouble()
    problem_data["start_time"].SetDouble(problem_data["end_time"].GetDouble())
    problem_data["end_time"].SetDouble(problem_data["end_time"].GetDouble() + duration)
    model_import_settings = result["solver_settings"]["model_import_settings"]
    model_import_settings.RemoveValue("input_filename")
    model_import_settings["input_type"].SetString("use_input_model_part")
    return result


def run_stages(test_directory, stage_parameters):
    model = Kratos.Model()
    simulations = []
    with context_managers.set_cwd_to(test_directory):
        for parameters in stage_parameters:
            simulations.append(analysis.GeoMechanicsAnalysis(model, parameters))
            simulations[-1].Run()

    return simulations


def get_computing_model_part(simulation):
    return simulation._GetSolver().GetComputingModelPart()


def get_stress_vectors(simulation):
    model_part = get_computing_model_part(simulation)
    return [stress_vector for element in model_part.Elements
            for stress_vector in element.CalculateOnIntegrationPoints(Kratos.CAUCHY_STRESS_VECTOR,
                                                                      model_part.ProcessInfo)]


class KratosGeoMechanicsGlobalForceErrorCriterionTests(KratosUnittest.TestCase):
    """
    Checks the global force error criterion in combination with the Newton-Raphson strategy of the
    GeoMechanicsApplication.
    """

    def assert_stresses_close(self, simulation, reference_simulation, relative_tolerance):
        reference_stresses = get_stress_vectors(reference_simulation)
        stresses = get_stress_vectors(simulation)
        stress_scale = max(abs(value) for stress_vector in reference_stresses for value in stress_vector)
        self.assertEqual(len(stresses), len(reference_stresses))
        for reference_stress_vector, stress_vector in zip(reference_stresses, stresses):
            for reference_value, value in zip(reference_stress_vector, stress_vector):
                self.assertAlmostEqual(value, reference_value, delta=relative_tolerance * stress_scale)

    def test_strict_criterion_reproduces_reference_results(self):
        test_directory = get_column_under_gravity_directory()
        reference = run_stages(test_directory, [read_parameters(test_directory)])[0]

        parameters = use_global_force_error_criterion(read_parameters(test_directory),
                                                      '{"force_relative_tolerance": 1.0e-6, "force_absolute_tolerance": 0.0}')
        simulation = run_stages(test_directory, [parameters])[0]

        self.assertIsNotNone(simulation._GetSolver().global_force_error_criterion)
        # The reference results are converged up to a relative residual of 1.0e-4
        self.assert_stresses_close(simulation, reference, 1.0e-4)

    def test_default_force_relative_tolerance_gives_results_within_tolerance(self):
        test_directory = get_column_under_gravity_directory()
        reference = run_stages(test_directory, [read_parameters(test_directory)])[0]

        simulation = run_stages(test_directory, [use_global_force_error_criterion(read_parameters(test_directory))])[0]

        self.assert_stresses_close(simulation, reference, 0.01)

    def test_unchanged_stage_converges_in_first_iteration(self):
        # The loads of the first stage are inactive in the second stage. Since these are already in
        # equilibrium, the second stage converges immediately.
        test_directory = get_column_under_gravity_directory()
        first_stage = use_global_force_error_criterion(read_parameters(test_directory))
        simulations = run_stages(test_directory, [first_stage, make_next_stage(first_stage)])

        process_info = get_computing_model_part(simulations[1]).ProcessInfo
        self.assertEqual(process_info[Kratos.NL_ITERATION_NUMBER], 1)
        self.assert_stresses_close(simulations[1], simulations[0], 0.01)

    def test_combination_with_local_error_criteria(self):
        test_directory = get_column_under_gravity_directory()
        reference = run_stages(test_directory, [read_parameters(test_directory)])[0]

        parameters = use_global_force_error_criterion(read_parameters(test_directory), use_local_error_criteria=True)
        simulation = run_stages(test_directory, [parameters])[0]

        self.assertEqual(type(simulation._GetSolver().convergence_criterion).__name__, "AndCriteria")
        self.assert_stresses_close(simulation, reference, 0.01)

    def test_criterion_requires_newton_raphson_strategy(self):
        test_directory = get_column_under_gravity_directory()
        parameters = use_global_force_error_criterion(read_parameters(test_directory))
        parameters["solver_settings"]["strategy_type"].SetString("linear")

        with self.assertRaisesRegex(RuntimeError, "can't be used with strategy type"):
            run_stages(test_directory, [parameters])


if __name__ == '__main__':
    KratosUnittest.main()

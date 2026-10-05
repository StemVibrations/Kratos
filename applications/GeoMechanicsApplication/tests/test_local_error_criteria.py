import os

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication import geomechanics_analysis as analysis

import test_helper


def run_simulation(test_directory, use_local_error_criteria=True, use_block_builder=None):
    """
    Runs a single stage, optionally combined with the local error criteria. No output files are
    written, such that the (reference) output of the test case is left untouched.
    """
    cwd = os.getcwd()
    os.chdir(test_directory)
    try:
        with open('ProjectParameters.json', 'r') as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())

        if parameters.Has("output_processes"):
            parameters.RemoveValue("output_processes")

        solver_settings = parameters["solver_settings"]
        if use_block_builder is not None:
            solver_settings["block_builder"].SetBool(use_block_builder)
        for name in ("use_local_error_criteria", "local_error_criteria_settings"):
            if solver_settings.Has(name):
                solver_settings.RemoveValue(name)
        solver_settings.AddBool("use_local_error_criteria", use_local_error_criteria)
        solver_settings.AddValue("local_error_criteria_settings", Kratos.Parameters('{"echo_level": 1}'))

        simulation = analysis.GeoMechanicsAnalysis(Kratos.Model(), parameters)
        simulation.Run()
    finally:
        os.chdir(cwd)

    return simulation


def get_computing_model_part(simulation):
    return simulation._GetSolver().GetComputingModelPart()


def get_on_integration_points(simulation, variable):
    model_part = get_computing_model_part(simulation)
    return [element.CalculateOnIntegrationPoints(variable, model_part.ProcessInfo)
            for element in model_part.Elements]


class KratosGeoMechanicsLocalErrorCriteriaTests(KratosUnittest.TestCase):
    """
    Checks that calculations in which the local error criteria are combined with the global criterion
    of a test case converge to the (regression) results of that test case.
    """

    def run_with_local_error_criteria(self, test_directory, **kwargs):
        simulation = run_simulation(test_directory, **kwargs)
        self.assertEqual(type(simulation._GetSolver().convergence_criterion).__name__, "AndCriteria")
        return simulation

    def assert_relatively_equal(self, value, expected_value, relative_tolerance=1.0e-3):
        self.assertAlmostEqual(value, expected_value, delta=relative_tolerance * abs(expected_value))

    def test_mohr_coulomb_regular_failure_zone(self):
        test_directory = test_helper.get_file_path(
            os.path.join('test_mohr_coulomb_with_tension_cutoff', 'test_dirichlet_regular_failure_zone_2d'))
        simulation = self.run_with_local_error_criteria(test_directory)

        stress_vector = get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)[0][0]
        self.assert_relatively_equal(stress_vector[0], 3.430712712091948)
        self.assert_relatively_equal(stress_vector[1], -25.759721409731483)
        plasticity_status = get_on_integration_points(simulation, KratosGeo.GEO_PLASTICITY_STATUS)[0][0]
        self.assertEqual(plasticity_status, 4)

    def test_mohr_coulomb_column_under_gravity(self):
        test_directory = test_helper.get_file_path(
            os.path.join('test_mohr_coulomb_with_tension_cutoff', 'test_column_under_gravity'))
        simulation = self.run_with_local_error_criteria(test_directory)

        stress_vector = get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)[0][0]
        self.assert_relatively_equal(stress_vector[0], -6300.238764425369)
        self.assert_relatively_equal(stress_vector[1], -22364.817908413854)

    def test_interface_coulomb(self):
        test_directory = test_helper.get_file_path(
            os.path.join('test_mohr_coulomb_with_tension_cutoff', 'interface_coulomb_2plus2'))
        # The internal force calculation of the Newton-Raphson strategy requires all degrees of
        # freedom to be part of the system, which is not the case for the elimination builder
        simulation = self.run_with_local_error_criteria(test_directory, use_block_builder=True)

        first_node = next(iter(get_computing_model_part(simulation).Nodes))
        self.assert_relatively_equal(first_node.GetSolutionStepValue(Kratos.REACTION)[0], 176.1686911267588)

    def test_mohr_coulomb_udsm(self):
        test_directory = test_helper.get_file_path(
            os.path.join('single_element_with_Mohr_Coulomb', 'UPwSmallStrainElement2D4N'))
        simulation = self.run_with_local_error_criteria(test_directory)

        for element_stress_vectors in get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR):
            for stress_vector in element_stress_vectors:
                self.assertVectorAlmostEqual(stress_vector, [1.5, -1.5, 0.0, 0.0], places=6)

    def test_mohr_coulomb_umat(self):
        test_directory = test_helper.get_file_path(
            os.path.join('Simple_Dike_Gravity_Loading', 'simple_dike_test_with_gravity_umat.gid'))
        reference = run_simulation(test_directory, use_local_error_criteria=False)
        simulation = self.run_with_local_error_criteria(test_directory)

        reference_stresses = get_on_integration_points(reference, Kratos.CAUCHY_STRESS_VECTOR)
        stresses = get_on_integration_points(simulation, Kratos.CAUCHY_STRESS_VECTOR)
        stress_scale = max(abs(value) for element_stresses in reference_stresses
                           for stress_vector in element_stresses for value in stress_vector)
        for element_reference_stresses, element_stresses in zip(reference_stresses, stresses):
            for reference_stress_vector, stress_vector in zip(element_reference_stresses, element_stresses):
                for reference_value, value in zip(reference_stress_vector, stress_vector):
                    self.assertAlmostEqual(value, reference_value, delta=1.0e-3 * stress_scale)


if __name__ == '__main__':
    KratosUnittest.main()

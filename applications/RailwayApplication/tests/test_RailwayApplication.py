
from KratosMultiphysics import *
from KratosMultiphysics.RailwayApplication import *

# Import Kratos "wrapper" for unittests
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Import the tests of test_classes to create the suits
from test_add_nodal_parameters import KratosRailwayAddNodalParametersTests
from test_call_uvec import KratosRailwayCallUvecTests
from test_fast_json_output_process import KratosRailwayFastJsonOutputProcessTests
from test_geomechanics_U_Pw_solver import KratosRailwayUPwSolverTests
from test_moving_load_multistage import KratosRailwayMovingLoadMultiStageTests
from test_set_parameter_field_process import KratosRailwaySetParameterFieldTests


def AssembleTestSuites():
    ''' Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    '''

    # Create an array with the selected tests (Small tests):
    small_test_cases = [
                        KratosRailwayAddNodalParametersTests,
                        KratosRailwayCallUvecTests,
                        KratosRailwayFastJsonOutputProcessTests,
                        KratosRailwayUPwSolverTests,
                        KratosRailwayMovingLoadMultiStageTests,
                        KratosRailwaySetParameterFieldTests
    ]

    night_test_cases = []
    night_test_cases.extend(small_test_cases)

    # Create an array with all long tests only for validations
    valid_test_cases = []

    # Create an array that contains all the tests from every testCase
    # in the list:

    all_test_cases = []
    all_test_cases.extend(night_test_cases)
    all_test_cases.extend(small_test_cases)
    all_test_cases.extend(valid_test_cases)
    suites = KratosUnittest.KratosSuites

    # add the tests to the corresponding suite,
    small_suite = suites['small']
    night_suite = suites['nightly']
    valid_suite = suites['validation']
    all_suite = suites['all']

    small_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(small_test_cases))
    night_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(night_test_cases))
    valid_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(valid_test_cases))
    all_suite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases(all_test_cases))

    return suites


if __name__ == '__main__':
    KratosMultiphysics.Logger.PrintInfo("Unittests", "\nRunning python tests ...")
    KratosUnittest.runTests(AssembleTestSuites())
    KratosMultiphysics.Logger.PrintInfo("Unittests", "Finished python tests!")

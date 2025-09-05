import os
import shutil
from pathlib import Path

import sys
sys.path.append(r"C:\software_development\KratosFork\bin\Release")
sys.path.append(r"C:\software_development\KratosFork\bin\Release\libs")

import KratosMultiphysics as Kratos
import KratosMultiphysics.RailwayApplication.geomechanics_analysis as analysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

from utils import Utils, assert_files_equal , assert_floats_in_files_almost_equal, RAILWAY_TEST_DIR

class KratosRailwayCallUvecTests(KratosUnittest.TestCase):

    def test_call_uvec(self):
        """
        Test the call of UVEC against benchmark
        """
        test_file_dir = RAILWAY_TEST_DIR / r"test_data/input_data_expected_moving_load_uvec"

        parameter_file_name = "ProjectParameters_stage1.json"

        Utils.run_multiple_stages(test_file_dir, [parameter_file_name])

        self.assertTrue(assert_files_equal(os.path.join(test_file_dir, "_output/porous_computational_model_part"),
                                           os.path.join(test_file_dir, "output/porous_computational_model_part")))

        shutil.rmtree(os.path.join(test_file_dir, "output"))
        os.remove(os.path.join(test_file_dir, "set_moving_load_process_moving_load_cloned_1.rest"))
        os.remove(os.path.join(test_file_dir, "set_moving_load_process_moving_load_cloned_2.rest"))

    def test_call_uvec_multi_stage(self):
        """
        Test the call of the UVEC in a multi-stage analysis
        """
        test_file_dir = RAILWAY_TEST_DIR / r"test_data/input_data_multi_stage_uvec"

        project_parameters = ["ProjectParameters_stage1.json", "ProjectParameters_stage2.json"]

        # run the analysis
        Utils.run_multiple_stages(test_file_dir, project_parameters)

        # calculated disp below first wheel
        calculated_disp_file = RAILWAY_TEST_DIR / r"test_data/input_data_multi_stage_uvec/output/calculated_disp"
        expected_disp_file = RAILWAY_TEST_DIR / r"test_data/input_data_multi_stage_uvec/_output/expected_disp"

        # check if calculated disp below first wheel is equal to expected disp
        are_files_equal, message = assert_floats_in_files_almost_equal(calculated_disp_file, expected_disp_file)

        # remove calculated disp file as data is appended
        calculated_disp_file.unlink()
        self.assertTrue(are_files_equal)

        expected_vtk_output_dir = RAILWAY_TEST_DIR / "test_data/input_data_multi_stage_uvec/_output/all"

        main_vtk_output_dir = (RAILWAY_TEST_DIR /
                               "test_data/input_data_multi_stage_uvec/output/porous_computational_model_part_1")
        stage_vtk_output_dir = (RAILWAY_TEST_DIR /
                                "test_data/input_data_multi_stage_uvec/output/porous_computational_model_part_2")

        # move all vtk files in stage vtk output dir to main vtk output dir
        for file in os.listdir(stage_vtk_output_dir):
            if file.endswith(".vtk"):
                os.rename(stage_vtk_output_dir / file, main_vtk_output_dir / file)

        # remove the stage vtk output dir if it is empty
        if not os.listdir(stage_vtk_output_dir):
            os.rmdir(stage_vtk_output_dir)

        # check if vtk files are equal
        self.assertTrue(assert_files_equal(expected_vtk_output_dir, main_vtk_output_dir))
        shutil.rmtree(main_vtk_output_dir)
        os.remove(RAILWAY_TEST_DIR / "test_data/input_data_multi_stage_uvec/set_moving_load_process_moving_load_cloned_1.rest")
        os.remove(RAILWAY_TEST_DIR / "test_data/input_data_multi_stage_uvec/set_moving_load_process_moving_load_cloned_2.rest")

    def test_call_uvec_multi_stage_expected_fail(self):
        """
        Test the call of the UVEC in a multi-stage analysis. This test is expected to fail, as extra DOFS are added in the
        second stage
        """
        test_file_dir = RAILWAY_TEST_DIR / r"test_data/input_data_multi_stage_uvec"

        project_parameters = ["ProjectParameters_stage1.json", "ProjectParameters_stage2.json"]

        cwd = os.getcwd()

        # initialize model
        model = Kratos.Model()

        # read first stage parameters
        os.chdir(test_file_dir)
        with open(project_parameters[0], 'r') as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())

        # remove rotation dofs from first stage
        parameters["solver_settings"]["rotation_dofs"].SetBool(False)

        # create and run first stage
        stage = analysis.StemGeoMechanicsAnalysis(model, parameters)
        stage.Run()

        # read second stage parameters
        with open(project_parameters[1], 'r') as parameter_file:
            parameters = Kratos.Parameters(parameter_file.read())

        # make sure rotation dofs are added in second stage
        parameters["solver_settings"]["rotation_dofs"].SetBool(True)

        # create second stage with expected fail
        with self.assertRaises(RuntimeError) as excinfo:
            analysis.StemGeoMechanicsAnalysis(model, parameters)

        # check if the error message is as expected
        self.assertTrue('Error: Attempting to add the variable "ROTATION" to the model part with name "PorousDomain"'
                in str(excinfo.exception))

        # change working directory back to original working directory
        os.chdir(cwd)

        # remove uvec disp file
        calculated_disp_file = RAILWAY_TEST_DIR / r"test_data/input_data_multi_stage_uvec/output/calculated_disp"
        calculated_disp_file.unlink()
        shutil.rmtree(os.path.join(test_file_dir, "output/porous_computational_model_part_1"))
        os.remove(RAILWAY_TEST_DIR / "test_data/input_data_multi_stage_uvec/set_moving_load_process_moving_load_cloned_1.rest")
        os.remove(RAILWAY_TEST_DIR / "test_data/input_data_multi_stage_uvec/set_moving_load_process_moving_load_cloned_2.rest")

if __name__ == '__main__':
    KratosUnittest.main()

import json
import math
import os
import tempfile

import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper


class _Mesh:
    """Minimal in-memory representation of a Kratos mdpa file"""

    def __init__(self):
        self.nodes = {}
        self.elements = {}
        self.conditions = {}
        self.sub_model_parts = {}

    def add_node(self, node_id, x, y, z=0.0):
        self.nodes[node_id] = (x, y, z)
        return node_id

    def add_element(self, element_name, element_id, properties_id, node_ids, sub_model_part):
        self.elements.setdefault(element_name, []).append((element_id, properties_id, list(node_ids)))
        self._add_to_sub_model_part(sub_model_part, node_ids, elements=[element_id])

    def add_condition(self, condition_name, condition_id, node_ids, sub_model_part):
        self.conditions.setdefault(condition_name, []).append((condition_id, 0, list(node_ids)))
        self._add_to_sub_model_part(sub_model_part, node_ids, conditions=[condition_id])

    def add_nodes_to_sub_model_part(self, sub_model_part, node_ids):
        self._add_to_sub_model_part(sub_model_part, node_ids)

    def _add_to_sub_model_part(self, name, node_ids, elements=(), conditions=()):
        part = self.sub_model_parts.setdefault(name, {"nodes": set(), "elements": set(), "conditions": set()})
        part["nodes"].update(node_ids)
        part["elements"].update(elements)
        part["conditions"].update(conditions)

    def renumber_nodes(self, key):
        """Renumbers the nodes in the order given by key(coordinates)"""
        mapping = {old: new for new, old in enumerate(sorted(self.nodes, key=lambda n: key(self.nodes[n])), start=1)}
        self.nodes = {mapping[n]: coordinates for n, coordinates in self.nodes.items()}
        for entities in (self.elements, self.conditions):
            for name, entity_list in entities.items():
                entities[name] = [(i, p, [mapping[n] for n in node_ids]) for i, p, node_ids in entity_list]
        for part in self.sub_model_parts.values():
            part["nodes"] = {mapping[n] for n in part["nodes"]}

    def find_node(self, x, y, z=0.0, tolerance=1.0e-9):
        for node_id, coordinates in self.nodes.items():
            if math.dist(coordinates, (x, y, z)) < tolerance:
                return node_id
        return None

    def write(self, file_path):
        properties_ids = {0} | {properties_id for entities in self.elements.values() for _, properties_id, _ in entities}
        with open(file_path, "w") as mdpa:
            for properties_id in sorted(properties_ids):
                mdpa.write(f"Begin Properties {properties_id}\nEnd Properties\n\n")

            mdpa.write("Begin Nodes\n")
            for node_id, (x, y, z) in sorted(self.nodes.items()):
                mdpa.write(f"  {node_id} {x:.12f} {y:.12f} {z:.12f}\n")
            mdpa.write("End Nodes\n\n")

            for block_name, entities in (("Elements", self.elements), ("Conditions", self.conditions)):
                for entity_name, entity_list in entities.items():
                    mdpa.write(f"Begin {block_name} {entity_name}\n")
                    for entity_id, properties_id, node_ids in entity_list:
                        mdpa.write(f"  {entity_id} {properties_id} {' '.join(str(n) for n in node_ids)}\n")
                    mdpa.write(f"End {block_name}\n\n")

            for name, part in self.sub_model_parts.items():
                mdpa.write(f"Begin SubModelPart {name}\n")
                for block_name, key in (("Nodes", "nodes"), ("Elements", "elements"), ("Conditions", "conditions")):
                    mdpa.write(f"  Begin SubModelPart{block_name}\n")
                    for entity_id in sorted(part[key]):
                        mdpa.write(f"    {entity_id}\n")
                    mdpa.write(f"  End SubModelPart{block_name}\n")
                mdpa.write("End SubModelPart\n\n")


def _create_structured_quadrilateral_mesh(mesh, x_coordinates, y_coordinates, sub_model_part, x_shift=None):
    """Creates UPwSmallStrainElement2D4N elements on a structured grid. The optional function
    x_shift(x, y) distorts the grid in horizontal direction."""
    node_ids = {}
    for j, y in enumerate(y_coordinates):
        for i, x in enumerate(x_coordinates):
            shifted_x = x + (x_shift(x, y) if x_shift else 0.0)
            node_ids[(i, j)] = mesh.add_node(len(mesh.nodes) + 1, shifted_x, y)

    element_id = 1
    for j in range(len(y_coordinates) - 1):
        for i in range(len(x_coordinates) - 1):
            nodes = [node_ids[(i, j)], node_ids[(i + 1, j)], node_ids[(i + 1, j + 1)], node_ids[(i, j + 1)]]
            mesh.add_element("UPwSmallStrainElement2D4N", element_id, 1, nodes, sub_model_part)
            element_id += 1
    return node_ids


def _create_structured_hexahedral_mesh(mesh, x_coordinates, y_coordinates, z_coordinates, sub_model_part):
    node_ids = {}
    for k, z in enumerate(z_coordinates):
        for j, y in enumerate(y_coordinates):
            for i, x in enumerate(x_coordinates):
                node_ids[(i, j, k)] = mesh.add_node(len(mesh.nodes) + 1, x, y, z)

    element_id = 1
    for k in range(len(z_coordinates) - 1):
        for j in range(len(y_coordinates) - 1):
            for i in range(len(x_coordinates) - 1):
                nodes = [node_ids[(i, j, k)], node_ids[(i + 1, j, k)], node_ids[(i + 1, j + 1, k)],
                         node_ids[(i, j + 1, k)], node_ids[(i, j, k + 1)], node_ids[(i + 1, j, k + 1)],
                         node_ids[(i + 1, j + 1, k + 1)], node_ids[(i, j + 1, k + 1)]]
                mesh.add_element("UPwSmallStrainElement3D8N", element_id, 1, nodes, sub_model_part)
                element_id += 1
    return node_ids


def _linspace(start, end, number_of_divisions):
    return [start + (end - start) * i / number_of_divisions for i in range(number_of_divisions + 1)]


def _soil_material(model_part_name, domain_size, young_modulus, poisson_ratio):
    variables = {
        "IGNORE_UNDRAINED": True,
        "YOUNG_MODULUS": young_modulus,
        "POISSON_RATIO": poisson_ratio,
        "DENSITY_SOLID": 0.0,
        "DENSITY_WATER": 1.0,
        "POROSITY": 0.3,
        "BULK_MODULUS_SOLID": 1.0e16,
        "BULK_MODULUS_FLUID": 2.0e-30,
        "PERMEABILITY_XX": 4.5e-30,
        "PERMEABILITY_YY": 4.5e-30,
        "PERMEABILITY_XY": 0.0,
        "DYNAMIC_VISCOSITY": 1.0e-3,
    }
    if domain_size == 2:
        law = "GeoLinearElasticPlaneStrain2DLaw"
        variables["THICKNESS"] = 1.0
    else:
        law = "GeoIncrementalLinearElastic3DLaw"
        variables.update({"PERMEABILITY_ZZ": 4.5e-30, "PERMEABILITY_YZ": 0.0, "PERMEABILITY_ZX": 0.0})
    return {"model_part_name": f"PorousDomain.{model_part_name}", "properties_id": 1,
            "Material": {"constitutive_law": {"name": law}, "Variables": variables, "Tables": {}}}


def _beam_material(model_part_name, properties_id, domain_size, variables):
    law = "TimoshenkoBeamElasticConstitutiveLaw" if domain_size == 2 else "TimoshenkoBeamElasticConstitutiveLaw3D"
    return {"model_part_name": f"PorousDomain.{model_part_name}", "properties_id": properties_id,
            "Material": {"constitutive_law": {"name": f"KratosMultiphysics.StructuralMechanicsApplication.{law}"},
                         "Variables": variables, "Tables": {}}}


def _vector_process(model_part_name, variable_name, active, value, is_fixed=None):
    parameters = {"model_part_name": f"PorousDomain.{model_part_name}", "variable_name": variable_name,
                  "active": active, "value": value, "table": [0, 0, 0]}
    if is_fixed is not None:
        parameters["is_fixed"] = is_fixed
    return {"python_module": "apply_vector_constraint_table_process",
            "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
            "process_name": "ApplyVectorConstraintTableProcess", "Parameters": parameters}


def _fixed_water_pressure_process(model_part_name):
    return {"python_module": "apply_scalar_constraint_table_process",
            "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
            "process_name": "ApplyScalarConstraintTableProcess",
            "Parameters": {"model_part_name": f"PorousDomain.{model_part_name}", "variable_name": "WATER_PRESSURE",
                           "active": True, "is_fixed": True, "value": 0.0, "table": 0}}


def _find_hosts_process(embedded_model_part_names, host_model_part_name):
    return {"python_module": "find_embedded_beam_host_elements_process",
            "kratos_module": "KratosMultiphysics.GeoMechanicsApplication",
            "process_name": "FindEmbeddedBeamHostElementsProcess",
            "Parameters": {"model_part_name_list": [f"PorousDomain.{name}" for name in embedded_model_part_names],
                           "host_model_part_name": f"PorousDomain.{host_model_part_name}"}}


def _project_parameters(domain_size, domain_parts, processes):
    return {
        "problem_data": {"problem_name": "embedded_beam", "start_time": 0.0, "end_time": 1.0, "echo_level": 0,
                         "parallel_type": "OpenMP", "number_of_threads": 1},
        "solver_settings": {
            "solver_type": "U_Pw", "model_part_name": "PorousDomain", "domain_size": domain_size,
            "model_import_settings": {"input_type": "mdpa", "input_filename": "embedded_beam"},
            "material_import_settings": {"materials_filename": "MaterialParameters.json"},
            "time_stepping": {"time_step": 1.0},
            "buffer_size": 2, "echo_level": 0, "clear_storage": False, "compute_reactions": True,
            "move_mesh_flag": False, "reform_dofs_at_each_step": False, "block_builder": True,
            "solution_type": "Quasi-Static", "scheme_type": "Newmark", "reset_displacements": True,
            "newmark_beta": 0.25, "newmark_gamma": 0.5, "newmark_theta": 0.5, "rayleigh_m": 0.0, "rayleigh_k": 0.0,
            "strategy_type": "newton_raphson", "convergence_criterion": "displacement_criterion",
            "displacement_relative_tolerance": 1.0e-9, "displacement_absolute_tolerance": 1.0e-12,
            "min_iterations": 6, "max_iterations": 30, "number_cycles": 1, "reduction_factor": 1.0,
            "quasi_newton_type": "broyden", "relaxation_factor": 1.0,
            "increase_factor": 1.0, "calculate_reactions": True, "rotation_dofs": True,
            "linear_solver_settings": {"solver_type": "LinearSolversApplication.sparse_lu", "scaling": True},
            "problem_domain_sub_model_part_list": domain_parts,
            "body_domain_sub_model_part_list": domain_parts
        },
        "output_processes": {},
        "processes": {"auxiliary_process_list": processes}
    }


def _run(mesh, materials, project_parameters):
    with tempfile.TemporaryDirectory() as directory:
        mesh.write(os.path.join(directory, "embedded_beam.mdpa"))
        with open(os.path.join(directory, "MaterialParameters.json"), "w") as f:
            json.dump({"properties": materials}, f, indent=2)
        with open(os.path.join(directory, "ProjectParameters.json"), "w") as f:
            json.dump(project_parameters, f, indent=2)
        working_directory = os.getcwd()
        try:
            return test_helper.run_kratos(directory)
        finally:
            os.chdir(working_directory)


def _get_displacement(simulation, node_id):
    return simulation.model.GetModelPart("PorousDomain").GetNode(node_id).GetSolutionStepValue(Kratos.DISPLACEMENT)


def _get_beam_moments(simulation, model_part_name, variable):
    model_part = simulation.model.GetModelPart(f"PorousDomain.{model_part_name}")
    process_info = model_part.ProcessInfo
    result = []
    for element in model_part.Elements:
        result.extend(element.CalculateOnIntegrationPoints(variable, process_info))
    return result


class KratosGeoMechanicsEmbeddedBeamTests(KratosUnittest.TestCase):
    """
    Tests of the embedded beam element according to Sadek & Shahrour (2004), "A three dimensional
    embedded beam element for reinforced geomaterials", Int. J. Numer. Anal. Meth. Geomech. 28:931-946.
    """

    # Pile of the verification example in section 3.1 of the paper: EA = 1155.2 MN, EI = 892.5 kNm2
    PILE_YOUNG_MODULUS = 2.1e8
    PILE_VARIABLES_2D = {"YOUNG_MODULUS": 2.1e8, "POISSON_RATIO": 0.3, "DENSITY": 0.0,
                         "CROSS_AREA": 1155.2e3 / 2.1e8, "I33": 892.5 / 2.1e8,
                         "AREA_EFFECTIVE_Y": 0.5 * 1155.2e3 / 2.1e8}

    def _create_2d_pile_problem(self, use_embedded_beam, pile_nodes_on_soil_nodes, inclination, element_size=0.5):
        """A pile embedded in a 20 m wide and 10 m deep soil layer, loaded by a horizontal force of
        50 kN at its head. The vertical extent of the pile is 5 m."""
        width, depth = 20.0, 10.0
        pile_head_x, pile_depth = 10.0, 5.0
        tan_inclination = math.tan(math.radians(inclination))

        def conforming_shift(x, y):
            # distort the grid such that the column of nodes at the pile head follows the pile
            weight = x / pile_head_x if x <= pile_head_x else (width - x) / (width - pile_head_x)
            return weight * tan_inclination * min(-y, pile_depth)

        mesh = _Mesh()
        _create_structured_quadrilateral_mesh(mesh, _linspace(0.0, width, round(width / element_size)),
                                              _linspace(-depth, 0.0, round(depth / element_size)), "soil",
                                              conforming_shift if pile_nodes_on_soil_nodes else None)

        number_of_pile_elements = round(pile_depth / element_size)
        pile_node_ids = []
        for k in range(number_of_pile_elements + 1):
            y = -k * element_size
            x = pile_head_x + tan_inclination * (-y)
            node_id = mesh.find_node(x, y)
            if node_id is None or (k > 0 and not pile_nodes_on_soil_nodes):
                node_id = mesh.add_node(len(mesh.nodes) + 1, x, y)
            pile_node_ids.append(node_id)

        element_name = "GeoEmbeddedBeamElement2D2N" if use_embedded_beam else "LinearTimoshenkoBeamElement2D2N"
        number_of_soil_elements = sum(len(entities) for entities in mesh.elements.values())
        for k in range(number_of_pile_elements):
            mesh.add_element(element_name, number_of_soil_elements + k + 1, 2, pile_node_ids[k:k + 2], "pile")

        head_node_id = pile_node_ids[0]
        mesh.add_condition("UPwForceCondition2D1N", 1, [head_node_id], "load")
        mesh.add_nodes_to_sub_model_part("bottom", [n for n, (x, y, z) in mesh.nodes.items() if abs(y + depth) < 1e-9])
        mesh.add_nodes_to_sub_model_part("sides", [n for n, (x, y, z) in mesh.nodes.items()
                                                   if abs(x) < 1e-9 or abs(x - width) < 1e-9])

        materials = [_soil_material("soil", 2, 1.0e5, 0.33), _beam_material("pile", 2, 2, self.PILE_VARIABLES_2D)]
        processes = [_vector_process("bottom", "DISPLACEMENT", [True, True, False], [0.0, 0.0, 0.0], [True, True, False]),
                     _vector_process("sides", "DISPLACEMENT", [True, False, False], [0.0, 0.0, 0.0], [True, False, False]),
                     _fixed_water_pressure_process("soil"),
                     _vector_process("load", "POINT_LOAD", [True, False, False], [50.0, 0.0, 0.0])]
        if use_embedded_beam:
            processes.append(_find_hosts_process(["pile"], "soil"))

        simulation = _run(mesh, materials, _project_parameters(2, ["soil", "pile"], processes))
        return simulation, head_node_id

    def test_embedded_beam_on_soil_nodes_reproduces_conventional_beam_2d(self):
        """When the nodes of the embedded beam coincide with soil nodes, the embedded beam must give
        exactly the same results as the conventional beam element connected to these nodes."""
        conventional, head_node_id = self._create_2d_pile_problem(False, True, 0.0)
        embedded, _ = self._create_2d_pile_problem(True, True, 0.0)

        conventional_nodes = conventional.model.GetModelPart("PorousDomain.soil").Nodes
        maximum_displacement = max(abs(node.GetSolutionStepValue(Kratos.DISPLACEMENT_X)) for node in conventional_nodes)
        self.assertGreater(maximum_displacement, 0.0)
        for node in conventional_nodes:
            expected = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            actual = _get_displacement(embedded, node.Id)
            for i in range(2):
                self.assertAlmostEqual(actual[i], expected[i], delta=1.0e-8 * maximum_displacement)

        expected_moments = _get_beam_moments(conventional, "pile", KratosStructural.BENDING_MOMENT)
        actual_moments = _get_beam_moments(embedded, "pile", KratosStructural.BENDING_MOMENT)
        maximum_moment = max(abs(moment) for moment in expected_moments)
        for actual, expected in zip(actual_moments, expected_moments):
            self.assertAlmostEqual(actual, expected, delta=1.0e-6 * maximum_moment)

    def test_inclined_embedded_beam_agrees_with_conventional_beam_2d(self):
        """Equivalent of the verification in section 3.1 of Sadek & Shahrour (2004): a pile inclined
        by 30 degrees, loaded by a horizontal force at its head. The conventional beam requires a
        mesh that follows the pile, the embedded beam is placed in a regular mesh. Both solutions
        converge to each other upon mesh refinement (the difference of the head displacement is
        8.9%, 6.0% and 2.7% for element sizes of 0.5, 0.25 and 0.125 m respectively)."""
        element_size = 0.25
        conventional, _ = self._create_2d_pile_problem(False, True, 30.0, element_size)
        embedded, _ = self._create_2d_pile_problem(True, False, 30.0, element_size)

        def lateral_displacement_profile(simulation):
            nodes = simulation.model.GetModelPart("PorousDomain.pile").Nodes
            return [node.GetSolutionStepValue(Kratos.DISPLACEMENT_X) for node in sorted(nodes, key=lambda n: -n.Y)]

        expected_profile = lateral_displacement_profile(conventional)
        actual_profile = lateral_displacement_profile(embedded)
        self.assertEqual(len(actual_profile), round(5.0 / element_size) + 1)
        self.assertEqual(len(actual_profile), len(expected_profile))

        head_displacement = expected_profile[0]
        self.assertGreater(head_displacement, 0.0)
        for actual, expected in zip(actual_profile, expected_profile):
            self.assertAlmostEqual(actual, expected, delta=0.1 * head_displacement)

    def _run_2d_dynamic_pile_problem(self, use_embedded_beam):
        """A vertical pile in a 10 m wide and 10 m deep soil layer, loaded by a horizontal step load
        at its head (dynamic analysis)."""
        width, depth, pile_depth = 10.0, 10.0, 6.0
        mesh = _Mesh()
        _create_structured_quadrilateral_mesh(mesh, _linspace(0.0, width, 10), _linspace(-depth, 0.0, 10), "soil")

        # Number the nodes from the top, such that the first degrees of freedom of the system are not fixed
        mesh.renumber_nodes(lambda coordinates: (-coordinates[1], coordinates[0]))

        # The embedded beam has its own nodes, which coincide with soil nodes
        pile_node_ids = []
        for k in range(round(pile_depth) + 1):
            soil_node_id = mesh.find_node(width / 2, -float(k))
            pile_node_ids.append(mesh.add_node(len(mesh.nodes) + 1, width / 2, -float(k))
                                 if use_embedded_beam else soil_node_id)
        element_name = "GeoEmbeddedBeamElement2D2N" if use_embedded_beam else "LinearTimoshenkoBeamElement2D2N"
        number_of_soil_elements = sum(len(entities) for entities in mesh.elements.values())
        for k in range(round(pile_depth)):
            mesh.add_element(element_name, number_of_soil_elements + k + 1, 2, pile_node_ids[k:k + 2], "pile")

        mesh.add_condition("UPwForceCondition2D1N", 1, [mesh.find_node(width / 2, 0.0)], "load")
        mesh.add_nodes_to_sub_model_part("bottom", [n for n, (x, y, z) in mesh.nodes.items() if abs(y + depth) < 1e-9])

        soil = _soil_material("soil", 2, 1.0e5, 0.3)
        soil["Material"]["Variables"]["DENSITY_SOLID"] = 2.0
        materials = [soil, _beam_material("pile", 2, 2, self.PILE_VARIABLES_2D)]
        processes = [_vector_process("bottom", "DISPLACEMENT", [True, True, False], [0.0, 0.0, 0.0], [True, True, False]),
                     _fixed_water_pressure_process("soil"),
                     _vector_process("load", "POINT_LOAD", [True, False, False], [10.0, 0.0, 0.0])]
        if use_embedded_beam:
            processes.append(_find_hosts_process(["pile"], "soil"))

        project_parameters = _project_parameters(2, ["soil", "pile"], processes)
        project_parameters["problem_data"]["end_time"] = 0.05
        project_parameters["solver_settings"].update({"solution_type": "Dynamic", "time_stepping": {"time_step": 0.005}})
        return _run(mesh, materials, project_parameters)

    def test_embedded_beam_reproduces_conventional_beam_in_dynamic_analysis_2d(self):
        """The nodes of the embedded beam coincide with (but are not the same as) soil nodes, such that
        the embedded beam must give the same results as a conventional beam connected to the soil
        nodes. The translations of the embedded beam nodes are not part of the system of equations."""
        conventional = self._run_2d_dynamic_pile_problem(False)
        embedded = self._run_2d_dynamic_pile_problem(True)

        conventional_nodes = conventional.model.GetModelPart("PorousDomain.soil").Nodes
        maximum_displacement = max(abs(node.GetSolutionStepValue(Kratos.DISPLACEMENT_X)) for node in conventional_nodes)
        self.assertGreater(maximum_displacement, 0.0)
        expected_displacements = {}
        for node in conventional_nodes:
            expected = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            expected_displacements[(round(node.X, 9), round(node.Y, 9))] = expected
            actual = _get_displacement(embedded, node.Id)
            for i in range(2):
                self.assertAlmostEqual(actual[i], expected[i], delta=1.0e-8 * maximum_displacement)

        # The displacements of the embedded beam nodes are interpolated from the soil
        for node in embedded.model.GetModelPart("PorousDomain.pile").Nodes:
            expected = expected_displacements[(round(node.X, 9), round(node.Y, 9))]
            actual = node.GetSolutionStepValue(Kratos.DISPLACEMENT)
            for i in range(2):
                self.assertAlmostEqual(actual[i], expected[i], delta=1.0e-8 * maximum_displacement)

    def test_embedded_pile_in_very_soft_soil_behaves_as_cantilever_3d(self):
        """A vertical pile, clamped at its base, embedded in a very soft soil (not coinciding with
        soil nodes) and loaded by a horizontal force at its head, deflects as a cantilever beam."""
        size, depth = 4.0, 10.0
        mesh = _Mesh()
        _create_structured_hexahedral_mesh(mesh, _linspace(0.0, size, 4), _linspace(0.0, size, 4),
                                           _linspace(-depth, 0.0, 10), "soil")

        pile_x = pile_y = 2.5
        pile_node_ids = [mesh.add_node(len(mesh.nodes) + 1, pile_x, pile_y, -float(k)) for k in range(11)]
        number_of_soil_elements = sum(len(entities) for entities in mesh.elements.values())
        for k in range(10):
            mesh.add_element("GeoEmbeddedBeamElement3D2N", number_of_soil_elements + k + 1, 2,
                             pile_node_ids[k:k + 2], "pile")
        mesh.add_nodes_to_sub_model_part("pile_base", [pile_node_ids[-1]])

        # The load on the pile head is transferred to the soil nodes around it
        force = 10.0
        for condition_id, (x, y) in enumerate([(2.0, 2.0), (3.0, 2.0), (3.0, 3.0), (2.0, 3.0)], start=1):
            mesh.add_condition("UPwForceCondition3D1N", condition_id, [mesh.find_node(x, y, 0.0)], "load")
        mesh.add_nodes_to_sub_model_part("bottom", [n for n, (x, y, z) in mesh.nodes.items()
                                                    if abs(z + depth) < 1e-9 and n not in pile_node_ids])

        diameter = 0.5
        area = math.pi * diameter ** 2 / 4.0
        second_moment_of_area = math.pi * diameter ** 4 / 64.0
        young_modulus, poisson_ratio = 3.0e7, 0.2
        pile_variables = {"YOUNG_MODULUS": young_modulus, "POISSON_RATIO": poisson_ratio, "DENSITY": 0.0,
                          "CROSS_AREA": area, "I22": second_moment_of_area, "I33": second_moment_of_area,
                          "AREA_EFFECTIVE_Y": 0.9 * area, "AREA_EFFECTIVE_Z": 0.9 * area}
        materials = [_soil_material("soil", 3, 1.0e-2, 0.0), _beam_material("pile", 2, 3, pile_variables)]
        processes = [_vector_process("bottom", "DISPLACEMENT", [True, True, True], [0.0, 0.0, 0.0], [True, True, True]),
                     _vector_process("pile_base", "ROTATION", [True, True, True], [0.0, 0.0, 0.0], [True, True, True]),
                     _fixed_water_pressure_process("soil"),
                     _vector_process("load", "POINT_LOAD", [True, False, False], [force / 4.0, 0.0, 0.0]),
                     _find_hosts_process(["pile"], "soil")]

        simulation = _run(mesh, materials, _project_parameters(3, ["soil", "pile"], processes))

        shear_modulus = young_modulus / (2.0 * (1.0 + poisson_ratio))
        expected_head_displacement = (force * depth ** 3 / (3.0 * young_modulus * second_moment_of_area) +
                                      force * depth / (shear_modulus * 0.9 * area))
        head_displacement = _get_displacement(simulation, pile_node_ids[0])
        self.assertAlmostEqual(head_displacement[0], expected_head_displacement, delta=0.01 * expected_head_displacement)
        self.assertAlmostEqual(head_displacement[1], 0.0, delta=1.0e-6 * expected_head_displacement)

        # The bending moment increases linearly with the distance to the head: M = F |z|
        model_part = simulation.model.GetModelPart("PorousDomain.pile")
        for element in model_part.Elements:
            geometry = element.GetGeometry()
            z_start, z_end = geometry[0].Z, geometry[1].Z
            moments = element.CalculateOnIntegrationPoints(KratosStructural.BENDING_MOMENT_Y, model_part.ProcessInfo)
            gauss_point_coordinates = [-math.sqrt(0.6), 0.0, math.sqrt(0.6)]
            self.assertEqual(len(moments), len(gauss_point_coordinates))
            for xi, moment in zip(gauss_point_coordinates, moments):
                z = 0.5 * (1.0 - xi) * z_start + 0.5 * (1.0 + xi) * z_end
                self.assertAlmostEqual(abs(moment), force * abs(z), delta=0.01 * force * depth)

    def _run_2d_pull_out(self, prescribed_displacement):
        """A pile (with its own degrees of freedom) in a fixed soil, connected to the soil by an
        embedded beam interface with a cohesion of 10 kPa, pulled out at its head."""
        mesh = _Mesh()
        _create_structured_quadrilateral_mesh(mesh, _linspace(0.0, 4.0, 4), _linspace(-6.0, 0.0, 6), "soil")
        mesh.add_nodes_to_sub_model_part("soil_nodes", list(mesh.nodes.keys()))

        pile_node_ids = [mesh.add_node(len(mesh.nodes) + 1, 2.5, -float(k)) for k in range(6)]
        number_of_soil_elements = sum(len(entities) for entities in mesh.elements.values())
        for k in range(5):
            mesh.add_element("LinearTimoshenkoBeamElement2D2N", number_of_soil_elements + k + 1, 2,
                             pile_node_ids[k:k + 2], "pile")
            mesh.add_element("GeoEmbeddedBeamInterfaceElement2D2N", number_of_soil_elements + k + 11, 3,
                             pile_node_ids[k:k + 2], "pile_interface")
        mesh.add_nodes_to_sub_model_part("pile_head", [pile_node_ids[0]])

        # a practically rigid pile
        pile_variables = {"YOUNG_MODULUS": 2.1e8, "POISSON_RATIO": 0.3, "DENSITY": 0.0, "CROSS_AREA": 100.0,
                          "I33": 1.0, "AREA_EFFECTIVE_Y": 50.0}
        interface = {"model_part_name": "PorousDomain.pile_interface", "properties_id": 3,
                     "Material": {"Variables": {"INTERFACE_SHEAR_STIFFNESS": 1.0e5, "INTERFACE_NORMAL_STIFFNESS": 1.0e5,
                                                "GEO_COHESION": 10.0}, "Tables": {}}}
        materials = [_soil_material("soil", 2, 1.0e5, 0.3), _beam_material("pile", 2, 2, pile_variables), interface]
        processes = [_vector_process("soil_nodes", "DISPLACEMENT", [True, True, False], [0.0, 0.0, 0.0], [True, True, False]),
                     _vector_process("pile_head", "DISPLACEMENT", [False, True, False],
                                     [0.0, prescribed_displacement, 0.0], [False, True, False]),
                     _fixed_water_pressure_process("soil"),
                     _find_hosts_process(["pile_interface"], "soil")]

        simulation = _run(mesh, materials, _project_parameters(2, ["soil", "pile", "pile_interface"], processes))
        model_part = simulation.model.GetModelPart("PorousDomain")
        return model_part.GetNode(pile_node_ids[0]).GetSolutionStepValue(Kratos.REACTION_Y)

    def test_embedded_beam_interface_pull_out_2d(self):
        """In the elastic range, the pull-out force equals K_t * u * L (the pile is practically
        rigid). Once the pile slides along its full length, the pull-out force equals the shear
        strength of the interface: 2 faces * cohesion * L."""
        pile_length, shear_stiffness, cohesion = 5.0, 1.0e5, 10.0

        elastic_reaction = self._run_2d_pull_out(1.0e-5)
        self.assertAlmostEqual(abs(elastic_reaction), shear_stiffness * 1.0e-5 * pile_length, delta=1.0e-3)

        plastic_reaction = self._run_2d_pull_out(1.0e-2)
        self.assertAlmostEqual(abs(plastic_reaction), 2.0 * cohesion * pile_length, delta=1.0e-6)


if __name__ == '__main__':
    KratosUnittest.main()

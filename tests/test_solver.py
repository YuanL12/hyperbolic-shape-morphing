import math
import unittest

from hyper_morph import HarmonicMapSolver
from hyper_morph.solver import hyperbolic_distance, log_map


class EdgeForceModelTest(unittest.TestCase):
    def test_hyperbolic_mean_value_scales_log_map_by_sinh_d_over_d(self) -> None:
        vertices = [[0.0, 0.0], [0.3, 0.1]]
        data = {
            "vertices": vertices,
            "edges": [[0, 1]],
            "directed_edge_weights": [[2.0, 3.0]],
            "edge_force_model": "hyperbolic_mean_value",
        }

        solver = HarmonicMapSolver(data)
        positions = [complex(*vertex) for vertex in vertices]
        _, gradient = solver.energy_and_raw_gradient(positions)

        distance = hyperbolic_distance(*positions)
        scale = math.sinh(distance) / distance
        expected = -2.0 * scale * log_map(*positions)
        self.assertAlmostEqual(gradient[0].real, expected.real)
        self.assertAlmostEqual(gradient[0].imag, expected.imag)


class StopReasonTest(unittest.TestCase):
    def test_result_records_stop_reason_without_position_history(self) -> None:
        data = {
            "vertices": [[0.0, 0.0], [0.3, 0.1]],
            "edges": [[0, 1]],
            "iterations": 1,
            "tolerance": 0.0,
            "line_search_objective": "none",
            "convergence_criterion": "relative_step",
        }

        result = HarmonicMapSolver(data).solve()

        self.assertEqual(result["stop_reason"], "max_iterations")
        self.assertEqual(result["iterations_recorded"], 1)
        self.assertNotIn("position_history", result)
        self.assertNotIn("vertices_history", result)

    def test_zero_iteration_run_has_explicit_stop_reason(self) -> None:
        result = HarmonicMapSolver(
            {
                "vertices": [[0.0, 0.0], [0.3, 0.1]],
                "edges": [[0, 1]],
                "iterations": 0,
            }
        ).solve()

        self.assertEqual(result["stop_reason"], "zero_iterations")
        self.assertEqual(result["iterations_recorded"], 0)


if __name__ == "__main__":
    unittest.main()

#
# Test for the operator class
#
import pybamm
import tests.shared as shared

import numpy as np
import unittest


class TestFiniteElementDiscretisation(unittest.TestCase):
    def test_grad_convergence_without_bcs(self):
        # Convergence
        whole_cell = ["negative electrode", "separator", "positive electrode"]
        var = pybamm.Variable("var", domain=whole_cell)
        grad_eqn = pybamm.grad(var)

        # Function for convergence testing
        def get_l2_error(n):
            # Set up discretisation
            n = 3 * round(n / 3)
            # create discretisation
            defaults = shared.TestDefaults1DMacro(submesh_type=pybamm.ScikitFemMesh1D)
            defaults.set_equal_pts(n)
            disc = pybamm.ScikitFiniteElementDiscretisation(
                defaults.mesh_type, defaults.submesh_pts, defaults.submesh_types
            )
            disc.mesh_geometry(defaults.geometry)
            mesh = disc.mesh

            combined_submesh = mesh.combine_submeshes(*whole_cell)
            # Define exact solutions
            y = np.sin(combined_submesh.nodes)
            grad_exact = np.cos(combined_submesh.edges[1:-1])

            # Discretise and evaluate
            y_slices = disc.get_variable_slices([var])
            grad_eqn_disc = disc.process_symbol(grad_eqn, y_slices, {})
            grad_approx = grad_eqn_disc.evaluate(None, y)

            # Calculate errors
            return np.linalg.norm(grad_approx - grad_exact) / np.linalg.norm(grad_exact)

        # Get errors
        ns = 100 * (2 ** np.arange(2, 7))
        errs = np.array([get_l2_error(int(n)) for n in ns])

        # Get rates: expect h**2 convergence
        rates = np.log2(errs[:-1] / errs[1:])
        np.testing.assert_array_less(1.99 * np.ones_like(rates), rates)


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

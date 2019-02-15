#
# Test for the Finite Volume Mesh class
#
import pybamm

import numpy as np
import unittest


class TestFiniteElementMesh(unittest.TestCase):
    def test_mesh_creation(self):
        param = pybamm.ParameterValues(
            base_parameters={"Ln": 0.1, "Ls": 0.2, "Lp": 0.3}
        )
        mesh = pybamm.FiniteElementMacroMesh(param, 50)

        np.testing.assert_array_equal(mesh["whole cell"].nodes,np.linspace(0,1,50))
        np.testing.assert_array_equal(mesh["negative electrode"].nodes,np.linspace(0,1,50))

        self.assertAlmostEqual(
            np.linalg.norm(
                mesh["whole cell"].nodes
                - np.concatenate(
                    [
                        mesh["negative electrode"].nodes,
                        mesh["separator"].nodes,
                        mesh["positive electrode"].nodes,
                    ]
                )
            ),
            0,
        )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

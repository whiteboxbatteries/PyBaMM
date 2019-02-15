#
# Tests for the Reaction diffusion model
#
import pybamm
import tests

import unittest


class TestReactionDiffusionModel(unittest.TestCase):
    def test_basic_processing(self):
        model = pybamm.ReactionDiffusionModel()

        modeltest = tests.StandardModelTest(model)
        modeltest.test_all()


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

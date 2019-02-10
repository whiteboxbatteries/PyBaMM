#
# Tests for the porosity submodels
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm
import tests

import unittest


class TestPorosity(unittest.TestCase):
    def test_make_tree(self):
        j = pybamm.Scalar(1)
        pybamm.porosity.Porosity(j)

    def test_basic_processing(self):
        j = pybamm.Scalar(0.001)
        model = pybamm.porosity.Porosity(j)

        modeltest = tests.StandardModelTest(model)
        modeltest.test_all()


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

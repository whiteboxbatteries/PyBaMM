#
# Tests for the electrolyte concentrationsubmodels
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm
import tests

import unittest


class TestStefanMaxwellDiffusion(unittest.TestCase):
    def test_make_tree(self):
        G = pybamm.Scalar(1)
        pybamm.electrolyte_concentration.StefanMaxwellDiffusion(G)

    def test_basic_processing(self):
        G = pybamm.Scalar(0.001)
        model = pybamm.electrolyte_concentration.StefanMaxwellDiffusion(G)

        modeltest = tests.StandardModelTest(model)
        # Either
        # 1. run the tests individually (can pass options to individual tests)
        # modeltest.test_processing_parameters()
        # modeltest.test_processing_disc()
        # modeltest.test_solving()

        # Or
        # 2. run all the tests in one go
        modeltest.test_all()


class TestStefanMaxwellDiffusionWithPorosity(unittest.TestCase):
    def test_make_tree(self):
        j = pybamm.Scalar(1)
        pybamm.electrolyte_concentration.StefanMaxwellDiffusionWithPorosity(j)

    # def test_basic_processing(self):
    #     j = pybamm.Scalar(0.001)
    #     model = pybamm.electrolyte_concentration.StefanMaxwellDiffusionWithPorosity(j)
    #
    #     modeltest = tests.StandardModelTest(model)
    #     modeltest.test_all()


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

#
# Tests for the electrode-electrolyte interface equations
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import unittest
import numpy as np


class TestHomogeneousReaction(unittest.TestCase):
    def test_set_parameters(self):
        param = pybamm.ParameterValues(
            "input/parameters/lithium-ion/parameters/LCO.csv"
        )
        rxn = pybamm.interface.homogeneous_reaction()

        processed_rxn = param.process_symbol(rxn)

        # rxn (a concatenation of functions of scalars and parameters) should get
        # discretised to a concantenation of functions of scalars
        self.assertIsInstance(processed_rxn, pybamm.Concatenation)
        self.assertFalse(
            any(
                [
                    isinstance(x, pybamm.Parameter)
                    for x in processed_rxn.children[0].pre_order()
                ]
            )
        )
        self.assertFalse(
            any(
                [
                    isinstance(x, pybamm.Parameter)
                    for x in processed_rxn.children[2].pre_order()
                ]
            )
        )
        self.assertIsInstance(processed_rxn.children[1], pybamm.Scalar)
        self.assertEqual(processed_rxn.children[0].domain, ["negative electrode"])
        self.assertEqual(processed_rxn.children[1].domain, ["separator"])
        self.assertEqual(processed_rxn.children[2].domain, ["positive electrode"])

        # test values
        ln = param.process_symbol(pybamm.standard_parameters.ln)
        lp = param.process_symbol(pybamm.standard_parameters.lp)
        self.assertEqual(processed_rxn.children[0].evaluate() * ln.evaluate(), 1)
        self.assertEqual(processed_rxn.children[2].evaluate() * lp.evaluate(), -1)

    def test_discretisation(self):
        param = pybamm.ParameterValues(
            "input/parameters/lithium-ion/parameters/LCO.csv"
        )
        mesh = pybamm.FiniteVolumeMacroMesh(param, 2)
        disc = pybamm.FiniteVolumeDiscretisation(mesh)

        rxn = pybamm.interface.homogeneous_reaction()

        param_rxn = param.process_symbol(rxn)
        processed_rxn = disc.process_symbol(param_rxn)

        # processed_rxn should be a vector with the right shape
        self.assertIsInstance(processed_rxn, pybamm.Vector)
        self.assertEqual(processed_rxn.shape, mesh["whole cell"].nodes.shape)


class TestButlerVolmerLeadAcid(unittest.TestCase):
    def setUp(self):
        self.cn = pybamm.Variable("concentration", domain=["negative electrode"])
        self.cs = pybamm.Variable("concentration", domain=["separator"])
        self.cp = pybamm.Variable("concentration", domain=["positive electrode"])
        self.c = pybamm.Concatenation(self.cn, self.cs, self.cp)
        self.phin = pybamm.Variable(
            "potential difference", domain=["negative electrode"]
        )
        self.phis = pybamm.Variable("potential difference", domain=["separator"])
        self.phip = pybamm.Variable(
            "potential difference", domain=["positive electrode"]
        )
        self.phi = pybamm.Concatenation(self.phin, self.phis, self.phip)

    def tearDown(self):
        del self.cn
        del self.cs
        del self.cp
        del self.c
        del self.phin
        del self.phis
        del self.phip
        del self.phi

    def test_creation(self):
        # negative electrode passes, returns Multiplication
        bv = pybamm.interface.butler_volmer_lead_acid(self.cn, self.phin)
        self.assertIsInstance(bv, pybamm.Multiplication)
        self.assertEqual(bv.domain, ["negative electrode"])

        # positive electrode passes, returns Multiplication
        bv = pybamm.interface.butler_volmer_lead_acid(self.cp, self.phip)
        self.assertIsInstance(bv, pybamm.Multiplication)
        self.assertEqual(bv.domain, ["positive electrode"])

        # whole cell domain passes, retruns concatenation
        bv = pybamm.interface.butler_volmer_lead_acid(self.c, self.phi)
        self.assertIsInstance(bv, pybamm.Concatenation)
        self.assertEqual(
            bv.domain, ["negative electrode", "separator", "positive electrode"]
        )

        # None converts to whole cell
        bv = pybamm.interface.butler_volmer_lead_acid()
        self.assertIsInstance(bv, pybamm.Concatenation)
        self.assertEqual(
            bv.domain, ["negative electrode", "separator", "positive electrode"]
        )

    def test_failure(self):
        with self.assertRaises(pybamm.DomainError):
            pybamm.interface.butler_volmer_lead_acid(
                None, None, domain=["not a domain"]
            )

    def test_set_parameters(self):
        bv = pybamm.interface.butler_volmer_lead_acid(self.c, self.phi)
        param = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )
        proc_bv = param.process_symbol(bv)
        [self.assertNotIsInstance(x, pybamm.Parameter) for x in proc_bv.pre_order()]

    def test_discretisation(self):
        bv_n = pybamm.interface.butler_volmer_lead_acid(self.cn, self.phin)
        bv_p = pybamm.interface.butler_volmer_lead_acid(self.cp, self.phip)
        bv_whole = pybamm.interface.butler_volmer_lead_acid(self.c, self.phi)

        # process parameters
        param = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )
        param_bv_n = param.process_symbol(bv_n)
        param_bv_p = param.process_symbol(bv_p)
        param_bv_whole = param.process_symbol(bv_whole)

        # discretise
        mesh = pybamm.FiniteVolumeMacroMesh(param, 10)
        disc = pybamm.FiniteVolumeDiscretisation(mesh)
        y_slices = disc.get_variable_slices([self.cn, self.cp, self.phin, self.phip])
        processed_bv_n = disc.process_symbol(param_bv_n, y_slices)
        processed_bv_p = disc.process_symbol(param_bv_p, y_slices)
        processed_bv_whole = disc.process_symbol(param_bv_whole, y_slices)

        submesh = np.concatenate(
            [mesh["negative electrode"].nodes, mesh["positive electrode"].nodes]
        )
        y = np.concatenate([submesh ** 2, submesh ** 3, -submesh])

        # should evaluate to vectors with the right shape
        self.assertEqual(
            processed_bv_n.evaluate(None, y).shape,
            mesh["negative electrode"].nodes.shape,
        )
        self.assertEqual(
            processed_bv_p.evaluate(None, y).shape,
            mesh["positive electrode"].nodes.shape,
        )
        self.assertEqual(
            processed_bv_whole.evaluate(None, y).shape, mesh["whole cell"].nodes.shape
        )


class TestExchangeCurrentDensity(unittest.TestCase):
    def setUp(self):
        self.cn = pybamm.Variable("concentration", domain=["negative electrode"])
        self.cs = pybamm.Variable("concentration", domain=["separator"])
        self.cp = pybamm.Variable("concentration", domain=["positive electrode"])
        self.c = pybamm.Concatenation(self.cn, self.cs, self.cp)

    def tearDown(self):
        del self.cn
        del self.cs
        del self.cp
        del self.c

    def test_creation(self):
        # Concentration without domain passes
        c = pybamm.Variable("c")
        pybamm.interface.exchange_current_density(c, ["negative electrode"])
        pybamm.interface.exchange_current_density(c, ["positive electrode"])

        # Concentration with correct domain passes
        j0n = pybamm.interface.exchange_current_density(self.cn)
        j0p = pybamm.interface.exchange_current_density(self.cp)
        self.assertEqual(j0n.domain, ["negative electrode"])
        self.assertEqual(j0p.domain, ["positive electrode"])

        # Concentration with wrong domain fails
        with self.assertRaises(pybamm.DomainError):
            pybamm.interface.exchange_current_density(self.cp, ["negative electrode"])
        with self.assertRaises(pybamm.DomainError):
            pybamm.interface.exchange_current_density(self.cn, ["positive electrode"])

    def test_failure(self):
        with self.assertRaises(pybamm.DomainError):
            pybamm.interface.exchange_current_density(None, ["not a domain"])

    def test_set_parameters(self):
        param = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )
        j0n = pybamm.interface.exchange_current_density(self.cn)
        j0p = pybamm.interface.exchange_current_density(self.cp)
        proc_j0n = param.process_symbol(j0n)
        proc_j0p = param.process_symbol(j0p)
        [self.assertNotIsInstance(x, pybamm.Parameter) for x in proc_j0n.pre_order()]
        [self.assertNotIsInstance(x, pybamm.Parameter) for x in proc_j0p.pre_order()]

    def test_discretisation(self):
        # create exchange-current densities
        j0n = pybamm.interface.exchange_current_density(self.cn)
        j0p = pybamm.interface.exchange_current_density(self.cp)

        # process parameters
        param = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )
        param_j0n = param.process_symbol(j0n)
        param_j0p = param.process_symbol(j0p)

        # discretise
        mesh = pybamm.FiniteVolumeMacroMesh(param, 10)
        disc = pybamm.FiniteVolumeDiscretisation(mesh)
        y_slices = disc.get_variable_slices([self.cn, self.cp])
        processed_j0n = disc.process_symbol(param_j0n, y_slices)
        processed_j0p = disc.process_symbol(param_j0p, y_slices)

        submesh = np.concatenate(
            [mesh["negative electrode"].nodes, mesh["positive electrode"].nodes]
        )
        y = submesh ** 2
        # should evaluate to vectors with the right shape
        self.assertEqual(
            processed_j0n.evaluate(y=y).shape, mesh["negative electrode"].nodes.shape
        )
        self.assertEqual(
            processed_j0p.evaluate(y=y).shape, mesh["positive electrode"].nodes.shape
        )


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

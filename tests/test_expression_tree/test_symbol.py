#
# Test for the Symbol class
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import unittest
import numpy as np


class TestSymbol(unittest.TestCase):
    def test_symbol_init(self):
        sym = pybamm.Symbol("a symbol")
        self.assertEqual(sym.name, "a symbol")
        self.assertEqual(str(sym), "a symbol")

    def test_symbol_domains(self):
        a = pybamm.Symbol("a", domain=pybamm.KNOWN_DOMAINS[0])
        self.assertEqual(a.domain, [pybamm.KNOWN_DOMAINS[0]])
        a = pybamm.Symbol("a", domain=pybamm.KNOWN_DOMAINS[:2])
        self.assertEqual(a.domain, pybamm.KNOWN_DOMAINS[:2])
        with self.assertRaises(TypeError):
            a = pybamm.Symbol("a", domain=1)
        with self.assertRaises(ValueError):
            a = pybamm.Symbol("a", domain=["unknown domain"])
        with self.assertRaises(ValueError):
            a = pybamm.Symbol(
                "a", domain=[pybamm.KNOWN_DOMAINS[1], pybamm.KNOWN_DOMAINS[0]]
            )

    def test_symbol_methods(self):
        a = pybamm.Symbol("a")
        b = pybamm.Symbol("b")

        # unary
        self.assertIsInstance(-a, pybamm.Negate)
        self.assertIsInstance(abs(a), pybamm.AbsoluteValue)

        # binary - two symbols
        self.assertIsInstance(a + b, pybamm.Addition)
        self.assertIsInstance(a - b, pybamm.Subtraction)
        self.assertIsInstance(a * b, pybamm.Multiplication)
        self.assertIsInstance(a / b, pybamm.Division)
        self.assertIsInstance(a ** b, pybamm.Power)

        # binary - symbol and number
        self.assertIsInstance(a + 2, pybamm.Addition)
        self.assertIsInstance(a - 2, pybamm.Subtraction)
        self.assertIsInstance(a * 2, pybamm.Multiplication)
        self.assertIsInstance(a / 2, pybamm.Division)
        self.assertIsInstance(a ** 2, pybamm.Power)

        # binary - number and symbol
        self.assertIsInstance(3 + b, pybamm.Addition)
        self.assertEqual((3 + b).children[1].id, b.id)
        self.assertIsInstance(3 - b, pybamm.Subtraction)
        self.assertEqual((3 - b).children[1].id, b.id)
        self.assertIsInstance(3 * b, pybamm.Multiplication)
        self.assertEqual((3 * b).children[1].id, b.id)
        self.assertIsInstance(3 / b, pybamm.Division)
        self.assertEqual((3 / b).children[1].id, b.id)
        self.assertIsInstance(3 ** b, pybamm.Power)
        self.assertEqual((3 ** b).children[1].id, b.id)

        # error raising
        with self.assertRaises(NotImplementedError):
            a + "two"
        with self.assertRaises(NotImplementedError):
            a - "two"
        with self.assertRaises(NotImplementedError):
            a * "two"
        with self.assertRaises(NotImplementedError):
            a / "two"
        with self.assertRaises(NotImplementedError):
            a ** "two"
        with self.assertRaises(NotImplementedError):
            "two" + a
        with self.assertRaises(NotImplementedError):
            "two" - a
        with self.assertRaises(NotImplementedError):
            "two" * a
        with self.assertRaises(NotImplementedError):
            "two" / a
        with self.assertRaises(NotImplementedError):
            "two" ** a

    def test_multiple_symbols(self):
        a = pybamm.Symbol("a")
        b = pybamm.Symbol("b")
        c = pybamm.Symbol("c")
        exp = a * c * (a * b * c + a - c * a)
        expected_preorder = [
            "*",
            "*",
            "a",
            "c",
            "-",
            "+",
            "*",
            "*",
            "a",
            "b",
            "c",
            "a",
            "*",
            "c",
            "a",
        ]
        for node, expect in zip(exp.pre_order(), expected_preorder):
            self.assertEqual(node.name, expect)

    def test_symbol_evaluation(self):
        a = pybamm.Symbol("a")
        with self.assertRaises(NotImplementedError):
            a.evaluate()

    def test_symbol_is_constant(self):
        a = pybamm.Variable("a")
        self.assertFalse(a.is_constant())

        a = pybamm.Parameter("a")
        self.assertTrue(a.is_constant())

        a = pybamm.Scalar(1) * pybamm.Variable("a")
        self.assertFalse(a.is_constant())

        a = pybamm.Scalar(1) * pybamm.Parameter("a")
        self.assertTrue(a.is_constant())

        a = pybamm.Scalar(1) * pybamm.StateVector(slice(10))
        self.assertFalse(a.is_constant())

        a = pybamm.Scalar(1) * pybamm.Vector(np.zeros(10))
        self.assertTrue(a.is_constant())

    def test_symbol_evaluates_to_number(self):
        a = pybamm.Scalar(3)
        self.assertTrue(a.evaluates_to_number())

        a = pybamm.Parameter("a")
        self.assertFalse(a.evaluates_to_number())

        a = pybamm.Scalar(3) * pybamm.Scalar(2)
        self.assertTrue(a.evaluates_to_number())
        # highlight difference between this function and isinstance(a, Scalar)
        self.assertNotIsInstance(a, pybamm.Scalar)

        a = pybamm.Variable("a")
        self.assertFalse(a.evaluates_to_number())

        a = pybamm.Scalar(3) - 2
        self.assertTrue(a.evaluates_to_number())

        a = pybamm.Vector(np.ones(5))
        self.assertFalse(a.evaluates_to_number())

        a = pybamm.Matrix(np.ones((4, 6)))
        self.assertFalse(a.evaluates_to_number())

    def test_symbol_repr(self):
        """
        test that __repr___ returns the string
        `__class__(id, name, parent expression)`
        """
        a = pybamm.Symbol("a")
        b = pybamm.Symbol("b")
        c = pybamm.Symbol("c", domain=["test"])
        d = pybamm.Symbol("d", domain=["test"])
        hex_regex = r"\-?0x[0-9,a-f]+"
        self.assertRegex(
            a.__repr__(),
            r"Symbol\(" + hex_regex + r", a, children\=\[\], domain\=\[\]\)",
        )
        self.assertRegex(
            b.__repr__(),
            r"Symbol\(" + hex_regex + r", b, children\=\[\], domain\=\[\]\)",
        )
        self.assertRegex(
            c.__repr__(),
            r"Symbol\(" + hex_regex + r", c, children\=\[\], domain\=\['test'\]\)",
        )
        self.assertRegex(
            d.__repr__(),
            r"Symbol\(" + hex_regex + r", d, children\=\[\], domain\=\['test'\]\)",
        )
        self.assertRegex(
            (a + b).__repr__(),
            r"Addition\(" + hex_regex + r", \+, children\=\['a', 'b'\], domain=\[\]\)",
        )
        self.assertRegex(
            (c * d).__repr__(),
            r"Multiplication\("
            + hex_regex
            + r", \*, children\=\['c', 'd'\], domain=\['test'\]\)",
        )
        self.assertRegex(
            pybamm.grad(a).__repr__(),
            r"Gradient\(" + hex_regex + ", grad, children\=\['a'\], domain=\[\]\)",
        )
        self.assertRegex(
            pybamm.grad(c).__repr__(),
            r"Gradient\("
            + hex_regex
            + ", grad, children\=\['c'\], domain=\['test'\]\)",
        )

    def test_symbol_visualise(self):
        G = pybamm.Symbol("G")
        model = pybamm.electrolyte_concentration.StefanMaxwellDiffusion(G)
        c_e = list(model.rhs.keys())[0]
        rhs = model.rhs[c_e]
        rhs.visualise("StefanMaxwell_test", test=True)

    def test_has_spatial_derivatives(self):
        var = pybamm.Variable("var")
        grad_eqn = pybamm.grad(var)
        div_eqn = pybamm.div(var)
        grad_div_eqn = pybamm.div(grad_eqn)
        algebraic_eqn = 2 * var + 3
        self.assertTrue(grad_eqn.has_spatial_derivatives())
        self.assertTrue(grad_eqn.has_gradient())
        self.assertFalse(grad_eqn.has_divergence())
        self.assertTrue(div_eqn.has_spatial_derivatives())
        self.assertFalse(div_eqn.has_gradient())
        self.assertTrue(div_eqn.has_divergence())
        self.assertTrue(grad_div_eqn.has_spatial_derivatives())
        self.assertTrue(grad_div_eqn.has_gradient())
        self.assertTrue(grad_div_eqn.has_divergence())
        self.assertFalse(algebraic_eqn.has_spatial_derivatives())
        self.assertFalse(algebraic_eqn.has_gradient())
        self.assertFalse(algebraic_eqn.has_divergence())


if __name__ == "__main__":
    print("Add -v for more debug output")
    import sys

    if "-v" in sys.argv:
        debug = True
    unittest.main()

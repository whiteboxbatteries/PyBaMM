#
# Equation classes for the electrolyte porosity
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm


class Porosity(pybamm.BaseModel):
    """A class that generates the expression tree for the electrolyte porosity model

    Parameters
    ----------
    j : :class:`pybamm.Symbol`
        An expression tree that represents the concentration flux at the
        electrode-electrolyte interface

    *Extends:* :class:`BaseModel`
    """

    def __init__(self, j):
        super().__init__()

        # Variables
        eps_n = pybamm.Variable("porosity_n", domain=["negative electrode"])
        eps_s = pybamm.Variable("porosity_s", domain=["separator"])
        eps_p = pybamm.Variable("porosity_p", domain=["positive electrode"])
        eps = pybamm.Concatenation(eps_n, eps_s, eps_p)

        # Parameters
        beta_surf = pybamm.standard_parameters_lead_acid.beta_surf
        # Initial conditions
        eps_init = pybamm.standard_parameters_lead_acid.eps_init

        # Model
        self.rhs = {eps: -beta_surf * j}
        self.initial_conditions = {eps: eps_init}
        self.boundary_conditions = {}
        self.variables = {"porosity": eps}

        # Overwrite default parameter values
        self.default_parameter_values = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )

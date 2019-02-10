#
# Equation classes for the electrolyte concentration
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm


class StefanMaxwellDiffusion(pybamm.BaseModel):
    """A class that generates the expression tree for Stefan-Maxwell Diffusion in the
    electrolyte.

    Parameters
    ----------
    G : :class:`pybamm.Symbol`
        An expression tree that represents the concentration flux at the
        electrode-electrolyte interface

    Attributes
    ----------

    rhs: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the rhs
    initial_conditions: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the initial conditions
    boundary_conditions: dict
        A dictionary that maps expressions (variables) to expressions that represent
        the boundary conditions
    variables: dict
        A dictionary that maps strings to expressions that represent
        the useful variables

    *Extends:* :class:`BaseModel`
    """

    def __init__(self, G):
        super().__init__()

        epsilon = pybamm.standard_parameters.epsilon_s  # make issue for spatially
        # dependent parameters
        b = pybamm.standard_parameters.b
        delta = pybamm.standard_parameters.delta
        nu = pybamm.standard_parameters.nu
        t_plus = pybamm.standard_parameters.t_plus
        ce0 = pybamm.standard_parameters.ce0

        whole_cell = ["negative electrode", "separator", "positive electrode"]

        c_e = pybamm.Variable("c_e", domain=whole_cell)

        N_e = -(epsilon ** b) * pybamm.grad(c_e)

        self.rhs = {c_e: -pybamm.div(N_e) / delta / epsilon + nu * (1 - t_plus) * G}
        self.initial_conditions = {c_e: ce0}
        self.boundary_conditions = {
            N_e: {"left": pybamm.Scalar(0), "right": pybamm.Scalar(0)}
        }
        self.variables = {"c_e": c_e, "N_e": N_e}


class StefanMaxwellDiffusionWithPorosity(pybamm.BaseModel):
    """A class that generates the expression tree for Stefan-Maxwell Diffusion in the
    electrolyte.

    Parameters
    ----------
    j : :class:`pybamm.Symbol`
        An expression tree that represents the concentration flux at the
        electrode-electrolyte interface

    *Extends:* :class:`BaseModel`
    """

    def __init__(self, j):
        super().__init__()

        # Domains
        whole_cell = ["negative electrode", "separator", "positive electrode"]

        # Variables
        c = pybamm.Variable("concentration", domain=whole_cell)
        eps_n = pybamm.Variable("porosity_n", domain=["negative electrode"])
        eps_s = pybamm.Variable("porosity_s", domain=["separator"])
        eps_p = pybamm.Variable("porosity_p", domain=["positive electrode"])
        eps = pybamm.Concatenation(eps_n, eps_s, eps_p)

        # Parameters
        D = pybamm.standard_parameters_lead_acid.D
        Cd = pybamm.standard_parameters_lead_acid.Cd
        s = pybamm.standard_parameters_lead_acid.s
        beta_surf = pybamm.standard_parameters_lead_acid.beta_surf
        # Initial conditions
        c_init = pybamm.standard_parameters_lead_acid.c_init

        # Model
        # flux
        N = -(D(c) * eps ** 1.5) * pybamm.grad(c)
        # cation conservation equation
        dcdt = 1 / eps * (1 / Cd * pybamm.div(N) + (s + beta_surf) * j)

        self.rhs = {c: dcdt}
        self.initial_conditions = {c: c_init}
        self.boundary_conditions = {N: {"left": 0, "right": 0}}
        self.variables = {"concentration": c, "flux": N, "porosity": eps}

        # Overwrite default parameter values
        self.default_parameter_values = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )

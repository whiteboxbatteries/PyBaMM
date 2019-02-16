#
# Simulation class for a battery model
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import argparse
import numpy as np


class Simulation(object):
    """
    The simulation class for a battery model.

    Parameters
    ---------
    model : pybamm.models.(modelname).(ModelName)() instance
       The model to be used for the simulation. (modelname) and (ModelName)
       refer to a module and class to be chosen.
    parameter_values : :class:`pybamm.ParameterValues.Parameters` instance
       The parameters to be used for the simulation.
    discretisation : :class:`pybamm.discretisation.Mesh` instance
       The discretisation to be used for the simulation.
    solver : :class:`pybamm.solver.Solver` instance
       The algorithm for solving the model.
    name : string, optional
       The simulation name.

    Examples
    --------
    >>> import pybamm
    >>> model = pybamm.lead_acid.LOQS()
    >>> sim = pybamm.Simulation(model)
    >>> sim.run()

    """

    def __init__(
        self,
        model,
        parameter_values=None,
        geometry=None,
        submesh_types=None,
        submesh_pts=None,
        spatial_methods=None,
        solver=None,
        name="unnamed",
    ):
        # Read defaults from model
        if parameter_values is None:
            parameter_values = model.default_parameter_values
        if geometry is None:
            geometry = model.default_geometry
        if submesh_types is None:
            submesh_types = model.default_submesh_types
        if submesh_pts is None:
            submesh_pts = model.default_submesh_pts
        if solver is None:
            solver = model.default_solver
        if solver is None:
            solver = model.default_solver

        # Assign attributes
        self.model = model
        self.parameter_values = parameter_values
        parameter_values.process_geometry(geometry)
        self.geometry = geometry
        self.mesh = pybamm.Mesh(geometry, submesh_types, submesh_pts)
        self.discretisation = pybamm.Discretisation(self.mesh, spatial_methods)
        self.solver = solver
        self.name = name

        # Default evaluation time
        self.default_t_eval = np.linspace(0, 1)

    def __str__(self):
        return self.name

    def run(self, t_eval=None):
        self.parameter_values.process_model(self.model)
        self.discretisation.process_model(self.model)
        if t_eval is None:
            t_eval = self.default_t_eval
        self.solver.solve(self.model, t_eval)

    def load(self):
        raise NotImplementedError

    def save(self):
        raise NotImplementedError

    def plot(self):
        raise NotImplementedError

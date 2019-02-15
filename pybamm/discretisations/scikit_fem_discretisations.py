#
# Finite Volume discretisation class
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import skfem
import numpy as np


class ScikitFiniteElementDiscretisation(pybamm.BaseDiscretisation):
    """Discretisation using Finite Element method using scikit.fem.

    Parameters
    ----------
    mesh : :class:`pybamm.BaseMesh` (or subclass)
        The underlying mesh for discretisation

    **Extends:** :class:`pybamm.BaseDiscretisation`
    """

    def __init__(self, mesh_type, submesh_pts, submesh_types):
        super().__init__(mesh_type, submesh_pts, submesh_types)

    def compute_diffusivity(self, symbol):
        """
        This is just for finite volume, return symbol unchanged

        Parameters
        ----------
        symbol : :class:`pybamm.Symbol`

        Returns
        -------
        symbol : :class:`pybamm.Symbol`
        """

        return symbol

    def gradient(self, symbol, y_slices, boundary_conditions):
        """Matrix-vector multiplication to implement the gradient operator.
        See :meth:`pybamm.BaseDiscretisation.gradient`
        """
        # Check that boundary condition keys are hashes (ids)
        for key in boundary_conditions.keys():
            assert isinstance(key, int), TypeError(
                "boundary condition keys should be hashes, not {}".format(type(key))
            )
        # Discretise symbol
        discretised_symbol = self.process_symbol(symbol, y_slices, boundary_conditions)
        domain = symbol.domain
        # Add Dirichlet boundary conditions, if defined
        if symbol.id in boundary_conditions:
            lbc = boundary_conditions[symbol.id]["left"]
            rbc = boundary_conditions[symbol.id]["right"]
            discretised_symbol = self.add_ghost_nodes(discretised_symbol, lbc, rbc)
            domain = (
                [domain[0] + "_left ghost cell"]
                + domain
                + [domain[-1] + "_right ghost cell"]
            )
        gradient_matrix = self.gradient_matrix(domain)
        return gradient_matrix * discretised_symbol

    def gradient_matrix(self, domain):

        @skfem.bilinear_form
        def fem_gradient(u, du, v, dv, w):
            return du[0]*v[0]


        # Create appropriate submesh by combining submeshes in domain
        submesh = self.mesh.combine_submeshes(*domain)

        e = skfem.ElementTriP1()
        basis = skfem.InteriorBasis(submesh.fem_mesh, e)
        A = skfem.asm(laplace, basis)

        return pybamm.Matrix(A)


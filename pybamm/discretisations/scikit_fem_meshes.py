#
# Mesh class for space and time discretisation
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm

import numpy as np

import skfem



class ScikitFemMesh1D:
    """
    A 1D Finite Element mesh based on scikit.fem

    Parameters
    ----------

    lims :class: `Geometry`
        contains the geometry of the problem
    npts: int
        number of points in the mesh

    """
    def __init__(self, lims, npts):
        # currently accept lims and npts as dicts. This may get changed at a future
        # date depending on the form of mesh we desire for 2D/3D

        spatial_lims = list(lims.values())[0]
        npts = list(npts.values())[0]
        edges = np.linspace(spatial_lims["min"], spatial_lims["max"], npts + 1)
        self._mesh = skfem.MeshLine(edges)

    @property
    def fem_mesh(self):
        return self._mesh



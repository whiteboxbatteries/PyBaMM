#
# Native PyBaMM Meshes
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals

import numpy as np


class SubMesh1D:
    """ Submesh class.
        Contains the position of the nodes and the number of mesh points.

        Parameters
        ----------
        domain : dict
            A dictionary that contains the limits of the spatial variables
        npts : dict
            A dictionary that contains the number of points to be used on each
            spatial variable
        """

    def __init__(self, edges):
        self.edges = edges
        self.nodes = (self.edges[1:] + self.edges[:-1]) / 2
        self.d_edges = np.diff(self.edges)
        self.d_nodes = np.diff(self.nodes)
        self.npts = self.nodes.size


class Uniform1DSubMesh(SubMesh1D):
    """
    A class to generate a uniform submesh on a 1D domain
    """

    def __init__(self, lims, npts):

        # currently accept lims and npts as dicts. This may get changed at a future
        # date depending on the form of mesh we desire for 2D/3D

        spatial_lims = list(lims.values())[0]
        npts = list(npts.values())[0]

        edges = np.linspace(spatial_lims["min"], spatial_lims["max"], npts + 1)

        super().__init__(edges)

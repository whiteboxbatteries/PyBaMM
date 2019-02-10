#
# Lead-acid full 1D DAE model
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm


class Full1DDae(pybamm.BaseModel):
    """Full 1D DAE model for lead-acid, from [1].

    .. math::
        \\frac{\\partial \\tilde{c}}{\\partial t}
        = \\frac{1}{\\varepsilon^{(0)}}\\left(
            \\frac{D^{\\text{eff}, (0)}}{\\mathcal{C}_\\text{d}}
            \\frac{\\partial^2 \\tilde{c}}{\\partial x^2}
            + \\left(
                s + \\beta^{\\text{surf}}c^{(0)}
            \\right)j^{(0)}
        \\right)


    [1] Sulzer, V., Chapman, S. J., Please, C. P., Howey, D. A., & Monroe, C. W.
    (2019). Faster Lead-Acid Battery Simulations from Porous-Electrode Theory: I.
    Physical Model. arXiv preprint arXiv:1902.01771. `full text
    <https://arxiv.org/pdf/1902.01771.pdf>`_.

    **Notation for variables and parameters:**

    * f_xy means :math:`f^{(x)}_\\text{y}` (x is the power for the asymptotic \
    expansion and y is the domain). For example c_1n means :math:`c^{(1)}_\\text{n}`, \
    the first-order concentration in the negative electrode
    * fbar_n means :math:`\\bar{f}_n`, the average value of f in that domain, e.g.

    .. math::
        \\text{cbar_n}
        = \\bar{c}_\\text{n}
        = \\frac{1}{\\ell_\\text{n}}
        \\int_0^{\\ell_\\text{n}} \\! c_\\text{n} \\, \\mathrm{d}x

    **Extends:** :class:`pybamm.BaseModel`

    """

    def __init__(self):
        super().__init__()

        # Interfacial current density: Butler-Volmer
        j = pybamm.interface.butler_volmer_lead_acid()

        # Submodels for the electrolyte
        conc = pybamm.electrolyte_concentration.StefanMaxwellDiffusionWithPorosity(j)
        porosity = pybamm.electrolyte_concentration.Porosity(j)

        # Create from submodels
        self.update(conc, porosity)

        # Overwrite default parameter values
        self.default_parameter_values = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )

#
# Lead-acid full 1D DAE model
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm


class Full1DDae(pybamm.BaseModel):
    """Full 1D DAE model for lead-acid.

    .. math::
        \\frac{\\partial \\tilde{c}}{\\partial t}
        = \\frac{1}{\\varepsilon^{(0)}}\\left(
            \\frac{D^{\\text{eff}, (0)}}{\\mathcal{C}_\\text{d}}
            \\frac{\\partial^2 \\tilde{c}}{\\partial x^2}
            + \\left(
                s + \\beta^{\\text{surf}}c^{(0)}
            \\right)j^{(0)}
        \\right)


    [1] Paper reference

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

        # Variables
        c = pybamm.Variable("c", domain=[])
        epsn = pybamm.Variable("epsn", domain=[])
        epss = pybamm.Variable("epss", domain=[])
        epsp = pybamm.Variable("epsp", domain=[])

        # Parameters
        ln = pybamm.standard_parameters.ln
        ls = pybamm.standard_parameters.ls
        lp = pybamm.standard_parameters.lp
        sn = pybamm.standard_parameters_lead_acid.sn
        sp = pybamm.standard_parameters_lead_acid.sp
        beta_surf_n = pybamm.standard_parameters_lead_acid.beta_surf_n
        beta_surf_p = pybamm.standard_parameters_lead_acid.beta_surf_p
        iota_ref_n = pybamm.standard_parameters_lead_acid.iota_ref_n
        iota_ref_p = pybamm.standard_parameters_lead_acid.iota_ref_p
        U_Pb = pybamm.standard_parameters_lead_acid.U_Pb_ref
        U_PbO2 = pybamm.standard_parameters_lead_acid.U_PbO2_ref
        # Initial conditions
        c_init = pybamm.standard_parameters_lead_acid.c_init
        epsn_init = pybamm.standard_parameters_lead_acid.epsn_init
        epss_init = pybamm.standard_parameters_lead_acid.epss_init
        epsp_init = pybamm.standard_parameters_lead_acid.epsp_init

        # ODEs
        jn = pybamm.standard_parameters_lead_acid.icell / ln
        jp = -pybamm.standard_parameters_lead_acid.icell / lp
        depsndt = -beta_surf_n * jn
        depspdt = -beta_surf_p * jp
        dcdt = (
            1
            / (ln * epsn + ls * epss + lp * epsp)
            * ((sn - sp) - c * (ln * depsndt + lp * depspdt))
        )
        self.rhs = {c: dcdt, epsn: depsndt, epss: pybamm.Scalar(0), epsp: depspdt}
        # Initial conditions
        self.initial_conditions = {
            c: c_init,
            epsn: epsn_init,
            epss: epss_init,
            epsp: epsp_init,
        }
        # ODE model -> no boundary conditions
        self.boundary_conditions = {}

        # Variables
        Phi = -U_Pb - jn / (2 * iota_ref_n * c)
        V = Phi + U_PbO2 - jp / (2 * iota_ref_p * c)
        # Phisn = pybamm.Scalar(0)
        # Phisp = V
        # Concatenate variables
        # eps = pybamm.Concatenation(epsn, epss, epsp)
        # Phis = pybamm.Concatenation(Phisn, pybamm.Scalar(0), Phisp)
        # self.variables = {"c": c, "eps": eps, "Phi": Phi, "Phis": Phis, "V": V}
        self.variables = {
            "c": pybamm.Broadcast(c, ["whole cell"]),
            "Phi": pybamm.Broadcast(Phi, ["whole cell"]),
            "V": V,
            "int(epsilon_times_c)dx": (ln * epsn + ls * epss + lp * epsp) * c,
        }

        # Overwrite default parameter values
        self.default_parameter_values = pybamm.ParameterValues(
            "input/parameters/lead-acid/default.csv", {"current scale": 1}
        )

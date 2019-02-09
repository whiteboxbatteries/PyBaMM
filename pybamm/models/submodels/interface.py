#
# Equations for the electrode-electrolyte interface
#
from __future__ import absolute_import, division
from __future__ import print_function, unicode_literals
import pybamm


def homogeneous_reaction():
    """
    Homogeneous reaction at the electrode-electrolyte interface
    """

    current_neg = (
        pybamm.Scalar(1, domain=["negative electrode"]) / pybamm.standard_parameters.ln
    )
    current_sep = pybamm.Scalar(0, domain=["separator"])
    current_pos = (
        -pybamm.Scalar(1, domain=["positive electrode"]) / pybamm.standard_parameters.lp
    )
    return pybamm.Concatenation(current_neg, current_sep, current_pos)


def butler_volmer_lead_acid(domain=None, system="dae"):
    """
    Butler-Volmer reactions for lead-acid chemistry

    .. math::
        j = j_0(c) * \\sinh(\\eta),

        \\text{where} \\eta = \\Phi_\\text{s} - \\Phi - U(c)

    Parameters
    ----------
    domain : iterable of strings
        The domain in which to calculate the interfacial current density
    system : string
        Whether the model being solved is in the DAE form (e.g.
        :class:`pybamm.lead_acid.Full1DDAE` or the PDE form (e.g.
        :class:`pybamm.lead_acid.Full1DPDE`). If "dae", the variables are c, Phi and
        Phis. If PDE, the variables are c and eta. Default is "dae".


    Returns
    -------
    :class:`pybamm.Symbol`
        The interfacial current density in the appropriate domain
    """
    # avoid passing list (mutable) as default argument
    if domain is None:
        domain = ["negative electrode", "separator", "positive electrode"]

    # Variables
    cn = pybamm.Variable("concentration", domain=["negative electrode"])
    cp = pybamm.Variable("concentration", domain=["positive electrode"])
    if system == "dae":
        Phin = pybamm.Variable("electrolyte potential", domain=["negative electrode"])
        Phip = pybamm.Variable("electrolyte potential", domain=["positive electrode"])
        Phisn = pybamm.Variable("solid potential", domain=["negative electrode"])
        Phisp = pybamm.Variable("solid potential", domain=["positive electrode"])
        etan = Phisn - Phin - pybamm.standard_parameters_lead_acid.U_Pb(cn)
        etap = Phisp - Phip - pybamm.standard_parameters_lead_acid.U_PbO2(cp)
    elif system == "pde":
        etan = pybamm.Variable("overpotential", domain=["negative electrode"])
        etap = pybamm.Variable("overpotential", domain=["positive electrode"])
    else:
        raise ValueError("system '{}' not recognised".format(system))
    # Exchange-current densities
    j0n = exchange_current_density(cn, ["negative electrode"])
    j0p = exchange_current_density(cp, ["positive electrode"])

    # Interfacial current density
    current_neg = j0n * etan  # Function(etan, np.sinh)
    current_pos = j0p * etap  # Function(etap, np.sinh)
    if domain == ["negative electrode"]:
        return current_neg
    elif domain == ["positive electrode"]:
        return current_pos
    elif domain == ["negative electrode", "separator", "positive electrode"]:
        current_sep = pybamm.Scalar(0, domain=["separator"])
        return pybamm.Concatenation(current_neg, current_sep, current_pos)
    else:
        raise pybamm.DomainError("domain '{}' not recognised".format(domain))


def exchange_current_density(c, domain):
    """The exchange current-density as a function of concentration

    Parameters
    ----------
    c : :class:`pybamm.Variable`
        A Variable representing the concentration
    domain : string
        Which domain to calculate the exchange current density in ("negative electrode"
        or "positive electrode")

    Returns
    -------
    :class:`pybamm.Symbol`
        The exchange-current density
    """
    if domain == ["negative electrode"]:
        # concentration domain should be empty or the same as domain
        if c.domain not in [["negative electrode"], []]:
            raise pybamm.DomainError("""concentration and domain do not match""")
        iota_ref_n = pybamm.standard_parameters_lead_acid.iota_ref_n
        return iota_ref_n * c
    elif domain == ["positive electrode"]:
        # concentration domain should be empty or the same as domain
        if c.domain not in [["positive electrode"], []]:
            raise pybamm.DomainError("""concentration and domain do not match""")
        iota_ref_p = pybamm.standard_parameters_lead_acid.iota_ref_p
        Ve = pybamm.standard_parameters_lead_acid.Ve
        Vw = pybamm.standard_parameters_lead_acid.Vw
        cw = (1 - c * Ve) / Vw
        return iota_ref_p * c ** 2 * cw
    else:
        raise pybamm.DomainError("domain '{}' not recognised".format(domain))

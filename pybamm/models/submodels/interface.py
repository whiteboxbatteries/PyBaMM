#
# Equations for the electrode-electrolyte interface
#
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

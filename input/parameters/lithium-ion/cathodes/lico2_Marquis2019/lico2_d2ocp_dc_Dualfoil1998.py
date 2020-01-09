from pybamm import cosh, sech, exp, tanh


def lico2_d2ocp_dc_Dualfoil1998(sto, cpmax):
    """
    Lithium Cobalt Oxide (LiCO2) Open Circuit Potential (OCP) as a a function of the
    stochiometry. The fit is taken from Dualfoil [1]. Dualfoil states that the data
    was measured by Oscar Garcia 2001 using Quallion electrodes for 0.5 < sto < 0.99
    and by Marc Doyle for sto<0.4 (for unstated electrodes). We could not find any
    other records of the Garcia measurements. Doyles fits can be found in his
    thesis [2] but we could not find any other record of his measurments.

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html
    .. [2] CM Doyle. Design and simulation of lithium rechargeable batteries,
           1995.

    Parameters
    ----------
    sto: double
       Stochiometry of material (li-fraction)

    """

    cp = sto * cpmax

    d2udc = (
        -(
            (
                511.847
                * sech(30.834 - (57.8584 * cp) / cpmax) ** 2
                * tanh(30.834 - (57.8584 * cp) / cpmax)
            )
            / cpmax ** 2
        )
        - (
            12313.5
            * sech(52.294 - (53.4122 * cp) / cpmax) ** 2
            * tanh(52.294 - (53.4122 * cp) / cpmax)
        )
        / cpmax ** 2
        + (
            125.987
            * sech(11.0923 - (21.0853 * cp) / cpmax) ** 2
            * tanh(11.0923 - (21.0853 * cp) / cpmax)
        )
        / cpmax ** 2
        - (
            13.938
            * sech(1.4684 - (5.82911 * cp) / cpmax) ** 2
            * tanh(1.4684 - (5.82911 * cp) / cpmax)
        )
        / cpmax ** 2
        - (
            32.9655
            * sech(7.59878*(0.56478 - (1.062 * cp) / cpmax)) ** 2
            * tanh(7.59878*(0.56478 - (1.062 * cp) / cpmax))
        )
        / cpmax ** 2
        + (
            1357.8
            * sech(166.667*(-0.525 + (1.062 * cp) / cpmax)) ** 2
            * tanh(166.667*(-0.525 + (1.062 * cp) / cpmax))
        )
        / cpmax ** 2
    )

    return d2udc

from pybamm import cosh, sech, exp


def lico2_docp_dc_Dualfoil1998(sto, cpmax):
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

    # stretch = 1.062
    # sto = stretch * sto

    # ct3 = 3.6e3 * 247 * 5010 / 96487

    cp = sto * cpmax

    # g1p = (
    #     0.07645 * (-54.4806 / ct3) * ((1.0 / cosh(30.834 - 54.4806 * sto)) ** 2)
    #     + 2.1581 * (-50.294 / ct3) * ((cosh(52.294 - 50.294 * sto)) ** (-2))
    #     + 0.14169 * (19.854 / ct3) * ((cosh(11.0923 - 19.8543 * sto)) ** (-2))
    #     - 0.2051 * (5.4888 / ct3) * ((cosh(1.4684 - 5.4888 * sto)) ** (-2))
    #     - 0.2531 / 0.1316 / ct3 * ((cosh((-sto + 0.56478) / 0.1316)) ** (-2))
    #     - 0.02167 / 0.006 / ct3 * ((cosh((sto - 0.525) / 0.006)) ** (-2))
    # )
    dudc = (
        -((4.42327 * sech(30.834 - (57.8584 * cp) / cpmax) ** 2) / cpmax)
        - (115.269 * sech(52.294 - (53.4122 * cp) / cpmax) ** 2) / cpmax
        + (2.98757 * sech(11.0923 - (21.0853 * cp) / cpmax) ** 2) / cpmax
        - (1.19555 * sech(1.4684 - (5.82911 * cp) / cpmax) ** 2) / cpmax
        - (2.04249 * sech(7.59878 * (0.56478 - (1.062 * cp) / cpmax)) ** 2) / cpmax
        - (3.83559 * sech(166.667 * (-0.525 + (1.062 * cp) / cpmax)) ** 2) / cpmax
    )

    return dudc

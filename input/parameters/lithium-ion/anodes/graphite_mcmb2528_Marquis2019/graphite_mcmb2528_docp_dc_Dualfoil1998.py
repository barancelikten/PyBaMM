from pybamm import cosh, exp


def graphite_mcmb2528_docp_dc_Dualfoil1998(sto):
    """
       Graphite MCMB 2528 Open Circuit Potential (OCP) as a function of the
       stochiometry. The fit is taken from Dualfoil [1]. Dualfoil states that the data
       was measured by Chris Bogatu at Telcordia and PolyStor materials, 2000. However,
       we could not find any other records of this measurment.

       References
       ----------
       .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html
       """

    ct1 = 3.6e3 * 372 * 1800 / 96487

    g1n = (
        -1.5 * (120.0 / ct1) * exp(-120 * sto)
        + (0.0351 / (0.083 * ct1)) * ((cosh((sto - 0.286) / 0.083)) ** (-2))
        - (0.0045 / (ct1 * 0.119)) * ((cosh((sto - 0.849) / 0.119)) ** (-2))
        - (0.035 / (ct1 * 0.05)) * ((cosh((sto - 0.9233) / 0.05)) ** (-2))
        - (0.0147 / (ct1 * 0.034)) * ((cosh((sto - 0.5) / 0.034)) ** (-2))
        - (0.102 / (ct1 * 0.142)) * ((cosh((sto - 0.194) / 0.142)) ** (-2))
        - (0.022 / (ct1 * 0.0164)) * ((cosh((sto - 0.9) / 0.0164)) ** (-2))
        - (0.011 / (ct1 * 0.0226)) * ((cosh((sto - 0.124) / 0.0226)) ** (-2))
        + (0.0155 / (ct1 * 0.029)) * ((cosh((sto - 0.105) / 0.029)) ** (-2))
    )

    return g1n

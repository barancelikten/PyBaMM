from pybamm import exp, sech


def graphite_mcmb2528_docp_dc_Dualfoil1998(sto, cnmax):
    """
       Graphite MCMB 2528 Open Circuit Potential (OCP) as a function of the
       stochiometry. The fit is taken from Dualfoil [1]. Dualfoil states that the data
       was measured by Chris Bogatu at Telcordia and PolyStor materials, 2000. However,
       we could not find any other records of this measurment.

       References
       ----------
       .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html
       """

    cn = sto * cnmax

    dudc = (
        -((180 * exp(-((120 * cn) / cnmax))) / cnmax)
        - (0.7 * sech(20.0 * (-0.9233 + cn / cnmax)) ** 2) / cnmax
        - (1.34146 * sech(60.9756 * (-0.9 + cn / cnmax)) ** 2) / cnmax
        - (0.0378151 * sech(8.40336 * (-0.849 + cn / cnmax)) ** 2) / cnmax
        - (0.432353 * sech(29.4118 * (-0.5 + cn / cnmax)) ** 2) / cnmax
        + (0.422892 * sech(12.0482 * (-0.286 + cn / cnmax)) ** 2) / cnmax
        - (0.71831 * sech(7.04225 * (-0.194 + cn / cnmax)) ** 2) / cnmax
        - (0.486726 * sech(44.2478 * (-0.124 + cn / cnmax)) ** 2) / cnmax
        + (0.534483 * sech(34.4828 * (-0.105 + cn / cnmax)) ** 2) / cnmax
    )

    return dudc

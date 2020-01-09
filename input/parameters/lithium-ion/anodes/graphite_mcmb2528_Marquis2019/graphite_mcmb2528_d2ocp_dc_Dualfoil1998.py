from pybamm import cosh, exp, sech, tanh


def graphite_mcmb2528_d2ocp_dc_Dualfoil1998(sto, cnmax):
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

    d2udc = (
        (21600 * exp(-((120 * cn) / cnmax))) / cnmax ** 2
        + (
            28
            * sech(20 * (-0.9233 + cn / cnmax)) ** 2
            * tanh(20.0 * (-0.9233 + cn / cnmax))
        )
        / cnmax ** 2
        + (
            163.593
            * sech(60.9756 * (-0.9 + cn / cnmax)) ** 2
            * tanh(60.9756 * (-0.9 + cn / cnmax))
        )
        / cnmax ** 2
        + (
            0.635548
            * sech(8.40336 * (-0.849 + cn / cnmax)) ** 2
            * tanh(8.40336 * (-0.849 + cn / cnmax))
        )
        / cnmax ** 2
        + (
            25.4325
            * sech(29.4118 * (-0.5 + cn / cnmax)) ** 2
            * tanh(29.4118 * (-0.5 + cn / cnmax))
        )
        / cnmax ** 2
        - (
            10.1902
            * sech(12.0482 * (-0.286 + cn / cnmax)) ** 2
            * tanh(12.0482 * (-0.286 + cn / cnmax))
        )
        / cnmax ** 2
        + (
            10.117
            * sech(7.04225 * (-0.194 + cn / cnmax)) ** 2
            * tanh(7.04225 * (-0.194 + cn / cnmax))
        )
        / cnmax ** 2
        + (
            43.0731
            * sech(44.2478 * (-0.124 + cn / cnmax)) ** 2
            * tanh(44.2478 * (-0.124 + cn / cnmax))
        )
        / cnmax ** 2
        - (
            36.8609
            * sech(34.4828 * (-0.105 + cn / cnmax)) ** 2
            * tanh(34.4828 * (-0.105 + cn / cnmax))
        )
        / cnmax ** 2
    )

    return d2udc

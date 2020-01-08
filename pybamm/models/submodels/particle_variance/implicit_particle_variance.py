#
# Class for calculating the variance the particle concentrations as
# already defined in the model
#

import pybamm
from .base_particle_variance import BaseModel


class Implicit(BaseModel):
    """Class for calculating the variance in the particle
    concentrations using the particle concentrations already
    provided by the "Particle" submodel.

    Parameters
    ----------
    param : parameter class
        The parameters to use for this submodel
    domain : str
        The domain of the model either 'Negative' or 'Positive'


    **Extends:** :class:`pybamm.BaseSubModel`
    """

    def __init__(self, param, domain):
        super().__init__(param, domain)

    def get_coupled_variables(self, variables):

        c = variables[self.domain + " particle concentration"]
        c_av = variables[
            "X-averaged " + self.domain.lower() + " particle concentration"
        ]

        c_surf = pybamm.surf(c)
        c_av_surf = pybamm.surf(c_av)

        self.get_standard_variance_variables(c_surf - c_av_surf)


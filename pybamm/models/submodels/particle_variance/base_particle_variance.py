#
# Base class for calculating the variance the particle concentrations
#

import pybamm


class BaseModel(pybamm.BaseSubModel):
    """Base class for calculating the variance in the particle
    concentrations.

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

    def get_standard_variance_variables(self, c):

        x_n = pybamm.standard_spatial_vars.x_n
        x_p = pybamm.standard_spatial_vars.x_p

        if self.domain == "Negative":
            c_var = pybamm.Integral(c ** 2, x_n) / self.param.l_n
        elif self.domain == "Positive":
            c_var = pybamm.Integral(c ** 2, x_p) / self.param.l_p

        variables = {
            self.domain + " particle concentration variance": c_var,
            self.domain + " particle surface concentration difference": c,
        }

        return variables

    def get_standard_basis_variables(self, xi_1, xi_2):

        variables = {
            "First " + self.domain.lower() + " particle basis function": xi_1,
            "Second " + self.domain.lower() + " particle basis function": xi_2,
            "Surface value of first "
            + self.domain.lower()
            + " particle basis function": pybamm.surf(xi_1),
            "Surface value of second "
            + self.domain.lower()
            + " particle basis function": pybamm.surf(xi_2),
        }

        return variables


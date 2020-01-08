#
# Class for calculating the variance the particle concentrations using
# the ad-hoc reduced fickian model
#

import pybamm
from .base_particle_variance import BaseModel


class AdHocFickian(BaseModel):
    """Class for calculating the variance in the particle
    concentrations using the ad-hoc reduced fickian model.

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

    def get_fundamental_variables(self):
        xi_1 = pybamm.Variable(
            "First " + self.domain.lower() + " particle basis function",
            domain=self.domain.lower() + " particle",
        )
        xi_2 = pybamm.Variable(
            "Second " + self.domain.lower() + " particle basis function",
            domain=self.domain.lower() + " particle",
        )

        self.get_standard_basis_variables(xi_1, xi_2)

    def get_coupled_variables(self, variables):

        xi_1_surf = variables[
            "Surface value of first " + self.domain.lower() + " particle basis function"
        ]
        xi_2_surf = variables[
            "Surface value of second "
            + self.domain.lower()
            + " particle basis function"
        ]

        # parameters and variables for c_e_k_hat, phi_e_hat, and phi_s_hat
        t_plus = self.param.t_plus
        gamma_e = self.param.gamma_e
        l_p = self.param.l_p
        l_s = self.param.l_s
        l_n = self.param.l_n

        tor_n_av = variables["Leading-order x-averaged negative electrolyte tortuosity"]
        tor_s_av = variables["Leading-order x-averaged separator tortuosity"]
        tor_p_av = variables["Leading-order x-averaged positive electrolyte tortuosity"]

        x_n = pybamm.standard_spatial_vars.x_n
        x_p = pybamm.standard_spatial_vars.x_p

        c_e_av = variables["X-averaged electrolyte concentration"]
        T_av = variables["X-averaged cell temperature"]

        D_e = self.param.D_e(c_e_av, T_av)

        kappa_n_av = self.param.kappa_e(c_e_av, T_av)
        kappa_p_av = self.param.kappa_e(c_e_av, T_av)

        tor_n_0 = variables["Leading-order x-averaged negative electrode tortuosity"]
        tor_p_0 = variables["Leading-order x-averaged positive electrode tortuosity"]

        c_e_n_hat = (
            (1 - t_plus)
            / (gamma_e * 6 * D_e)
            * (
                2 * ((l_p ** 2 / tor_p_av) - (l_n ** 2 / tor_n_av))
                + 3 * l_s / tor_s_av * (1 - l_p - l_n)
                + 3 / (tor_n_0 * l_n) * (l_n ** 2 - x_n ** 2)
            )
        )
        c_e_p_hat = (
            (1 - t_plus)
            / (gamma_e * 6 * D_e)
            * (
                2 * ((l_p ** 2 / tor_p_av) - (l_n ** 2 / tor_n_av))
                + 3 * l_s / tor_s_av * (l_p - l_n - 1)
                + 3 / (tor_p_0 * l_p) * ((x_p - 1) ** 2 - l_p ** 2)
            )
        )

        phi_e_n_hat = (
            1
            / gamma_e
            / kappa_n_av
            * ((x_n ** 2 - l_n ** 2) / (2 * tor_n_av * l_n) + l_n / tor_s_av)
        )
        phi_e_p_hat = (
            1
            / gamma_e
            / kappa_p_av
            * (
                (x_p * (2 - x_p) + l_p ** 2 - 1) / (2 * tor_p_av * l_p)
                + (1 - l_p) / tor_s_av
            )
        )

        phi_s_n_hat = (
            x_n * (x_n - 2 * l_n) / (2 * self.param.sigma_n_prime * tor_n_0 * l_n)
        )
        phi_s_p_hat = x_p + (x_p - 1) ** 2 / (
            2 * self.param.sigma_p_prime * tor_p_0 * l_p
        )

        if self.domain == "Negative":
            beta_1 = c_e_n_hat
            beta_2 = phi_s_n_hat - phi_e_n_hat
            beta_1_av = pybamm.Integral(beta_1, x_n) / l_n
            beta_2_av = pybamm.Integral(beta_2, x_n) / l_n
        elif self.domain == "Positive":
            beta_1 = c_e_p_hat
            beta_2 = phi_s_p_hat - phi_e_p_hat
            beta_1_av = pybamm.Integral(beta_1, x_p) / l_p
            beta_2_av = pybamm.Integral(beta_2, x_p) / l_p

        c_surf = (beta_1 - beta_1_av) * xi_1_surf + (beta_2 - beta_2_av) * xi_2_surf

        self.get_standard_variance_variables(c_surf)

    def set_rhs(self, variables):

        c_s_av = variables[
            "X-averaged " + self.domain.lower() + " particle concentration"
        ]
        xi_1 = variables["First " + self.domain.lower() + " particle basis function"]
        xi_2 = variables["Second " + self.domain.lower() + " particle basis function"]
        T_av = variables["X-averaged cell temperature"]

        if self.domain == "Negative":
            D = self.param.D_n(c_s_av, T_av)
            C = self.param.C_n

        elif self.domain == "Positive":
            D = self.param.D_p(c_s_av, T_av)
            C = self.param.C_p

        self.rhs = {
            xi_1: -(1 / C) * pybamm.div(-D * pybamm.grad(xi_1)),
            xi_2: -(1 / C) * pybamm.div(-D * pybamm.grad(xi_2)),
        }

    def set_boundary_conditions(self, variables):

        xi_1 = variables["First " + self.domain.lower() + " particle basis function"]
        xi_2 = variables["Second " + self.domain.lower() + " particle basis function"]

        # define a(t), alpha_1(t), and alpha_2(t)

        # just writing out explicit expressions as not clear what
        # getting otherwise (interface submodel is a mess...)

        I = variables["Current collector current density"]
        l_p = self.param.l_p
        l_n = self.param.l_n

        if self.domain == "Negative":
            j = I / l_n
        elif self.domain == "Positive":
            j = -I / l_p

        c_s_surf_av = variables[
            "X-averaged " + self.domain.lower() + " particle surface concentration"
        ]
        T_av = variables["X-averaged " + self.domain.lower + " cell temperature"]

        if self.domain == "Negative":
            prefactor = self.param.m_n(T_av) / self.param.C_r_n

        elif self.domain == "Positive":
            prefactor = self.param.gamma_p * self.param.m_p(T_av) / self.param.C_r_p

        j0 = prefactor * (c_s_surf_av ** (1 / 2) * (1 - c_s_surf_av) ** (1 / 2))

        if self.domain == "Negative":
            U = self.param.U_n(c_s_surf_av, T_av)
        elif self.domain == "Positive":
            U = self.param.U_p(c_s_surf_av, T_av)

        U_prime = U.diff(c_s_surf_av)

        a_k = (
            1
            / 2
            * (
                j / c_s_surf_av
                + j / (1 - c_s_surf_av)
                - U_prime * ((j0 ** 2 + j ** 2) ** (0.5))
            )
        )
        alpha_1 = j * I / 2
        alpha_2 = I / 2 * ((j0 ** 2 + j ** 2) ** (0.5))

        self.boundary_conditions = {
            xi_1: {"left": 0, "right": a_k * pybamm.surf(xi_1) + alpha_1},
            xi_2: {"left": 0, "right": a_k * pybamm.surf(xi_2) + alpha_2},
        }

    def get_initial_conditions(self, variables):

        xi_1 = variables["First " + self.domain.lower() + " particle basis function"]
        xi_2 = variables["Second " + self.domain.lower() + " particle basis function"]

        self.initial_conditions = {xi_1: pybamm.Scalar(0), xi_2: pybamm.Scalar(0)}


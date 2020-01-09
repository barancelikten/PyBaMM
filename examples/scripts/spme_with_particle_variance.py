import pybamm
import numpy as np

t_eval = np.linspace(0, 0.16, 200)

c_rate = 3

options = {"particle variance": "Fickian"}
spme = pybamm.lithium_ion.SPMe(options)
spme.name = "With variance"
spme_sim = pybamm.Simulation(spme, C_rate=c_rate)
spme_sim.solve(t_eval=t_eval)


options = {"particle variance": None}
spme_2 = pybamm.lithium_ion.SPMe(options)
spme_2.name = "Without variance"
spme_2_sim = pybamm.Simulation(spme_2, C_rate=c_rate)
spme_2_sim.solve(t_eval=t_eval)

dfn = pybamm.lithium_ion.DFN()
dfn_sim = pybamm.Simulation(dfn, C_rate=c_rate)
dfn_sim.solve(t_eval=t_eval)

plot = pybamm.QuickPlot(
    [dfn_sim.built_model, spme_sim.built_model, spme_2_sim.built_model],
    dfn_sim.mesh,
    [dfn_sim.solution, spme_sim.solution, spme_2_sim.solution],
    output_variables=[
        "Negative particle concentration variance",
        "Positive particle concentration variance",
        "Negative particle surface concentration difference",
        "Positive particle surface concentration difference",
        "X-averaged open circuit voltage",
        "Corrector term",
        "Terminal voltage",
        "Terminal voltage [V]",
        "Electrolyte concentration",
    ],
)
plot.dynamic_plot()


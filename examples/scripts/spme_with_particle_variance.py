import pybamm
import numpy as np

options = {"particle variance": "Fickian"}
# options = {"particle variance": None}
model = pybamm.lithium_ion.SPMe(options)

sim = pybamm.Simulation(model)

t_eval = np.linspace(0, 0.14, 100)
sim.solve(t_eval=t_eval)
sim.plot()


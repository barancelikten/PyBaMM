import pybamm

options = {"particle variance": "Fickian"}
model = pybamm.lithium_ion.SPMe(options)

sim = pybamm.Simulation(model)
sim.solve()
sim.plot()

import pybamm

model = pybamm.lithium_ion.DFN()

sim = pybamm.Simulation(model)
sim.solve()
sim.plot()

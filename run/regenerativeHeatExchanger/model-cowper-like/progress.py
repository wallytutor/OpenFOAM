# -*- coding: utf-8 -*-
import os
import heat_recovery.cowper_like as cl
import matplotlib.pyplot as plt

plt.ion()
os.chdir("model-cowper-like")

origin = "0.000000e+00"

p = cl.plot_temperature("Initialization", origin, loc=3)

p = cl.plot_pressure("Initialization", origin)

p = cl.plot_flowrate("Initialization", origin, loc=3)

p = cl.plot_table("solid", "solidTemperature", "volFieldValue", time=origin)

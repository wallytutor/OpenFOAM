# -*- coding: utf-8 -*-
import os
from pathlib import Path
from heat_recovery.cowper_like import PostProcessing

try:
    HERE = Path(__file__).parent
except NameError:
    HERE = Path(os.getcwd()) / "01-fluid-only"
    os.chdir(HERE)

post = PostProcessing(HERE, mode="charging", domain=None)

post.plot_inlet_mass_flow_rate()
post.plot_imbalance_mass_flow_rate()
post.plot_total_pressure()
post.plot_pressure_drop()
post.plot_total_temperature()
post.plot_y_plus()
post.plot_convergence()

post.plot_field_temperature()
post.plot_field_pressure()
post.plot_field_density()
post.plot_field_velocity()
# post.plot_field_turb_kinetic_energy()
# post.plot_field_turb_dissipation_rate()
post.plot_field_buoyancy_pressure()

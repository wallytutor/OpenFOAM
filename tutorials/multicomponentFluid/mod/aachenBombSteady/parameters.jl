# -*- coding: utf-8 -*-
""" Setup for total combustion of C7H16 as C7H16 + 11O2 = 7CO2 + 8H2O. """

##############################################################################
# USER DEFINED
##############################################################################

mdot_fuel = 0.01
P_air = 101325.0
T_air = 800.0

##############################################################################
# CONSTANTS
##############################################################################

# Molecular masses of elements
const mH = 0.001
const mC = 0.012
const mO = 0.016

# Atomic composition of fuel C7H16
const nC = 7
const nH = 16

# Mass fraction of O2 in air
const YO2_air = 0.234

##############################################################################
# BALANCES
##############################################################################

# Molecular mass of fuel C7H16
mfuel = nC * mC + nH * mH
moxid = 2mO

# Mass of O2 per mole of C7H16
m_oxid_per_mole_fuel = 11moxid

# Mass of air per mole of C7H16 
m_air_per_mole_fuel = m_oxid_per_mole_fuel / YO2_air

# Mass flow rate [kg/s]
mdot_air = (mdot_fuel / mfuel) * m_air_per_mole_fuel

# Air reference density [kg/m³]
rho_air = (P_air * 0.02896) / (8.31446261815324 * T_air)

# Mean air inlet velocity [m/s]
U_air = mdot_air / (rho_air * 0.2^2)

##############################################################################
# DUMP RESULTS
##############################################################################

open("parameters", "w") do fp
    write(fp, """\
    YN2_air       $(1.0 - YO2_air);\n
    YO2_air       $(YO2_air);\n
    mdot_fuel     $(mdot_fuel);\n
    mdot_air      $(mdot_air);\n
    rho_air       $(rho_air);\n
    U_air         $(U_air);\n
    T_air         $(T_air);\n
    P_air         $(P_air);
    """)
end
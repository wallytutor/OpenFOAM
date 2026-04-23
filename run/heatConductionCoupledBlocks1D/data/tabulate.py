# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def conductivity(T):
    return 1.9847 - T * (2.531209E-03 - 1.43E-06 * T)

T_new = np.arange(200, 1301, 100)

df = pd.read_csv("cp.dat", sep=r"\s+", header=None)
cp_new = np.interp(T_new, df[0].to_numpy(),  df[1].to_numpy())

# only rhoConst is supported, not worth printint...
# df = pd.read_csv("rho.dat", sep=r"\s+", header=None)
# rho_new = np.interp(T_new, df[0].to_numpy(),  df[1].to_numpy())

kappa_new = conductivity(T_new)

for T, cp in zip(T_new, cp_new):
    print(f"({T:6.1f} {cp:6.1f})")

for T, kappa in zip(T_new, kappa_new):
    print(f"({T:6.1f} {kappa:6.1f})")

import numpy as np

def implicit_surface(X, Y, Z, **kwargs):
    a = np.sin(X) * np.cos(Y)
    b = np.sin(Y) * np.cos(Z)
    c = np.sin(Z) * np.cos(X)
    return a + b + c

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


def PowerLawFunction(f0, a, E0, Operation="Null"):

    E = sp.Symbol('E')

    f = f0 * (E / E0) ** (-a)

    if "Diff" in Operation:
        f = sp.diff(f, E)
    elif "Int" in Operation:
        f = sp.integrate(f, E)

    print("Function is ", f)

    f = sp.lambdify(E, f, "numpy")

    return f


if __name__ == "__main__":

    F = PowerLawFunction(10**4.3, 0.11, 10)

    print(F(10))


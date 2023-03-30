from labtool_ex2 import Project
from sympy import exp, pi, sqrt, Abs, pi

# from sympy.physics.units.systems.si import elementary_charge, boltzmann_constant
from scipy.constants import elementary_charge, Boltzmann, h, m_e
import numpy as np

# from numpy.typing import NDArray
import pandas as pd
import matplotlib.pyplot as plt  # noqa
import os
from uncertainties import ufloat


def test_compton_protokoll():
    gm = {
        "Z": r"Z",
        "S": r"S",
        "phi": r"\phi",
        "E": r"E_\text{Char}",
    }
    gv = {
        "Z": r"1",
        "S": r"1",
        "phi": r"\si{\degree}",
        "E": r"\si{\kilo\eV}",
    }

    pd.set_option("display.max_columns", None)
    plt.rcParams["axes.axisbelow"] = True
    P = Project("Compton", global_variables=gv, global_mapping=gm, font=13)
    P.output_dir = "./"
    P.figure.set_size_inches((10, 4))
    ax: plt.Axes = P.figure.add_subplot()

    # Aufgabe 1
    file = "../data/compton.csv"
    filepath = os.path.join(os.path.dirname(__file__), file)
    P.load_data(filepath, loadnew=True)


if __name__ == "__main__":
    test_compton_protokoll()

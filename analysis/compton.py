from labtool_ex2 import Project
from sympy import exp, pi, sqrt, Abs, cos

# from sympy.physics.units.systems.si import elementary_charge, boltzmann_constant
from scipy.constants import elementary_charge, Boltzmann, h, m_e, speed_of_light
import numpy as np

# from numpy.typing import NDArray
import pandas as pd
import matplotlib.pyplot as plt  # noqa
import os
from uncertainties import ufloat


def degreeUncertainty(x):
    return 0.1


def test_compton_protokoll():
    gm = {
        "Z": r"Z",
        "S": r"S",
        "phi": r"\phi",
        "E": r"E_\text{Char}",
        "Er": r"E_\text{Ruhe}",
    }
    gv = {
        "Z": r"1",
        "S": r"1",
        "phi": r"\si{\degree}",
        "E": r"\si{\kilo\eV}",
        "Er": r"\si{\kilo\eV}",
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

    P.data.loc[:, "dphi"] = phi.data.apply(degreeUncertainty)

    E_0 = 17.44

    P.plot_data(
        ax,
        phi,
        E,
        label="Gemessene Daten",
        style="#f49004",
        errors=True,
    )

    E = E_0 / (1 + E_0 * (1 - cos(phi / 180 * pi)) / Er)

    vals = P.plot_fit(
        axes=ax,
        x=phi,
        y=E,
        eqn=E,
        style=r"#f49004",
        label=r"Gut",
        offset=[0, 20],
        use_all_known=False,
        guess={
            "Er": 511,
        },
        bounds=[
            {"name": "Er", "min": 0.01, "max": 1e3},
        ],
        add_fit_params=True,
        granularity=10000,
        # gof=True,
        scale_covar=True,
    )
    Er = ufloat(vals["Er"].value, vals["Er"].stderr)
    print(Er * 1e3 * elementary_charge / speed_of_light**2)

    P.figure.suptitle("Fit Spezifischer Elektronenmasse")
    P.figure.tight_layout()
    P.ax_legend_all(loc=4)
    ax = P.savefig("Klein-Nishina.pdf")

    # Aufgabe 2
    file = "../data/alleSpektren.csv"
    filepath = os.path.join(os.path.dirname(__file__), file)
    P.load_data(filepath, loadnew=True)

    for elem in ["Ti", "Fe", "Ni", "Cu", "Zn", "Zr", "Mo", "Ag", "Mag", "Ring"]:
        print(P.data[f"E{elem}"].values[P.data[f"N{elem}"].idxmax()])

    # print(P.data)


if __name__ == "__main__":
    test_compton_protokoll()

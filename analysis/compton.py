from labtool_ex2 import Project
from sympy import exp, pi, sqrt, Abs, cos

# from sympy.physics.units.systems.si import elementary_charge, boltzmann_constant
from scipy.constants import (
    elementary_charge,
    Boltzmann,
    h,
    m_e,
    speed_of_light,
    Rydberg,
)
import numpy as np
from scipy import signal
from scipy.signal import find_peaks, peak_prominences

# from numpy.typing import NDArray
import pandas as pd
import matplotlib.pyplot as plt  # noqa
import os
from uncertainties import ufloat

# pyright: reportUnboundVariable=false
# pyright: reportUndefinedVariable=false


def degreeUncertainty(x):
    return 0.1


def test_compton_protokoll():
    gm = {
        "Z": r"Z",
        "S": r"S",
        "phi": r"\varphi",
        "E": r"E_\text{S}",
        "Er": r"E_\text{Ruhe}",
        "M": r"M",
        "ERyd": r"\sqrt{\frac{\tilde{E}}{R_y}}",
    }
    gv = {
        "Z": r"1",
        "S": r"1",
        "phi": r"\si{\degree}",
        "E": r"\si{\kilo\eV}",
        "Er": r"\si{\kilo\eV}",
        "M": r"\si{\kg}",
        "ERyd": r"1",
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

    P.print_table(
        phi,
        E,
        name="werte_compton",
        inline_units=True,
    )

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
        label=r"",
        offset=[0, 20],
        use_all_known=False,
        guess={
            "Er": 553,
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

    phi = np.linspace(30, 150, 1000)
    ax.plot(
        phi,
        E_0 / (1 + E_0 * (1 - np.cos(phi / 180 * np.pi)) / 510.998950),
        label="Theoretische Kurve",
    )

    P.figure.suptitle(
        "Bestimmung der Elektronenruheenergie\n mittels dem Klein-Nishina Wirkungsquerschnitt"
    )
    P.figure.tight_layout()
    P.ax_legend_all(loc=5)
    ax = P.savefig("Klein-Nishina.pdf")

    # Aufgabe 2
    file = "../data/alleSpektren.csv"
    filepath = os.path.join(os.path.dirname(__file__), file)
    P.load_data(filepath, loadnew=True)

    # elements = ["Ti", "Fe", "Ni", "Cu", "Zn", "Zr", "Mo", "Ag", "Mag", "Ring"]
    elements = ["Ti", "Fe", "Ni", "Cu", "Zn", "Zr", "Mo", "Ag"]
    k = np.asarray(
        [P.data[f"E{elem}"].values[P.data[f"N{elem}"].idxmax()] for elem in elements]
    )

    alpha = pd.DataFrame()
    Za = [22, 26, 28, 29, 30, 40, 42, 47]
    Ma = np.asarray([47.867, 55.84, 58.693, 63.55, 65.4, 91.22, 95.95, 107.868])
    alpha["E"] = k * elementary_charge * 1e3
    alpha["dE"] = (k * 0 + 0.2) * elementary_charge * 1e3
    alpha["Z"] = Za
    alpha["dZ"] = 0
    alpha["M"] = Ma * 1.6605390666e-27
    alpha["dM"] = alpha["M"] * 0.00001
    P.data = alpha
    P.vload()

    P.print_table(
        E,
        Z,
        name="roentgenKalpha",
        inline_units=True,
    )

    ERyd = sqrt(4 * E / (3 * h * speed_of_light * Rydberg) * (1 + m_e / M))
    P.resolve(ERyd)

    S = Z - ERyd
    P.resolve(S)
    print(S.data)

    P.plot_data(
        ax,
        Z,
        ERyd,
        label=r"$K_\alpha$ Daten",
        style="#f49004",
        errors=True,
    )

    ERyd = Z - S
    vals = P.plot_fit(
        axes=ax,
        x=Z,
        y=ERyd,
        eqn=ERyd,
        style=r"#f49004",
        label=r"$K_\alpha$",
        offset=[0, 20],
        use_all_known=False,
        guess={
            "S": 1.0,
        },
        bounds=[
            {"name": "S", "min": 0.8, "max": 1.4},
        ],
        add_fit_params=True,
        granularity=10000,
        # gof=True,
        scale_covar=True,
    )
    beta = pd.DataFrame()
    Zb = [28, 29, 30, 40, 42, 47]
    Mb = np.asarray([58.693, 63.55, 65.4, 91.22, 95.95, 107.868])
    kb = np.asarray([8.2, 8.7, 9.5, 17.5, 19.4, 24.4])
    beta["E"] = kb * elementary_charge * 1e3
    beta["dE"] = (kb * 0 + 0.2) * elementary_charge * 1e3
    beta["Z"] = Zb
    beta["dZ"] = 0
    beta["M"] = Mb * 1.6605390666e-27
    beta["dM"] = beta["M"] * 0.00001
    P.data = beta
    P.vload()

    P.print_table(
        E,
        Z,
        name="roentgenKbeta",
        inline_units=True,
    )

    ERyd = sqrt(9 * E / (8 * h * speed_of_light * Rydberg) * (1 + m_e / M))

    P.resolve(ERyd)

    S = Z - ERyd
    P.resolve(S)
    print(P.data)

    P.plot_data(
        ax,
        Z,
        ERyd,
        label=r"$K_\beta$ Daten",
        style="#BC2C1A",
        errors=True,
    )

    ERyd = Z - S

    vals = P.plot_fit(
        axes=ax,
        x=Z,
        y=ERyd,
        eqn=ERyd,
        style=r"#BC2C1A",
        label=r"$K_\beta$",
        offset=[0, 40],
        use_all_known=False,
        guess={
            "S": 2.0,
        },
        bounds=[
            {"name": "S", "min": 0.8, "max": 4},
        ],
        add_fit_params=True,
        granularity=10000,
        # gof=True,
        scale_covar=True,
    )
    P.figure.suptitle(
        "Bestimmung der Abschirmungskonstante \n mittels dem Moseyleysches Gesetzt"
    )
    P.figure.tight_layout()
    P.ax_legend_all(loc=4)
    ax = P.savefig("moseley.pdf")


if __name__ == "__main__":
    test_compton_protokoll()

import numpy as np
import pandas as pd
from common import *

L=2
N=L*L

e, m, C, chi = readfile_num_vals("datasets/numeric_sol.txt")

Z = 2 * np.exp(8.0) + 12 + 2*np.exp(-8.0)
analytic_E = - 1.0 / Z * (16 * np.exp(8.0) - 16 * np.exp(-8.0))
analytic_e = analytic_E / N
analytic_E2 = 128.0 / Z * (np.exp(8.0) + np.exp(-8.0))
analytic_magnetisation = 8.0 / Z * (np.exp(8.0) + 2)
analytic_magnetisation2 = 32.0 / Z * (np.exp(8.0) + 1)

analytic_heat_capacity = 1.0 / N * (analytic_E2 - analytic_E * analytic_E)
analytic_susceptibility = 1.0 / N * (analytic_magnetisation2 - analytic_magnetisation*analytic_magnetisation)

analytic_magnetisation = 8.0 / Z * (np.exp(8.0) + 2) / N

e.append(analytic_e)
m.append(analytic_magnetisation)
C.append(analytic_heat_capacity)
chi.append(analytic_susceptibility)

e_error = []
m_error = []
C_error = []
chi_error = []

for i in range(3):
    e_error.append(np.abs( abs(e[i] - analytic_e) / analytic_e * 100))
    m_error.append( abs(m[i] - analytic_magnetisation) / analytic_magnetisation * 100)
    C_error.append( abs(C[i] - analytic_heat_capacity) / analytic_heat_capacity *100)
    chi_error.append( abs(chi[i] - analytic_susceptibility) / analytic_susceptibility *100)

e_error.append(0)
m_error.append(0)
C_error.append(0)
chi_error.append(0)

df = pd.DataFrame(
    {
        "N" : [5000, 50000, 500000, "Analytic"],
        "epsilon" : e,
        "epsilon error" : e_error,
        "m" : m,
        "m error" : m_error,
        "C_V" : C,
        "C_V error": C_error,
        "chi" : chi,
        "chi error" : chi_error
    })

print(df.to_latex(index=False, float_format="%.5f"))
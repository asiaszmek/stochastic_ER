import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
from scipy import stats

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

plt.rcParams["text.latex.preamble"]+=r"\usepackage{sfmath}"
plt.rcParams["text.latex.preamble"]+=r"\usepackage{siunitx}"
plt.rcParams["text.latex.preamble"]+=r"\DeclareSIUnit{\molar}{M}"
plt.rcParams["text.latex.preamble"]+=r"\DeclareSIUnit{\Molar}{M}"


colors = {"Ca_wave_normal_SERCA_aging": "tab:orange",
          "Ca_wave_no_RyR_simple_SERCA_SOCE": "tab:purple",
          "Ca_wave_RyR2CaM_simple_SERCA_SOCE": "tab:green",
          "Ca_wave_simple_SERCA_SOCE": "tab:blue",
          "Ca_wave_RyR2CaM_aging": "tab:olive",
}
labels = {"Ca_wave_normal_SERCA_aging": r"old age",
          "Ca_wave_no_RyR_simple_SERCA_SOCE": r"no RyR2",
          "Ca_wave_RyR2CaM_simple_SERCA_SOCE": r"ctrl",
          "Ca_wave_simple_SERCA_SOCE": r"RyR2 no CaM",
          "Ca_wave_RyR2CaM_aging": r"Old age with fully inhibited RyR2",
}
if __name__ == "__main__":
    fname = sys.argv[1]
    my_data = pd.read_csv(fname)
    conditions = ["Ca_wave_RyR2CaM_simple_SERCA_SOCE",
                  "Ca_wave_no_RyR_simple_SERCA_SOCE",
                  "Ca_wave_simple_SERCA_SOCE", "Ca_wave_normal_SERCA_aging",
                  "Ca_wave_RyR2CaM_aging"]
    models = set(my_data["model"])
    my_data["stim"] = np.log2(my_data["stim"])
    formula = "y ~ stim + C(condition)"
    for model in models:
        data = my_data[my_data["model"] == model]
        my_linear_model = sm.ols(formula, data).fit()
      
        print(model)
        print(my_linear_model.summary())

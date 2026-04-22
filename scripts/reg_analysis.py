import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
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
    fname = "spatial_extent_4_cond.csv"
    my_data = pd.read_csv(fname)
    conditions = ["Ca_wave_RyR2CaM_simple_SERCA_SOCE", "Ca_wave_no_RyR_simple_SERCA_SOCE",
                  "Ca_wave_simple_SERCA_SOCE", "Ca_wave_normal_SERCA_aging",
                  "Ca_wave_RyR2CaM_aging"]
    models = set(my_data["model"])
    my_data["stim"] = np.log2(my_data["stim"])
    for model in models:
        data = my_data[my_data["model"] == model]
        fig, ax = plt.subplots(1,1)
        for condition in conditions:
            data_cond = data[data["condition"]==condition]
            x_val = data_cond["stim"]
            ax.plot(x_val, data_cond["y"], linewidth=0, marker="o",
                    label=labels[condition], color=colors[condition])
            x_val = data_cond["stim"]
            x_val_ols = sm.add_constant(x_val)
            y_val = data_cond["y"]
            ra_mod = sm.OLS(endog=y_val, exog=x_val_ols)
            res = ra_mod.fit()
            print(model, labels[condition])
            print(res.summary())
            b = res.params["const"]
            a = res.params["stim"]
            #print(a, b, r_value, p_value, std_err)
            x = np.linspace(x_val.min(), x_val.max(), 100)
            y = a*x+b
            ax.plot(x, y, linewidth=1, color=colors[condition])
        ax.legend()
        ax.set_xlabel(r"$\mathrm{\log_2 Ca_{i}^{2+}}$ $\left(\log_2\unit{\micro\Molar}\right)$")
        ax.set_ylabel(r"Spatial spread ($\unit{\micro\metre}$)")
        ax.set_title(r"Linear regression dend diam = %4.2f $\unit{\micro\metre}$" % model)

        fig.savefig("spatial_spread_lin_regression_%4.2f.png" % model)
    #plt.show()

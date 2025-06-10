import pandas as pd
import numpy as np


fname_basal = "Spatial_specificity_data_%s.csv"
fname_err_basal = "Spatial_specificity_data_error_%s.csv"
diams = ["1.2", "2.4", "6.0"]

if __name__ == "__main__":
    for diam in diams:
        data = pd.read_csv(fname_basal % diam)
        data_err = pd.read_csv(fname_basal % diam)
        ctrl = data.loc[0]
        for i in range(1, 7):
            row = data.loc[i]
            err = data_err.loc[i]
            print(diam, row[0])
            print(np.average(row.to_numpy()[1:]-ctrl.to_numpy()[1:],
                  weights=err.to_numpy()[1:]))

        ctrl = data.loc[1]
        for i in range(2, 7):
            row = data.loc[i]
            err = data_err.loc[i]
            print(diam, row[0])
            print(np.average(row.to_numpy()[1:]-ctrl.to_numpy()[1:],
                  weights=err.to_numpy()[1:]))

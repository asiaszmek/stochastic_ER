import numpy as np
import matplotlib.pyplot as plt

fname_list = {
    "TG": [
        "model_one_comp_1.8_CaOut_TG_trial0_dend.txt",
        "model_one_comp_1.8_CaOut_TG_trial1_dend.txt",
        "model_one_comp_1.8_CaOut_TG_trial2_dend.txt",
        "model_one_comp_1.8_CaOut_TG_trial3_dend.txt"
    ],
    "no_TG": [
        "model_one_comp_1.8_CaOut_trial0_dend.txt",
        "model_one_comp_1.8_CaOut_trial1_dend.txt",
        "model_one_comp_1.8_CaOut_trial2_dend.txt",
        "model_one_comp_1.8_CaOut_trial3_dend.txt"
    ],
}
specie = "Fura2Ca"

if __name__ == "__main__":
    # read in data
    data_TG = []
    data_no_TG = []
    for fname in fname_list["TG"]:
        f = open(fname, "r")
        header = f.readline().split(" ")
        data = np.loadtxt(f)
        idx = header.index(specie)

        time = data[:, 0]
        data_TG.append(data[:, idx])

    fig, ax = plt.subplots(1, 1)
    ax.plot(time/1000, np.mean(data_TG, axis=0), label="TG")

    for fname in fname_list["no_TG"]:
        f = open(fname, "r")

        header = f.readline().split(" ")
        data = np.loadtxt(f)
        idx = header.index(specie)
        time = data[:, 0]

        data_no_TG.append(data[:, idx])

    ax.plot(time/1000, np.mean(data_no_TG, axis=0), label="no TG")
    ax.set_ylabel("Fura2 Ca")
    ax.set_xlabel("Time (s)")
    ax.legend()

    plt.show()
    

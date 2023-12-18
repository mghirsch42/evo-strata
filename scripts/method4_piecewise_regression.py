import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import piecewise_regression

# data_file = "P_curvifrons_LG19_XY_SNPs" # best n breakpoints = 2
data_file = "P_babaulti_LG19_XY_SNPs" # best n breakpoints = 3-4
run_breakpoint_test = False # True to run the model selection to determine best number of breakpoints
n_breakpoints = 4
save_figure = True

figure_save_file = "figures/method4/{}_{}.png".format(data_file, n_breakpoints)

with open("data/"+data_file+".csv", "r") as f:
    data = f.readlines()
    data = [l.strip().split(",") for l in data]

df = pd.DataFrame(data[1:], columns=data[0])
df["count"] = df["count"].astype(int)
# print(df)

if run_breakpoint_test:
    ms = piecewise_regression.ModelSelection(list(df.index), list(df["count"]), max_breakpoints=10, max_iterations=300)
    exit()

pw_fit = piecewise_regression.Fit(list(df.index), list(df["count"]), n_breakpoints=n_breakpoints)
print(pw_fit.summary())

plt.plot(np.asarray(df.index), np.asarray(df["count"]), color="blue",linewidth=1)
pw_fit.plot_fit(color="red", linewidth=2)
pw_fit.plot_breakpoints(color="black")
pw_fit.plot_breakpoint_confidence_intervals(color="black")
plt.title("{} piecewise regression with {} breakpoints".format(data_file, n_breakpoints))

if save_figure:
    plt.savefig(figure_save_file)
plt.show()
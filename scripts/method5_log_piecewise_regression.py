import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import piecewise_regression
import math

# data_file = "P_curvifrons_LG19_XY_SNPs" # best n breakpoints = 2
data_file = "P_babaulti_LG19_XY_SNPs" # best n breakpoints = 3-4
run_breakpoint_test = True # True to run the model selection to determine best number of breakpoints
n_breakpoints = 2
save_figure = True

figure_save_file = "figures/method5/{}_{}.png".format(data_file, n_breakpoints)

with open("data/"+data_file+".csv", "r") as f:
    data = f.readlines()
    data = [l.strip().split(",") for l in data]

df = pd.DataFrame(data[1:], columns=data[0])
df["count"] = df["count"].astype(int)
df["ln_count"] = df["count"]
df = df.replace({"ln_count": {0: 10**-10}})
df["ln_count"] = np.log(df["ln_count"])
# print(df[df["ln_count"] == min(df["ln_count"])])
# print(df["ln_count"])
# plt.scatter(np.asarray(df.index), np.asarray(df["ln_count"]), color="blue",s=7)
# plt.show()
# exit()

if run_breakpoint_test:
    ms = piecewise_regression.ModelSelection(list(df.index), list(df["ln_count"]), max_breakpoints=10, max_iterations=300)
    exit()

pw_fit = piecewise_regression.Fit(list(df.index), list(df["ln_count"]), n_breakpoints=n_breakpoints)
print(pw_fit.summary())

pw_results = pw_fit.get_results()

plt.scatter(np.asarray(df.index), np.asarray(df["count"]), color="blue", s=7)
x = np.asarray(df.index)
breakpoint_1 = pw_results["estimates"]["breakpoint1"]["estimate"]
const = pw_results["estimates"]["const"]["estimate"]
alpha_1 = pw_results["estimates"]["alpha1"]["estimate"]
beta_1 = pw_results["estimates"]["beta1"]["estimate"]
beta_2 = pw_results["estimates"]["beta2"]["estimate"]
breakpoint_2 = pw_results["estimates"]["breakpoint2"]["estimate"]
print(breakpoint_1, const, alpha_1, beta_1)
# y = np.exp(const + alpha_1*x + beta_1 * np.maximum(x - breakpoint_1, 0))
y = np.exp(const + alpha_1*x + beta_1 * np.maximum(x - breakpoint_1, 0)) + beta_2 * np.maximum(x - breakpoint_2, 0)
plt.plot(x, y, color="black", linewidth=2)
plt.vlines(x=breakpoint_1, ymin=min(df["count"]), ymax=max(df["count"]), color="black", linewidth=2)
plt.vlines(x=breakpoint_2, ymin=min(df["count"]), ymax=max(df["count"]), color="black", linewidth=2)
plt.title("{}\nexponential piecewise regression\nwith {} breakpoints".format(data_file, n_breakpoints))

if save_figure:
    plt.savefig(figure_save_file)
plt.show()
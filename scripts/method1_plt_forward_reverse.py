from matplotlib import pyplot as plt
import pandas as pd

# Plots the forward and reverse data from the merge-as-you-go method on the same plot
# Requires splits to be saved for each method

# species = "P_curvifrons_LG19_XY_SNPs"
species = "P_babaulti_LG19_XY_SNPs"
block_range = "8-12"
plot_type = "line" # line or scatter
save_figure = True

forward_data = "results/method1/{}_{}_forward.csv".format(species, block_range)
reverse_data = "results/method1/{}_{}_reversed.csv".format(species, block_range)

with open("data/"+species+".csv", "r") as f:
    data = f.readlines()
    data = [l.strip().split(",") for l in data]
df = pd.DataFrame(data[1:], columns=data[0])
df["count"] = df["count"].astype(int)

with open(forward_data, "r") as f:
    forward_splits = f.readline()
    forward_splits = forward_splits.split(",")
    forward_splits = list(map(int, forward_splits))

with open(reverse_data, "r") as f:
    reverse_splits = f.readline()
    reverse_splits = reverse_splits.split(",")
    reverse_splits = list(map(int, reverse_splits))

forward_bin_vals = [forward_splits.count(x) for x in range(len(df))]
reverse_bin_vals = [reverse_splits.count(x) for x in range(len(df))]
bins = [x for x in range(len(df))]

y_max = max(df["count"])
bin_max = max(forward_bin_vals)
bin_scale = y_max / bin_max
forward_bin_vals = [bin_scale * bin_val for bin_val in forward_bin_vals]
bin_max = max(reverse_bin_vals)
bin_scale = y_max / bin_max
reverse_bin_vals = [bin_scale * bin_val for bin_val in reverse_bin_vals]


if plot_type == "scatter":
    plt.scatter(x=df.index, y=df["count"], color="black", s=7)
if plot_type == "line":
    plt.plot(df.index, df["count"], color="black")
plt.bar(bins, forward_bin_vals, width=1, alpha=0.3, color="blue", label="forward")
plt.bar(bins, reverse_bin_vals, width=1, alpha=0.3, color="orange", label="reverse")
plt.legend()


plt.title("Merge-as-you-go forward and reverse\n{}_block_sizes_{},".format(
    species, block_range))
save_file = "figures/method1/{}_{}_{}.png".format(species, block_range, plot_type)

if save_figure:
    plt.savefig(save_file)
plt.show()
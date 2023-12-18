from scipy import stats
from scipy.interpolate import make_interp_spline
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from statsmodels.stats.multitest import multipletests

data_file = "P_curvifrons_LG19_XY_SNPs"
# data_file = "P_babaulti_LG19_XY_SNPs"
block_sizes = [8, 9, 10, 11, 12]
plot_type = "line" # line or scatter -- how to show the SNP data
add_means = True # True if you want the plot to have the means of the splits
save_splits = True # True if you want to save a file with a list of the split locations
save_figure = True # True to save the figure showing splits

if save_splits:
    split_save_file = "results/method3/{}_{}-{}.csv".format(data_file, block_sizes[0], block_sizes[-1])
if add_means:
    figure_save_file = "figures/method3/{}_{}-{}_{}_means.png".format(data_file, block_sizes[0], block_sizes[-1], plot_type)
else:
    figure_save_file = "figures/method3/{}_{}-{}_{}.png".format(data_file, block_sizes[0], block_sizes[-1], plot_type)

with open("data/"+data_file+".csv", "r") as f:
    data = f.readlines()
    data = [l.strip().split(",") for l in data]

df = pd.DataFrame(data[1:], columns=data[0])
df["count"] = df["count"].astype(int)
# print(df)

blocks = {} # {block size 1: {block start 1: [splits], ...}, ...}
for block_size in block_sizes:
    blocks_b = {} # dictionary {block start 1: [splits], ...}
    for start in range(block_size):
        block_starts = np.arange(start, len(df), block_size) # start location of each block
        if block_starts[0] != 0: block_starts = np.insert(block_starts, 0, 0)
        
        merge_cont = True
        while(merge_cont):
            p_raw = [] # index 0 = pvalue between blocks 0 and 1; index 1 = pvalue between blocks 1 and 2
            for i in range(len(block_starts)-1):
                x = df.iloc[block_starts[i]:block_starts[i+1]]["count"].astype(int)
                if i == len(block_starts)-2:
                    y = df.iloc[block_starts[i+1]:len(df)]["count"].astype(int)
                else:
                    y = df.iloc[block_starts[i+1]:block_starts[i+2]]["count"].astype(int)
                mwu_res = stats.mannwhitneyu(x, y, method="exact")
                p_raw.append(mwu_res.pvalue)
            fdr_p=multipletests(pvals=p_raw, alpha=0.05, method="fdr_bh")[1]
            if any(p > 0.05 for p in fdr_p):
                merge_loc = np.argmax(fdr_p)
                block_starts = np.delete(block_starts, merge_loc+1)
            else:
                merge_cont = False
        blocks_b[start] = [block_starts[i+1] for i in range(len(fdr_p)) if fdr_p[i] < 0.05 ]
    blocks[block_size] = blocks_b
print(blocks)

splits = []
for b in blocks.keys():
    for s in blocks[b].keys():
        splits = splits + blocks[b][s]

splits.sort()

if save_splits:
    with open(split_save_file, "w") as f:
        f.write(",".join(list(map(str, splits))))

bin_vals = [splits.count(x) for x in range(len(df))]
bins = [x for x in range(len(df))]

if add_means:
    split_sum_cutoff = 15
    split_means = []
    curr_start = 0
    curr_end = curr_start
    while curr_end < len(splits) - 1:
        while splits[curr_end] == splits[curr_end + 1] or splits[curr_end] + 1 == splits[curr_end + 1]:
            curr_end = curr_end + 1
            if curr_end + 1 >= len(splits):
                break
        if curr_end - curr_start >= split_sum_cutoff:
            curr_mean = np.mean([splits[x] for x in range(curr_start, curr_end)])
            split_means = split_means + [curr_mean]
        curr_start = curr_end + 1
        curr_end = curr_start
        if curr_end >= len(splits):
            continue
    split_means = [round(x) for x in split_means]
    block_means = []
    df_idx = 0
    for split_idx in range(len(split_means)):
        block_means = block_means + [np.mean(df[df_idx:split_means[split_idx]]["count"])]
        df_idx = split_means[split_idx]
    block_means = block_means + [np.mean(df[df_idx:len(df)]["count"])]
print(block_means)

y_max = max(df["count"])
bin_max = max(bin_vals)
bin_scale = y_max / bin_max
bin_vals = [bin_scale * bin_val for bin_val in bin_vals]

fig = plt.figure(figsize=(10,6))
if plot_type == "scatter":
    plt.scatter(x=df.index, y=df["count"], color="blue", s=7)
if plot_type == "line":
    plt.plot(np.asarray(df.index), np.asarray(df["count"]), color="blue",linewidth=1)
plt.bar(bins, bin_vals, width=1, alpha=.25, color="orange")

if add_means:
    plt.vlines(split_means, 0, y_max, color="black")
    plt.hlines(block_means[0], 0, split_means[0], color="black")
    for i in range(len(split_means)-1):
        plt.hlines(block_means[i+1], split_means[i], split_means[i+1], color="black")
    plt.hlines(block_means[-1], split_means[-1], len(df), color="black")

plt.title("{}_block_sizes_{}-{},\neach repeated started at each possible location\n(pairwise merge)".format(
    data_file, block_sizes[0], block_sizes[-1]))

if save_figure:
    plt.savefig(figure_save_file)
plt.show()

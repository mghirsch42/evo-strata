from scipy import stats
from scipy.interpolate import make_interp_spline
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

# data_file = "P_curvifrons_LG19_XY_SNPs"
data_file = "P_babaulti_LG19_XY_SNPs"
block_sizes = [8, 9, 10, 11, 12]
reverse = True # False to run merging from the begining, False to run from the end
plot_type = "scatter" # line or scatter -- how to show the SNP data
save_splits = True # True if you want to save a file with a list of the split locations
save_figure = True # True to save the figure showing splits

if reverse:
    split_save_file = "results/method1/{}_{}-{}_reversed.csv".format(data_file, block_sizes[0], block_sizes[-1])
else:
    split_save_file = "results/method1/{}_{}-{}_forward.csv".format(data_file, block_sizes[0], block_sizes[-1])
    
with open("data/"+data_file+".csv", "r") as f:
    data = f.readlines()
    data = [l.strip().split(",") for l in data]

df = pd.DataFrame(data[1:], columns=data[0])
df["count"] = df["count"].astype(int)
# print(df)

blocks = {}
for b in range(len(block_sizes)):
    block_size = block_sizes[b]
    blocks_b = {}
    blocks_bs = []
    for start in range(block_size):
        if reverse:
            curr_block = [len(df) - start - block_size, len(df) - 1 - start]
            while curr_block[1] > 0:
                next_block = [curr_block[0]-block_size, curr_block[0]-1]
                if next_block[0] < 0:
                    if next_block[1] > 0:
                        next_block[0] = 0
                    else:
                        blocks_bs.append(curr_block)
                        break
                x = df.iloc[curr_block[0]:curr_block[1]+1]["count"].astype(int)
                y = df.iloc[next_block[0]:next_block[1]+1]["count"].astype(int)
                mwu_res = stats.mannwhitneyu(x, y)
                if mwu_res.pvalue < 0.01:
                    if curr_block[1] != len(df) - 1: # don't include the first block 
                        blocks_bs.append(curr_block)
                    curr_block = next_block
                else:
                    curr_block[0] = next_block[0]
        else: 
            curr_block = [start,start+block_size-1]
            while curr_block[1] < len(df):
                next_block = [curr_block[1]+1, curr_block[1]+block_size]
                if next_block[1] > len(df):
                    if next_block[0] < len(df):
                        next_block[1] = len(df)
                    else:
                        blocks_bs.append(curr_block)
                        break
                x = df.iloc[curr_block[0]:curr_block[1]+1]["count"].astype(int)
                y = df.iloc[next_block[0]:next_block[1]+1]["count"].astype(int)
                mwu_res = stats.mannwhitneyu(x, y)
                if mwu_res.pvalue < 0.01:
                    if curr_block[0] != start:
                        blocks_bs.append(curr_block) # don't include first block
                    curr_block = next_block
                else:
                    curr_block[1] = next_block[1]
        blocks_b[start] = blocks_bs
    blocks[b] = blocks_b

splits = []
for b in blocks.keys():
    for s in blocks[b].keys():
        for block in blocks[b][s]:
            splits.append(block[0])

splits = [x for x in splits if x != 0 and x != 259]

splits.sort()

if save_splits:
    with open(split_save_file, "w") as f:
        f.write(",".join(list(map(str, splits))))

bin_vals = [splits.count(x) for x in range(len(df))]
bins = [x for x in range(len(df))]

y_max = max(df["count"])
bin_max = max(bin_vals)
bin_scale = y_max / bin_max
bin_vals = [bin_scale * bin_val for bin_val in bin_vals]

if plot_type == "scatter":
    plt.scatter(x=df.index, y=df["count"], color="blue", s=7)
if plot_type == "line":
    plt.plot(df.index, df["count"], color="blue")
plt.bar(bins, bin_vals, width=1, alpha=.5, color="orange")

if reverse:
    plt.title("Merge-as-you-go\n{}_block_sizes_{}-{},\nreversed".format(
        data_file, block_sizes[0], block_sizes[-1]))
    save_file = "figures/method1/{}_{}-{}_{}_reversed.png".format(data_file, block_sizes[0], block_sizes[-1], plot_type)
else:
    plt.title("Merge-as-you-go\n{}_block_sizes_{}-{}\nforward".format(
        data_file, block_sizes[0], block_sizes[-1]))
    save_file = "figures/method1/{}_{}-{}_{}_forward.png".format(data_file, block_sizes[0], block_sizes[-1], plot_type)

if save_figure:
    plt.savefig(save_file)
plt.show()
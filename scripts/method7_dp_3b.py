import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from matplotlib import pyplot as plt
import time
from scipy.special import binom
import argparse
import os

def plot_results(x, reg_y, plt_y, bps, log, title, save_path):

    bps = np.concatenate([bps, [len(x)-1]])
    coef, ssr = peice_wise_regression(x, reg_y, bps)
    bps = np.concatenate([[0], bps])
    plt.figure()
    plt.scatter(x, plt_y, color="blue", s=5)
    for i in range(len(bps)-1):
        if i != 0: plt.vlines(x=bps[i], ymin=min(plt_y), ymax=max(plt_y), color="red", alpha=.5)
        if log:
            plt.plot(x[bps[i]:bps[i+1]+1], [np.exp(coef[i][0]*k+coef[i][1]) for k in x[bps[i]:bps[i+1]+1]], color="black")
        else:  
            plt.plot(x[bps[i]:bps[i+1]+1], [coef[i][0]*k+coef[i][1] for k in x[bps[i]:bps[i+1]+1]], color="black")
    plt.title(title + "\nSSR = " + f"{round(ssr):,}")
    plt.tight_layout()
    plt.savefig(save_path)
    # plt.show()

def regression(x, y):
    LR = LinearRegression().fit(x, y)
    pred = LR.predict(x)
    return round(((y-pred)**2).sum(), 3)

def peice_wise_regression(x, y, bps):
    ssr = 0
    params = []
    LR = LinearRegression()
    # 0 to first breakpoint
    reg = LR.fit(x[0:bps[0]+1].reshape(-1, 1), y[0:bps[0]+1])
    params.append((reg.coef_[0], reg.intercept_))
    pred = LR.predict(x[0:bps[0]+1].reshape(-1, 1))
    ssr += ((y[0:bps[0]+1] - pred)**2).sum()
    # Internal breakpoints
    for i in range(len(bps)-1):
        reg = LR.fit(x[bps[i]+1:bps[i+1]+1].reshape(-1, 1), y[bps[i]+1:bps[i+1]+1])
        params.append((reg.coef_[0], reg.intercept_))
        pred = LR.predict(x[bps[i]+1:bps[i+1]+1].reshape(-1, 1))
        ssr += ((y[bps[i]+1:bps[i+1]+1] - pred)**2).sum()
    return params, ssr

def precompute_ssr_mat(x, y, save_file):
    n = len(x)
    ssr_mat = np.zeros((n, n))
    ssr_mat[:] = float("inf")
    for start_point in range(0, n):
        for end_point in range(start_point+1, n):
            curr_x = x[start_point:end_point+1]
            curr_y = y[start_point:end_point+1]
            ssr = regression(curr_x.reshape(-1,1), curr_y)
            ssr_mat[start_point, end_point] = ssr
    if save_file:
        np.savetxt(save_file, ssr_mat, fmt="%1.3f")
    return ssr_mat

def load_ssr_mat(save_file):
    return np.loadtxt(save_file)

def base_case(ssr_mat, min_len, shift):
    n = len(ssr_mat)
    n_vals = n-2*min_len-shift+1
    curr_min = np.inf
    curr_idx = -1
    for i in range(0, n_vals):
        this_idx = shift+min_len-1+i
        this_min = ssr_mat[shift,shift+min_len-1+i] + ssr_mat[shift+min_len+i,n-1]
        if this_min < curr_min:
            curr_min = this_min
            curr_idx = this_idx
    return curr_min, [curr_idx]

def recurse(ssr_mat, n_bp, curr_n_bp, min_len, shift):
    if curr_n_bp == 1:
        return base_case(ssr_mat, min_len, shift)
    else:
        n = len(ssr_mat)
        curr_min = (np.inf)
        curr_idxs = []
        for bp1 in range(shift+min_len-1, n-(curr_n_bp*min_len)):
            bp1_val = ssr_mat[shift,bp1]
            return_min, return_idxs = recurse(ssr_mat, n_bp, curr_n_bp-1, min_len, bp1+1)
            if bp1_val+return_min < curr_min:
                curr_min = bp1_val+return_min
                curr_idxs = [bp1] + return_idxs
    return curr_min, curr_idxs

def get_seg_vals(ssr_mat, n_bp, min_len):
    n = len(ssr_mat)
    x = n-1-min_len-(n_bp*min_len-1)
    n_options = int(np.sum([binom(x-2+i, i)*(n_bp+1-i) for i in range(0, n_bp+1)]))
    # result_idxs = np.zeros((n_options, n_bp))
    # result_vals = np.zeros((n_options, n_bp+1))
    min_ssr, idxs = recurse(ssr_mat, n_bp, n_bp, min_len, 0)
    return min_ssr, idxs

def get_segs(ssr_mat, n_bp, min_len):
    val, idxs = get_seg_vals(ssr_mat, n_bp, min_len)
    return np.asarray(idxs), val

def main(data_path, data_name, log, n_bp, min_len, mat_save_path, result_save_path, fig_save_path):
    with open(data_path+data_name+".txt", "r") as f:
        data = f.readlines()
        data = [l.strip().split("\t") for l in data]
    df = pd.DataFrame(data, columns=["id","start","stop","count"])
    df["count"] = df["count"].astype(int)
    x = df.index.to_numpy()
    
    if not os.path.exists("{}{}/".format(result_save_path, data_name)):
        os.makedirs("{}{}/".format(result_save_path, data_name))
    if not os.path.exists("{}{}/ssr_mats/".format(result_save_path, data_name)):
        os.makedirs("{}{}/ssr_mats/".format(result_save_path, data_name))
    if not os.path.exists("{}{}/bps/".format(result_save_path, data_name)):
        os.makedirs("{}{}/bps/".format(result_save_path, data_name))
    if not os.path.exists("{}{}/".format(fig_save_path, data_name)):
        os.makedirs("{}{}/".format(fig_save_path, data_name))
    
    if log:
        df["ln_count"] = df["count"]
        df = df.replace({"ln_count": {0: 1}})
        df["ln_count"] = np.log(df["ln_count"])
        y = df["ln_count"].to_numpy()
        mat_save_file = "{}{}/ssr_mats/{}_log.txt".format(mat_save_path, data_name, data_name) 
        save_file = "{}{}/bps/{}_log_{}_{}.txt".format(result_save_path, data_name, data_name, n_bp, min_len)
        plot_title = "{}\n log, bp: {}, min len: {}".format(data_name, n_bp, min_len)
        fig_save_file = "{}{}/{}_log_{}_{}.png".format(fig_save_path, data_name, data_name, n_bp, min_len)
    else:
        y = df["count"].to_numpy()
        mat_save_file = "{}{}/ssr_mats/{}.txt".format(mat_save_path, data_name, data_name) 
        save_file = "{}{}/bps/{}_{}_{}.txt".format(result_save_path, data_name, data_name, n_bp, min_len)
        plot_title = "{}\nbp: {}, min len: {}".format(data_name, n_bp, min_len)
        fig_save_file = "{}{}/{}_{}_{}.png".format(fig_save_path, data_name, data_name, n_bp, min_len)

    if os.path.isfile(mat_save_file):
        print("Loading SSR matrix...")
        ssr_mat = load_ssr_mat(mat_save_file)
    else:
        print("Computing SSR matrix...")
        start_time = time.time()
        ssr_mat = precompute_ssr_mat(x, y, mat_save_file)
        end_time = time.time()
        print("Time:", end_time-start_time)
    
    print("Running breakpoint estimation...")
    start_time = time.time()
    bps, ssr = get_segs(ssr_mat, n_bp, min_len)
    end_time = time.time()
    print("Time:", end_time-start_time)
    print("Breakpoints:", bps)
    print("SSR:", [ssr])

    print("Saving to", save_file)
    with open(save_file, "w") as f:
        f.writelines("Breakpoints: " + str(bps))
        f.writelines("\n")
        f.writelines("SSR: " + str([ssr]))

    if log:
        plot_results(x, df["ln_count"], df["count"], bps, log, plot_title, fig_save_file)
    else:
        plot_results(x, df["count"], df["count"], bps, log, plot_title, fig_save_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--data_path", type=str, action="store", default="data/")
    parser.add_argument("-d", "--data_name", type=str, action="store", default=None)
    parser.add_argument("-l", "--log", action="store_true", default=False)
    parser.add_argument("-b", "--n_bp", type=int, action="store", default=1)
    parser.add_argument("-k", "--min_len", type=int, action="store", default=5)
    parser.add_argument("-m", "--mat_save_path", type=str, action="store", default="results/method7_dp/")
    parser.add_argument("-r", "--result_save_path", type=str, action="store", default="results/method7_dp_3b/")
    parser.add_argument("-f", "--fig_save_path", type=str, action="store", default="figures/method7_dp_3b/")
    args = parser.parse_args()
    if args.data_name == None:
        args.data_name = "P_curvifrons_LG19_XY_SNPs"
        # args.data_name = "P_babaulti_LG19_XY_SNPs"
    main(args.data_path, args.data_name, args.log, args.n_bp, args.min_len, args.mat_save_path, 
        args.result_save_path, args.fig_save_path)
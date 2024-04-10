import os
import argparse
import numpy as np
from sklearn.linear_model import LinearRegression
from scipy import stats
from matplotlib import pyplot as plt


def read_results(results_path, data_name, log, min_len, min_bp, max_bp):
    results = {}
    for bp in range(1, max_bp-min_bp+2):
        if log:
            path = "{}{}/bps/{}_log_{}_{}.txt".format(results_path, data_name, data_name, bp, min_len)
        else:
            path = "{}{}/bps/{}_{}_{}.txt".format(results_path, data_name, data_name, bp, min_len)
        with open(path, "r") as f:
            lines = f.readlines()
            bps = lines[0]
            ssrs = lines[1]
            bps = list(map(int, bps[14:-2].split()))
            ssrs = list(map(float, ssrs[6:-1].split()))
            results[bp] = {"bps": bps, "ssrs": ssrs}   
    return results

def read_data(data_path, data_name):
    with open("{}{}.txt".format(data_path, data_name), "r") as f:
        data = f.readlines()
        data = np.asarray([l.strip().split("\t") for l in data])
    y = data[:, 3].astype(float)
    x = np.linspace(0, len(y)-1, len(y))
    return x, y

def regression(x, y):
    LR = LinearRegression().fit(x, y)
    pred = LR.predict(x)
    return [y-pred]


def peice_wise_regression(x, y, bps):
    reds = []
    params = []
    LR = LinearRegression()
    # 0 to first breakpoint
    reg = LR.fit(x[0:bps[0]+1].reshape(-1, 1), y[0:bps[0]+1])
    params.append((reg.coef_[0], reg.intercept_))
    pred = LR.predict(x[0:bps[0]+1].reshape(-1, 1))
    reds = np.concatenate([reds, y[0:bps[0]+1] - pred])
    # Internal breakpoints
    for i in range(len(bps)-1):
        reg = LR.fit(x[bps[i]+1:bps[i+1]+1].reshape(-1, 1), y[bps[i]+1:bps[i+1]+1])
        params.append((reg.coef_[0], reg.intercept_))
        pred = LR.predict(x[bps[i]+1:bps[i+1]+1].reshape(-1, 1))
        reds = np.concatenate([reds, y[bps[i]+1:bps[i+1]+1] - pred])
    return params, reds

def t_tests(x, y, bps_ssrs):
    print("T-tests")
    bps = [bps_ssrs[i]["bps"] for i in range(1, len(bps_ssrs.keys())+1)]
    # print(bps)
    reds = regression(x.reshape(-1, 1), y)
    for i in range(len(bps)):
        curr_bps = np.concatenate([bps[i], [len(x)-1]])
        # print(i, curr_bps)
        _, curr_reds = peice_wise_regression(x, y, curr_bps)
        reds += [curr_reds.tolist()]
    tests = []
    for i in range(len(reds)-1):
        curr_test = stats.ttest_ind(np.abs(reds[i]), np.abs(reds[i+1]))
        print("{} vs {}: p={}".format(i, i+1, round(curr_test.pvalue, 4)))

def loglikelihoods(n, bps_ssrs):
    print("Log-Likelihoods")
    for i in bps_ssrs.keys():
        ssr = sum(bps_ssrs[i]["ssrs"])
        ll = -1*n/2*np.log(2*np.pi*ssr/n) - 1/(2*ssr/n)*ssr
        print("-ln(L|bp={}) = {}".format(i, round(ll, 4)))

def aics(n, bps_ssrs):
    print("Akaike Information Critiria")
    for i in bps_ssrs.keys():
        ssr = sum(bps_ssrs[i]["ssrs"])
        ll = -1*n/2*np.log(2*np.pi*ssr/n) - 1/(2*ssr/n)*ssr
        aic = 2*(2*i) - 2*ll
        print("AIC of {} breakpoints = {}".format(i, round(aic, 4)))

def main(data_path, data_name, results_path, log, min_len, min_bp, max_bp):
    bps_ssrs = read_results(results_path, data_name, log, min_len, min_bp, max_bp)
    x, y = read_data(data_path, data_name)
    t_tests(x, y, bps_ssrs)
    loglikelihoods(len(x), bps_ssrs)
    aics(len(x), bps_ssrs)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--data_path", type=str, action="store", default="data/")
    parser.add_argument("-d", "--data_name", type=str, action="store")
    parser.add_argument("-r", "--results_path", type=str, default="results/method7_dp_3b/")
    parser.add_argument("-l", "--log", action="store_true", default=False)
    parser.add_argument("-k", "--min_len", type=int, action="store", default=5)
    parser.add_argument("-i", "--min_bp", type=int, action="store", default=1)
    parser.add_argument("-j", "--max_bp", type=int, action="store", default=4)
    args = parser.parse_args()
    main(args.data_path, args.data_name, args.results_path, args.log,
         args.min_len, args.min_bp, args.max_bp)
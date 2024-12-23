#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <numeric>
#include <utility>
#include <algorithm>
#include <limits>
#include <math.h>
#include <Rcpp.h>
// #include "strata_finder.h"
using namespace std;
using namespace Rcpp;


Rcpp::NumericVector transform_y(Rcpp::NumericVector y, std::string fit_type) {
    Rcpp::NumericVector y_ (y.size());
    if (fit_type == "exponential") {
        transform(y.begin(), y.end(), y_.begin(), [](int i){
            if (i == 0) i = 1;
            return log(i);
        });
    }
    else if (fit_type == "linear") {
        y_ = Rcpp::NumericVector(y.begin(), y.end());
    }
    return y_;
};

std::vector<double> regression(std::vector<int> x, std::vector<double> y) {
    int n = x.size();
    const double sum_x = std::accumulate(x.begin(), x.end(), 0.0);
    const double sum_y = std::accumulate(y.begin(), y.end(), 0.0);
    const double sum_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const double sum_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    double m = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    double b = (sum_y * sum_xx - sum_x * sum_xy) / (n * sum_xx - sum_x * sum_x);
    double ssr = 0.0;
    for (int i=0; i<n; i++) {
        double ypred = m*x[i]+b;
        double sr = (y[i] - ypred) * (y[i] - ypred);
        ssr += sr;
    }
    std::vector<double> results {m, b, ssr};
    return results;
};

std::vector<double> get_residuals(Rcpp::IntegerVector x, Rcpp::NumericVector y, std::vector<double> ms, std::vector<double> bs, Rcpp::IntegerVector bps) {
    std::vector<double> preds = std::vector<double>();
    std::vector<double> reds = std::vector<double>();
    std::vector<double> y_sub;
    std::vector<double> x_sub;
    // bps.insert(bps.begin(), 0);
    // bps.push_back(this->x.size());
    // cout << "getting residuals " << bps.size() << endl;
    // cout << "breakpoints" << endl;
    for (int i=0; i<bps.size(); i++) {
        // cout << bps[i] << " ";
    }
    // cout << endl;

    // First segment goes from 0 up to and including the breakpoint.
    // cout << "first one " << 0 << " " << 0 << " " << bps[0]+1 << endl;
    // cout << "m0 " << ms[0] << endl;
    x_sub = std::vector<double> (x.begin(), x.begin()+bps[0]+1);
    std::transform(x_sub.begin(), x_sub.end(), back_inserter(preds), [ms, bs](int j){return ms[0]*j+bs[0];});
    y_sub = std::vector<double> (y.begin(), y.begin()+bps[0]+1);
    
    for(int j=0; j<x_sub.size(); j++) {
        reds.push_back(y_sub[j] - preds[j]);
    }
    // The rest start at the index after the breakpoint and include the next breakpoint.
    for (int i=0; i<bps.size()-1; i++) {
        preds.clear();
        x_sub = std::vector<double> (x.begin()+bps[i]+1, x.begin()+bps[i+1]+1);
        std::transform(x_sub.begin(), x_sub.end(), back_inserter(preds), [ms, bs, i](int j){return ms[i+1]*j+bs[i+1];});
        y_sub = std::vector<double> (y.begin()+bps[i]+1, y.begin()+bps[i+1]+1);
        for(int j=0; j<x_sub.size(); j++) {
            reds.push_back(y_sub[j] - preds[j]);
        }
    }
    return reds;
}

std::tuple<std::vector<int>, double> base_case(int min_len, int shift, std::vector<std::vector<double>> ssr_mat) {
    // cout << "In base case" << endl;
    int n = ssr_mat.size();
    int n_vals = n-2*min_len-shift+1;
    double curr_min = std::numeric_limits<int>::max();
    int curr_idx = -1;
    for (int i=0; i<n_vals; i++) {
        // cout << "i " << i << endl;
        int this_idx = shift+min_len-1+i;
        // double this_min = ssr_mat(shift, shift+min_len-1+i) + ssr_mat(shift+min_len+i, n-1);
        double this_min = ssr_mat[shift][shift+min_len-1+i] + ssr_mat[shift+min_len+i][n-1];
        if (this_min < curr_min) {
            curr_min = this_min;
            curr_idx = this_idx;
        }
    }
    return std::tuple<std::vector<int>, double> {std::vector<int> {curr_idx}, curr_min};
};

std::tuple<std::vector<int>, double> recurse(int curr_n_bp, int min_len, int shift, std::vector<std::vector<double>> ssr_mat) {
// std::tuple<std::vector<int>, double> strata_finder::recurse(int curr_n_bp, int min_len, int shift, Rcpp::NumericMatrix ssr_mat) {
    // cout << "ssr mat rows: " << ssr_mat.nrow() << " cols: " << ssr_mat.ncol()) <<endl;
    // cout << "ssr mat val 5, 5: " << ssr_mat(5, 5);
    
    if (curr_n_bp == 1) {
        // cout << "In recurse: going to base case" << endl;
        return base_case(min_len, shift, ssr_mat);
    }
    else {
        // cout << "In recurse: " << curr_n_bp << " bp " << shift << " shift" <<  endl;
        std::vector<double> curr_mins;
        int n = ssr_mat.size();
        double curr_min = std::numeric_limits<int>::max();
        std::vector<int> curr_idxs;
        for (int bp1=shift+min_len-1; bp1<n-(curr_n_bp*min_len); bp1++) {
            // cout << "In if statement; bp1: " << bp1 << endl;
            // double bp1_val = ssr_mat(shift, bp1);
            double bp1_val = ssr_mat[shift][bp1];
            std::tuple<std::vector<int>, double> recurse_return = recurse(curr_n_bp-1, min_len, bp1+1, ssr_mat);
            double return_min = std::get<1>(recurse_return);
            if (bp1_val + return_min < curr_min) {
                curr_min = bp1_val + return_min;
                std::vector<int> return_idxs = std::get<0>(recurse_return);
                curr_idxs = std::vector<int> {bp1};  
                curr_idxs.insert(curr_idxs.end(), return_idxs.begin(), return_idxs.end());            
            }     
        }
        return std::tuple<std::vector<int>, double> {curr_idxs, curr_min};
    }
};


// [[Rcpp::export]]
Rcpp::List piecewise_regression(Rcpp::IntegerVector x, Rcpp::NumericVector y, Rcpp::IntegerVector bps, std::string fit_type) {
    // A note on the ranges for each regression segment: the vector initializer includes the first value and excludes
    // the last value. So the range x.begin(), x.begin()+bps[0]+1 includes 0 up to and including the first breakpoint.
    // The range x.begin()+bps[i], x.begin()+bps[i+1]+1 includes breakpoint i up to and including the next breakpoint.
    std::vector<double> ssrs;
    std::vector<double> ms;
    std::vector<double> bs;
    bps.push_back(x.size()-1);
    y = transform_y(y, fit_type);
    std::vector<double> reg_results = regression(
        std::vector<int> (x.begin(), x.begin()+bps[0]+1),
        std::vector<double> (y.begin(), y.begin()+bps[0]+1)
    );
    // cout << reg_results[0] << endl;
    ms.push_back(reg_results[0]);
    bs.push_back(reg_results[1]);
    ssrs.push_back(reg_results[2]);

    for (int i=0; i<bps.size()-1; i++) {        
        reg_results = regression(
            std::vector<int> (x.begin()+bps[i]+1, x.begin()+bps[i+1]+1),
            std::vector<double> (y.begin()+bps[i]+1, y.begin()+bps[i+1]+1)
        );
        ms.push_back(reg_results[0]);
        bs.push_back(reg_results[1]);
        ssrs.push_back(reg_results[2]);
    }
    

    std::vector<double> residuals = get_residuals(x, y, ms, bs, bps);
    // std::vector<std::vector<double>> results = std::vector<std::vector<double>> {ms, bs, ssrs};
    return Rcpp::List::create(
        Named("m_vals") = ms,
        Named("b_vals") = bs,
        Named("ssr_vals") = ssrs,
        Named("residuals") = residuals
    );
};

// [[Rcpp::export]]
std::vector<std::vector<double>> precompute_ssr_mat(Rcpp::IntegerVector x, Rcpp::NumericVector y, std::string fit_type) {
// Rcpp::NumericMatrix strata_finder::precompute_ssr_mat(Rcpp::IntegerVector x, Rcpp::NumericVector y, string fit_type) {
    // check_data("precompute_ssr_mat");
    // cout << "In precompute ssr mat" << endl;
    int n = x.size();
    Rcpp::NumericVector y_ = transform_y(y, fit_type);
    
    // cout << "initializing ssr mat to maxes" << endl;
    std::vector<std::vector<double>> ssr_mat_ (n);
    // Rcpp::NumericMatrix ssr_mat_ (n, n);
    for (int i=0; i<n; i++) {
        // for (int j=0; j<n; j++) {
            // ssr_mat_(i, j) = std::numeric_limits<int>::max();
            ssr_mat_[i] = std::vector<double> (n, std::numeric_limits<int>::max());
        // }
    }
    // cout << "calculating ssr mat values" << endl;
    for (int start_point=0; start_point<n; start_point++) {
        for (int end_point=start_point+1; end_point<n; end_point++) {
            std::vector<int> curr_x = std::vector<int> (x.begin()+start_point, x.begin()+end_point+1);
            std::vector<double> curr_y = std::vector<double> (y_.begin()+start_point, y_.begin()+end_point+1);
            std::vector<double> reg_results = regression(curr_x, curr_y);
            // ssr_mat_(start_point,end_point) = reg_results[2];
            ssr_mat_[start_point][end_point] = reg_results[2];
        }
    }
    // cout << "returning ssr mat" << endl;
    return ssr_mat_;
};

// [[Rcpp::export]]
Rcpp::List get_breakpoints(Rcpp::IntegerVector x, Rcpp::NumericVector y, int n_bp, std::string fit_type, int min_len, std::vector<std::vector<double>> ssr_mat) {
    // cout << "In get_breakpoints" << endl;
    // check_data("get_breakpoints");
    if (ssr_mat.empty()) {
    // if (is_true(any(is_na(ssr_mat)))) {
        cout << "SSR matrix is not properly set. Call load_ssr or precompute_ssr before calling get_breakpoints." << endl;
        return Rcpp::List::create();
    }
    else if (n_bp < 1) {
        cout << "Number of breakpoints must be at least 1." << endl;
        return Rcpp::List::create();
    }
    else if (min_len < 1) {
        cout << "Minimum length must be at least 1." << endl;
        return Rcpp::List::create();
    }
    else {
        // cout << "Calling recurse first time" << endl;
        std::tuple<std::vector<int>, double> result = recurse(n_bp, min_len, 0, ssr_mat);
        return Rcpp::List::create(
            Named("bps") = std::get<0>(result),
            Named("ssr") = std::get<1>(result)
        );
    }
};


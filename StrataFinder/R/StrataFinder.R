StrataFinder <- setRefClass("StrataFinder",
    methods = list (
        findNStrata = function(x, y, n, fit_type, min_len, ssr_mat=NULL) {
            if (is.null(ssr_mat)) {
                ssr_mat <- precompute_ssr_mat(x, y, fit_type)
            }
            get_breakpoints(x, y, n, fit_type, min_len, ssr_mat)
        },
        findStrata = function(x, y, fit_type, min_len, min_strata, max_strata, ssr_mat=NULL) {
            if (is.null(ssr_mat)) {
                ssr_mat <- precompute_ssr_mat(x, y, fit_type)
            }
            if (min_strata == 0) {
                curr_residuals <- piecewise_regression(x, y, numeric(0), fit_type)
            }
            else {
                curr_results <- get_breakpoints(x, y, min_strata, fit_type, min_len, ssr_mat)
                curr_residuals <- piecewise_regression(x, y, curr_results$bps, fit_type)
            }
            for (n in min_strata+1:max_strata) {
                next_results <- get_breakpoints(x, y, n, fit_type, min_len, ssr_mat)
                print(paste("Breakpoints for", n, "breakpoints =", toString(next_results$bps)))
                next_residuals <- piecewise_regression(x, y, next_results$bps, fit_type)
                p <- ttest(curr_residuals$residuals, next_residuals$residuals)
                print(paste("P-value between", n-1, "and", n, "=", p))
                curr_residuals <- next_residuals
            }
        },
        precomputeSSRMat = function(x, y, fit_type) {
            precompute_ssr_mat(x, y, fit_type)
        },
        ttest = function(residuals1, residuals2) {
            t <- t.test(abs(residuals1), abs(residuals2), var.equal=TRUE)
            t$p.value
        }
    )
)
# Identifying evolutionary strata

## Data

The code expects the data files in a folder named "data" at the root of the directory.
Data files should be tab-deliminated, with the first column labeled as the linkage group, the second column as the first index of the block, the third as the end index of the block, and the fourth as the counts in that block. See examples currently in the "data/" folder.

## Method 7: Dynamic Programming

Breakpoint estimation algorithms are found in the scripts/ folder.

Method 7 uses a dynamic programming approach based on "Estimating and Testing Linear Models with Multiple Structural Changes" by Bai and Perron in 1998.

You can run Method 7 using scripts/method7_dp_3b.py, specifying the following flags:
-p/--data_path: the path to the folder containing the data files (default="data/");
-d/--data_name: the name of the data file (default=None);
-l/--log: binary flag indicating if you want to use exponential curves rather than linear fits (default=False);
-b/--n_bp: number of breakpoints (default=1);
-k/--min_len: minimum length of the segments (default=5);
-m/--mat_save_path: the path to save the calculated sum squared residuals (intermediate file) (default="results/method7_dp_3b/");
-r/--result_save_path: path to save the estimated breakpoints (default="results/method7_dp_3b/");
-f/--fig_save_path: path to save the figures (default="figures/method7_dp_3b/")

You can run ``python scripts/method7_dp_3b.py -h`` to see these options listed.

For example, to estimate breakpoint locations on the P_curvifrons_LG19_XY_SNPs data with 3 breakpoints, a minimum segment length of 10, linear fits, a matrix save path of "temp/", and other default parameters, you would run:
``python scripts/method7_dp_3b.py -d P_curvifrons_LG19_XY_SNPs -b 3 -k 10 -m temp/``

Or, to estimate breakpoint locations on the P_curvifrons_LG19_XY_SNPs data with exponential fits and otherwise default parameters, you would run:
``python scripts/method7_dp_3b.py -d P_curvifrons_LG19_XY_SNPs -l``

If you want to loop over all data files in a folder with a given set of parameters, you can use scripts/method7_dp_3b_run_all.sh. Modify the parameters as you wish, and then run the script using:
``./scripts/method7_dp_3b_run_all.sh``

To calculate statistics over the results, use scripts/method7_dp_stats.py. This script has the following flags:
-p/--data_path: the path to the folder containing the data files (default="data/");
-d/--data_name: the name of the data file (default=None);
-r/--results_path: the path to the folder containing the results from breakpoint estimation (default=results/method7_dp_3b/);
-l/--log: binary flag indicating if you want to use exponential curves rather than linear fits (default=False);
-k/--min_len: minimum length of the segments (default=5);
-i/--min_bp: the minimum number of breakpoints to use for statistics (default=1);
-j/--max_bp: the maximum number of breakpoints to use for statistics (default=4)

You can run ``python scripts/method7_dp_stats.py -h`` to see these options listed.

Note: the input parameters specify the results files to pull, so method 7 must have been run and saved for the parameter combination specified.

For example, to get the statistics for P_curvifrons_LG19_XY_SNPs with the default parameters, 
``python scripts/method7_dp_stats.py -d P_curvifrons_LG19_XY_SNPs``

Or for the exponential fits with a minimum segment length of 10 and only for 2-3 breakpoints:
``python scripts/method7_dp_stats.py -d P_curvifrons_LG19_XY_SNPs -l -k 10 -i 2 -j 3``

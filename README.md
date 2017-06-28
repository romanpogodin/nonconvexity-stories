# nonconvexity-stories
Research project of the [Yury Maximov's group](http://faculty.skoltech.ru/people/yurymaximov) at Skoltech. 

## Test scripts
1. `algorithms_convergence*` -- runs all algorithms and record objective function values. Main script is `*_setup.m`, it is used to avoid nested loops in parfor
2. `algorithms_rank_and_cut*` -- runs all test for rank and mean cut. The code is a bit odd due to inital problems with parfor
3. `run_gset_test` -- downloads a graph from Gset repository and finds its cut with SDP and a chosen method
## Functions
1. `solve_*` -- solvers for different approaches
2. Other functions are used in the solvers
## Bachelor thesis results
Roman Pogodin's data and a Jupyter notebook with some scripts for plotting. Contains maxcut and psd testы, and also a test of first 21 Gset problems (from https://web.stanford.edu/%7Eyyye/yyye/Gset/)

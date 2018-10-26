# nonconvexity-stories
Research project of the [Yury Maximov's group](http://faculty.skoltech.ru/people/yurymaximov) at Skoltech. 

Cite as:

Pogodin, R., Krechetov, M., & Maximov, Y. (2018). Efficient rank minimization to tighten semidefinite programming for unconstrained binary quadratic optimization. In 55th Annual Allerton Conference on Communication, Control, and Computing, Allerton 2017 (Vol. 2018–Janua, pp. 1153–1159). IEEE. https://doi.org/10.1109/ALLERTON.2017.8262867

Arxiv version: https://arxiv.org/abs/1708.01690

## Test scripts
1. `algorithms_convergence*` -- runs all algorithms and record objective function values. Main script is `*_setup.m`, it is used to avoid nested loops in parfor
2. `algorithms_rank_and_cut*` -- runs all test for rank and mean cut. The code is a bit odd due to inital problems with parfor
3. `run_gset_test` -- downloads a graph from Gset repository and finds its cut with SDP and a chosen method
## Functions
1. `solve_*` -- solvers for different approaches
2. Other functions are used in the solvers
## Bachelor thesis results
Roman Pogodin's data and a Jupyter notebook with some scripts for plotting. Contains maxcut and psd tests, and also a test of first 21 Gset problems (from https://web.stanford.edu/%7Eyyye/yyye/Gset/)

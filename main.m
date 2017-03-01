laplacian_matrix = get_laplacian('random', 0.5, 10);

num_cut_finder_trials = 10;
is_quiet = true;
is_cvx_quiet = true;
p = 0.05;
eps = 0.1;
num_iter = 10;
precision = 0.001;
methods = {'schatten'};
rank_tolerance = 1e-3;

return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
    p, eps, num_iter, precision, num_cut_finder_trials, ...
    is_quiet, is_cvx_quiet, rank_tolerance);
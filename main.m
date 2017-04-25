laplacian_matrix = get_laplacian('triangle');

num_cut_finder_trials = 10;
is_quiet = true;
is_cvx_quiet = true;
p = 0.01;
eps = 0.1;
num_iter = 5000;
precision = 0.001;
methods = {'langevin'};
rank_tolerance = 1e-3;

get_schatten_statistics(0.5, 30, 30, 'lang - 0.5 - 30 - norm - eta0.1 - ksi100 - best - 10iter.csv') 

%return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
%    p, eps, num_iter, precision, num_cut_finder_trials, ...
%    is_quiet, is_cvx_quiet, rank_tolerance);
%return_values_map.keys
%return_values_map.values
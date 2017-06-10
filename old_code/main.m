% laplacian_matrix = get_laplacian('random', 0.5, 50);
% 
% num_cut_finder_trials = 10;
% is_quiet = true;
% is_cvx_quiet = true;
% p = 0.01;
% eps = 0.1;
% num_iter = 5000;
% precision = 0.001;
% methods = {'langevin'};
% rank_tolerance = 1e-3;
% 
% get_schatten_statistics(0.5, 30, 30, 'lang - 0.5 - 40 - norm-poiss-2 - eta0.1 - ksi0 - best - 5000iter.csv') 

%return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
 %   p, eps, num_iter, precision, num_cut_finder_trials, ...
  %  is_quiet, is_cvx_quiet, rank_tolerance);
%return_values_map.keys
%return_values_map.values

tic;
get_schatten_statistics(0.5, 50, 100, './results/langevin_05_50_100iter_bestvals_ksi10_eta01', {'langevin'}, false, 100)
disp(toc);
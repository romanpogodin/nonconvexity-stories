function algorithms_gset_test(graph_numbers, method)
% A script to test algorithms on Gset graphs 
% from https://web.stanford.edu/%7Eyyye/yyye/Gset/G1
if nargin < 2
    method = 'singval';
elseif nargin < 1
    graph_numbers = [1];
end
    
    
%% Preliminaries
eps = 0.005;
q = 0.8;
p = 0.01;
num_iter = 100;
precision = 1e-6;
write_precision = 6;
num_cut_finder_trials = 10000;
    
rank_tol_one = 1e-4;
rank_tol_two = 1e-6; 

colnames = {'rank', 'max_cut', 'cut_mean', 'cut_std'};
results_sdp = zeros(1, 5);  
results_relaxed = zeros(1, 5);  

folder = 'gset_results/';

%% Solutions
for n_graph = graph_numbers
    general_name = strcat('gset', int2str(n_graph));
    
    %% Downloading
    disp(strcat('Downloading G', int2str(n_graph)));
    graph_laplacian = read_gset_laplacian(n_graph); 
    
    disp('Solving...');
    %% SDP
    disp('...SDP');
    [sdp_matrix, cut, sdp_optval, cut_optval] = ...
        solve_maxcut_sdp(graph_laplacian, num_cut_finder_trials, true);
    
    results_sdp(1, 1) = rank(sdp_matrix, rank_tol_one);
    results_sdp(1, 2) = rank(sdp_matrix, rank_tol_two);
    [~, results_sdp(1, 3), results_sdp(1, 4), results_sdp(1, 5)] = ...
        compute_cut_randomized(graph_laplacian, ...
        sdp_matrix, num_cut_finder_trials);
    
    %% Relaxation
    if strcmp(method, 'singval')
        disp('...singular values');
        [matrix, ~] = solve_maxcut_singval(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, q, eps, num_iter, precision, false, true);
    elseif strcmp(method, 'logdet')
        disp('...logdet');
        [matrix, ~] = solve_maxcut_logdet(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, eps, num_iter, precision, false, true);
    elseif strcmp(method, 'grad')
        disp('...singular values');
        [matrix, ~] = solve_maxcut_grad(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, p, eps, num_iter, precision, false, true);
    elseif strcmp(method, 'irls')
        disp('...singular values');
        [matrix, ~] = solve_maxcut_irls(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, p, eps, num_iter, precision, false, true); 
    end
        
    results_relaxed(1, 1) = rank(matrix, rank_tol_one);
    results_relaxed(1, 2) = rank(matrix, rank_tol_two);
    [~, results_relaxed(1, 3), results_relaxed(1, 4), results_relaxed(1, 5)] = ...
        compute_cut_randomized(graph_laplacian, ...
        matrix, num_cut_finder_trials);
    
    disp('...done');
    
    %% Writing results
    disp('Writing results...');
    
    dlmwrite(strcat(folder, 'sdp_', general_name), ...
        results_sdp, 'precision', write_precision); 
    
    if strcmp(method, 'irls') || strcmp(method, 'grad')
        method = strcat(method, '_p', int2str(100 * p));
    elseif strcmp(method, 'singval')
        method = strcat(method, '_q', int2str(100 * q));
    end
    
    dlmwrite(strcat(folder, method, '_', general_name), ...
        results_relaxed, ...
        'precision', write_precision); 
    
    disp('...done');
end
end

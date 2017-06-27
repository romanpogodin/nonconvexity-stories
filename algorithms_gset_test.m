% A script to test algorithms on Gset graphs 
% from https://web.stanford.edu/%7Eyyye/yyye/Gset/G1
%% Preliminaries
eps = 0.005;
q = 0.8;
num_iter = 100;
precision = 1e-6;
write_precision = 6;
num_cut_finder_trials = 1000;
    
rank_tol_one = 1e-4;
rank_tol_two = 1e-6; 

graph_numbers = [1];

colnames = {'rank', 'max_cut', 'cut_mean', 'cut_std'};
results_sdp = zeros(1, 5);  
results_singval = zeros(1, 5);  

folder = 'gset_results';

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
    
    %% Singular values
    disp('...singular values');
    [matrix, ~] = solve_maxcut_singval(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, q, eps, num_iter, precision, false, true);
        
    results_singval(1, 1) = rank(matrix, rank_tol_one);
    results_singval(1, 2) = rank(matrix, rank_tol_two);
    [~, results_singval(1, 3), results_singval(1, 4), results_singval(1, 5)] = ...
        compute_cut_randomized(graph_laplacian, ...
        matrix, num_cut_finder_trials);
    
    disp('...done');
    
    %% Writing results
    disp('Writing results...');
    
    dlmwrite(strcat(folder, 'sdp_', general_name), ...
        results_sdp, 'precision', write_precision); 
    
    dlmwrite(strcat(folder, 'singval_', general_name, ...
        '_q', int2str(10 * q)), ...
        results_singval, ...
        'precision', write_precision); 
    
    disp('...done');
end
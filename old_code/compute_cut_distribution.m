function compute_cut_distribution(...
    graph_size, graph_prob, n_cut_finder_trials)
num_iter = 100;
laplacian_matrix = get_laplacian('random', graph_prob, graph_size);

[sdp_matrix, ~, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(laplacian_matrix, 1, true);

[~, ~, logdet_matrix] = solve_maxcut_logdet(...
    laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, ...
    1.0, 0.001, num_iter, 0.0001, 1, true);

[~, ~, irls_matrix] = solve_maxcut_irls(...
    laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, ...
    0.001, 0.0001, num_iter, 0.0001, 1, true);

results = table(zeros(n_cut_finder_trials, 1), ...
    zeros(n_cut_finder_trials, 1), zeros(n_cut_finder_trials, 1), ...
    'VariableNames', {'SDP', 'LogDet', 'IRLS'});

cut_val = 0;

% normal_vector = randn(size(decomposed_matrix, 1), 1);

for i = 1:n_cut_finder_trials
    [~, cut_val] = compute_cut_randomized(laplacian_matrix, sdp_matrix, 1);
    results.SDP(i) = cut_val;
    
    [~, cut_val] = compute_cut_randomized(laplacian_matrix, logdet_matrix, 1);
    results.LogDet(i) = cut_val;
    
    [~, cut_val] = compute_cut_randomized(laplacian_matrix, irls_matrix, 1);
    results.IRLS(i) = cut_val;
end

writetable(results, strcat('./results/cut_res_logdet_irls_s', ...
    int2str(graph_size), '_p', int2str(100 * graph_prob), '.csv'));
end
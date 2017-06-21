%% Run the same matrix for each method
% Laplacian
% objective_matrix = get_laplacian('random', prob, graph_size); 

objective_matrix = randn(graph_size, graph_size);
% PSD
objective_matrix = objective_matrix * transpose(objective_matrix);
% indef
% objective_matrix = objective_matrix + transpose(objective_matrix);

[sdp_matrix, cut, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(objective_matrix, 10, true);

%% IRLS
if do_irls
    [irls_matrix, ~] = solve_maxcut_irls(...
        objective_matrix, sdp_optval, cut_optval, ...
        sdp_matrix, p, eps, num_iter, precision, false, true); 
end

%% Schatten grad
if do_grad
    [grad_matrix, ~] = solve_maxcut_grad(...
        objective_matrix, sdp_optval, cut_optval, ...
        sdp_matrix, p, eps, num_iter, precision, false, true);
end

%% Singular values grad
if do_singval
    [singval_matrix, ~] = solve_maxcut_singval(...
        objective_matrix, sdp_optval, cut_optval, ...
        sdp_matrix, q, eps, num_iter, precision, false, true);
end

%% Log-det
if do_logdet
    [logdet_matrix, ~] = solve_maxcut_logdet(...
        objective_matrix, sdp_optval, cut_optval, ...
        sdp_matrix, eps, num_iter, precision, false, true);
end

results = table(zeros(n_cut_finder_trials, 1), ...
    zeros(n_cut_finder_trials, 1), zeros(n_cut_finder_trials, 1), ...
    'VariableNames', {'sdp', 'logdet', strcat('irls_', int2str(10 * p)), ...
    strcat('grad_', int2str(10 * p)), strcat('singval', int2str(10 * q))});

cut_val = 0;

for i = 1:n_cut_finder_trials
    [~, cut_val] = compute_cut_randomized(objective_matrix, sdp_matrix, 1);
    results.SDP(i) = cut_val;
    
    [~, cut_val] = compute_cut_randomized(objective_matrix, logdet_matrix, 1);
    results.LogDet(i) = cut_val;
    
    [~, cut_val] = compute_cut_randomized(objective_matrix, irls_matrix, 1);
    results.IRLS(i) = cut_val;
end

if save_file
    writetable(results, strcat('./results/cut_res_sdp', ...
        int2str(graph_size), '_p', int2str(100 * graph_prob), '.csv'));
end
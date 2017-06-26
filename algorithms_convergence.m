% Saves results for mean convergence of proposed algorithms
% Tries different parameters and size. Requires Parallel Computing Toolbox

%% Preliminaries
disp(strcat('Running for size=', int2str(graph_size), ...
    ' and prob=', int2str(100 * prob)));
general_name = strcat('s', int2str(graph_size), '_pr', ...
    int2str(100 * prob));

% methods are hard-coded below due to parfor restrictions
if do_irls
    results_irls = zeros(num_trials, num_iter + 1);  
end

if do_grad
    results_grad = zeros(num_trials, num_iter + 1);  
end

if do_singval
    results_singval = zeros(num_trials, num_iter + 1);  
end
   
if do_logdet
    results_logdet = zeros(num_trials, num_iter + 1);  
end

%% (No) parallel simulations
for i = 1:num_trials 
    %% Run the same matrix for each method
    laplacian_matrix = get_laplacian('random', prob, graph_size);

    [sdp_matrix, cut, sdp_optval, cut_optval] = ...
        solve_maxcut_sdp(laplacian_matrix, 10, true);

    %% IRLS
    if do_irls
        [~, results_irls(i, :)] = solve_maxcut_irls(...
            laplacian_matrix, sdp_optval, cut_optval, ...
            sdp_matrix, p, eps, num_iter, precision, true, true); 
    end
    
    %% Schatten grad
    if do_grad
        [~, results_grad(i, :)] = solve_maxcut_grad(...
            laplacian_matrix, sdp_optval, cut_optval, ...
            sdp_matrix, p, eps, num_iter, precision, true, true);
    end
        
    %% Singular values grad
    if do_singval
        [~, results_singval(i, :)] = solve_maxcut_singval(...
            laplacian_matrix, sdp_optval, cut_optval, ...
            sdp_matrix, q, eps, num_iter, precision, true, true);
    end
        
    %% Log-det
    if do_logdet
        [~, results_logdet(i, :)] = solve_maxcut_logdet(...
            laplacian_matrix, sdp_optval, cut_optval, ...
            sdp_matrix, eps, num_iter, precision, true, true);
    end
    
    disp(strcat(int2str(i), ' done'));
end

disp('Write started');
%% Write results
if do_irls
    dlmwrite(strcat(folder, 'irls_', general_name, ...
        '_p', int2str(10 * p)), ...
        transpose(results_irls(:, :)), ...
        'precision', write_precision); 
end

if do_grad
    dlmwrite(strcat(folder, 'grad_', general_name, ...
        '_p', int2str(10 * p)), ...
        transpose(results_grad(:, :)), ...
        'precision', write_precision); 
end
   
if do_singval
    dlmwrite(strcat(folder, 'singval_', general_name, ...
        '_q', int2str(10 * q)), ...
        transpose(results_singval(:, :)), ...
        'precision', write_precision); 
end

if do_logdet
    dlmwrite(strcat(folder, 'logdet_', general_name), ...
        transpose(results_logdet), 'precision', write_precision); 
end

disp('Results written');        

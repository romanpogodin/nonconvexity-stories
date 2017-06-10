% Saves results for mean convergence of proposed algorithms
% Tries different parameters and size. Requires Parallel Computing Toolbox

sizes = [50, 100];
probabilities = [0.3, 0.5, 0.8];
p_norm = [0.1, 0.5, 0.8];
q_param = [0.5, 0.8];
eps = 0.001;
num_trials = 100;
num_iter = 100;
precision = 1e-8;
write_precision = 8;
folder = 'convergence_results/';

mkdir(folder);

for size = sizes
    for prob = probabilities
        disp(strcat('Running for size=', int2str(size), ...
            ' and prob=', int2str(100 * prob)));
        general_name = strcat('s', int2str(size), '_pr', ...
            int2str(100 * prob));
        parfor i = 1:num_trials           
            %% Run the same matrix for each method
            laplacian_matrix = get_laplacian('random', prob, size);
            
            [sdp_matrix, cut, sdp_optval, cut_optval] = ...
                solve_maxcut_sdp(laplacian_matrix, 10, true);
            
            %% IRLS
            for p = p_norm
               [~, optvals] = solve_maxcut_irls(...
                   laplacian_matrix, sdp_optval, cut_optval, ...
                   sdp_matrix, p, eps, num_iter, precision, true, true) 

               dlmwrite(strcat(folder, 'irls_', general_name, ...
                   '_p', int2str(10 * p)), transpose(optvals), ...
                   '-append','precision', write_precision); 
            end

            %% Schatten grad
            for p = p_norm
               [~, optvals] = solve_maxcut_grad(...
                   laplacian_matrix, sdp_optval, cut_optval, ...
                   sdp_matrix, p, eps, num_iter, precision, true, true)

               dlmwrite(strcat(folder, 'grad_', general_name, ...
                   '_p', int2str(10 * p)), transpose(optvals), ...
                   '-append','precision', write_precision); 
            end

            %% Singular values grad
            for q = q_param
               [~, optvals] = solve_maxcut_singval(...
                   laplacian_matrix, sdp_optval, cut_optval, ...
                   sdp_matrix, q, eps, num_iter, precision, true, true)

               dlmwrite(strcat(folder, 'singval_', general_name, ...
                   '_q', int2str(10 * q)), transpose(optvals), ...
                   '-append','precision', write_precision); 
            end

            %% Log-det
            [~, optvals] = solve_maxcut_logdet(...
                laplacian_matrix, sdp_optval, cut_optval, ...
                sdp_matrix, eps, num_iter, precision, true, true)

            dlmwrite(strcat(folder, 'logdet_', general_name), ...
                transpose(optvals), '-append','precision', write_precision);
        end
    end
end
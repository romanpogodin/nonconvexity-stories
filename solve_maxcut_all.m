function return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
    p, eps, num_iter, precision, num_cut_finder_trials, ...
    is_quiet, is_cvx_quiet, rank_tolerance, is_rankone_start)
if nargin < 11
   is_rankone_start = false; 
end
if nargin < 10
    rank_tolerance = 1e-4;
end

if nargin < 9
    is_cvx_quiet = true;
end

if nargin < 8
    is_quiet = false;
end

if nargin < 7
    num_cut_finder_trials = 10;
end

if nargin < 6
    precision = 0.001;
end
                                           
if nargin < 5
    num_iter = 10;
end

if nargin < 4
    eps = 0.1;
end

if nargin < 3
    p = 1.0;
end

if ~is_quiet
    disp('Solving SDP...')
end
[sdp_matrix, cut, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(laplacian_matrix, num_cut_finder_trials, is_cvx_quiet);

return_values_map = containers.Map(...
    {'sdp_optval', 'sdp_cut_optval', 'sdp_rank'}, ...
    [sdp_optval, cut_optval, rank(full(sdp_matrix), rank_tolerance)]);

if ismember('schatten', methods)
    if ~is_quiet
        disp('Solving Schatten...')
    end
    
    if is_rankone_start
       x_start = cut * transpose(cut); 
    else
       x_start = sdp_matrix;
    end
    
    [schatten_cut, schatten_cut_optval, schatten_matrix] = ...
        solve_maxcut_irls(laplacian_matrix, sdp_optval, cut_optval, x_start, ...
        p, eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);
    
    return_values_map = [return_values_map; containers.Map(...
        {'schatten_optval', 'schatten_cut_optval', 'schatten_rank'}, ...
        [(0.25 * trace(laplacian_matrix * schatten_matrix)), ...
        schatten_cut_optval, rank(full(schatten_matrix), rank_tolerance)])];
end

if ismember('grad', methods)
    if ~is_quiet
        disp('Solving grad...')
    end
    
    if is_rankone_start
       x_start = cut * transpose(cut); 
    else
       x_start = sdp_matrix;
    end
    
    [grad_cut, grad_cut_optval, grad_matrix] = ...
        solve_maxcut_grad(laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
        eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);
    
    return_values_map = [return_values_map; containers.Map(...
        {'grad_optval', 'grad_cut_optval', 'grad_rank'}, ...
        [(0.25 * trace(laplacian_matrix * grad_matrix)), ...
        grad_cut_optval, rank(full(grad_matrix), rank_tolerance)])];
end

if ismember('langevin', methods)
    if ~is_quiet
        disp('Solving langevin...')
    end
    
    if is_rankone_start
       x_start = cut * transpose(cut); 
    else
       x_start = sdp_matrix;
    end
    
    [langevin_cut, langevin_cut_optval, langevin_matrix] = ...
        solve_maxcut_langevin(laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
        eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);
    
    return_values_map = [return_values_map; containers.Map(...
        {'langevin_optval', 'langevin_cut_optval', 'langevin_rank'}, ...
        [(0.25 * trace(laplacian_matrix * langevin_matrix)), ...
        langevin_cut_optval, rank(full(langevin_matrix), rank_tolerance)])];
end

if ismember('langevin_vector', methods)
    if ~is_quiet
        disp('Solving vectorized langevin...')
    end
    
    if is_rankone_start
       x_start = cut * transpose(cut); 
    else
       x_start = sdp_matrix;
    end
    
    [langevin_cut, langevin_cut_optval, langevin_matrix] = ...
        solve_maxcut_sgld_vector(laplacian_matrix, sdp_optval, ...
        cut_optval, x_start, p, eps, num_iter, precision, ...
        num_cut_finder_trials);
    
    return_values_map = [return_values_map; containers.Map(...
        {'slgd_vec_optval', 'slgd_vec_cut_optval', 'slgd_vec_rank'}, ...
        [(0.25 * trace(laplacian_matrix * langevin_matrix)), ...
        langevin_cut_optval, rank(full(langevin_matrix), rank_tolerance)])];
end

if ismember('greedy', methods)
    if ~is_quiet
        disp('Solving greedy...')
    end
    [greedy_cut, greedy_cut_optval] = ...
        solve_maxcut_greedy(laplacian_matrix, cut);
   
    return_values_map = [return_values_map; containers.Map(...
        {'greedy_cut_optval'}, greedy_cut_optval)];
end


if ~is_quiet
    disp('SDP cut optval:')
    disp(cut_optval)
    
    if ismember('schatten', methods)
        disp('Schatten optval:')
        disp(schatten_cut_optval)
    end
    
    if ismember('grad', methods)
        disp('Grad optval:')
        disp(grad_cut_optval)
    end
    
    if ismember('greedy', methods)
        disp('Greedy optval:')
        disp(greedy_cut_optval)
    end
    
    % transpose(svd(sdp_matrix))
    % transpose(svd(schatten_matrix))
    % svd(grad_matrix)
    % 0.25 * trace(laplacian_matrix * schatten_matrix)
    % 0.25 * t
end
end
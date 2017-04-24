function [cut, new_cut_optval, curr_x] = solve_maxcut_sgld_vector(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
    eps, num_iter, precision, num_cut_finder_trials)
if nargin < 9
    num_cut_finder_trials = 10;
end

if nargin < 8
    precision = 0.001;
end
                                           
if nargin < 7
    num_iter = 10;
end

if nargin < 6
    eps = 0.1;
end

if nargin < 5
    p = 1.0;
end

curr_x = x_start;
temp_x = x_start;
best_x = x_start;
best_optval_x = x_start;
best_rank_x = x_start;

eta = 0.1;
ksi = 10;

diagonal_index = logical(eye(size(temp_x)));

for n = 1:num_iter
    w = triu(randn(size(curr_x)), 1);
    w = w + transpose(w);
    
    temp_x = curr_x - 2 * eta * compute_schatten_grad(curr_x, p, eps) + ...
        sqrt(2 * eta / ksi) * w;
    temp_x(diagonal_index) = 1;

    if norm_schatten(temp_x - curr_x, p, eps) < precision
        break
    end
    
    [~, is_negative] = cholcov(temp_x);
    optval_temp_x = trace(laplacian_matrix * temp_x) / 4;
    
    assert(is_negative || (~is_negative && optval_temp_x <= sdp_optval), ...
        'Temporary optval is bigger than SDP');
    
    if (~is_negative && cut_optval <= optval_temp_x)
        curr_x = temp_x;
    
        if (norm_schatten(curr_x) < norm_schatten(best_x))
            best_x = curr_x;
        end

        if (optval_temp_x > 4 * trace(laplacian_matrix * best_optval_x))
            best_optval_x = curr_x;
        end

        if (rank(curr_x) < rank(best_rank_x))
            best_rank_x = curr_x;
        end
    end
end
curr_x = best_x;
    
[cut, new_cut_optval] = compute_cut_randomized(laplacian_matrix, curr_x, ...
    num_cut_finder_trials);
[cut_best_optval, new_cut_optval_best_optval] = ...
    compute_cut_randomized(laplacian_matrix, best_optval_x, ...
    num_cut_finder_trials);
[cut_best_rank, new_cut_optval_best_rank] = ...
    compute_cut_randomized(laplacian_matrix, best_rank_x, ...
    num_cut_finder_trials);

if (new_cut_optval_best_optval > new_cut_optval)
    cut = cut_best_optval;
    new_cut_optval = new_cut_optval_best_optval;
end

if (new_cut_optval_best_rank > new_cut_optval)
    cut = cut_best_rank;
    new_cut_optval = new_cut_optval_best_rank;
end

end
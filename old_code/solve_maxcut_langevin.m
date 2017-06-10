function [cut, new_cut_optval, curr_x] = solve_maxcut_langevin(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
    eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet)
if nargin < 10
    is_cvx_quiet = true;
end

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
best_x = x_start;
best_optval_x = x_start;
best_rank_x = x_start;

eta = 0.1;
ksi = 10;

for n = 1:num_iter
    w = randn(size(curr_x));
    curr_x = curr_x - eta * compute_schatten_grad(curr_x, p, eps) + ...
        sqrt(2 * eta / ksi) * w;
    
    if is_cvx_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
        
    variable Y(size(curr_x)) symmetric
        minimize norm(Y - curr_x, 1)
        subject to
            Y >= 0
            diag(Y) == 1
            4 * cut_optval <= trace(laplacian_matrix * Y) <= 4 * sdp_optval
    cvx_end

    if norm_schatten(Y - curr_x, p, eps) < precision
        break
    end
    
    if (norm_schatten(Y) < norm_schatten(best_x))
        best_x = Y;
    end
    
    if (trace(laplacian_matrix * Y) > trace(laplacian_matrix * best_optval_x))
        best_optval_x = Y;
    end
    
    if (rank(Y) < rank(best_rank_x))
        best_rank_x = Y;
    end
    
    curr_x = Y;
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
function [cut, new_cut_optval, curr_x] = solve_maxcut_logdet(...
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

for n = 1:num_iter
    curr_weight = mpower(transpose(curr_x) * curr_x + ...
        eps * eye(size(curr_x, 1)), -1);
    
    if is_cvx_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
        
    variable Y(size(curr_x)) symmetric
        minimize trace(curr_weight * Y)
        subject to
            Y >= 0
            diag(Y) == 1
            4 * cut_optval <= trace(laplacian_matrix * Y) <= 4 * sdp_optval
    cvx_end

    if norm_schatten(Y - curr_x, p, eps) < precision
        break
    end
    
    curr_x = Y;
end
    
[cut, new_cut_optval] = compute_cut_randomized(laplacian_matrix, curr_x, ...
num_cut_finder_trials);
end
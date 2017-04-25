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

eta = 0.1;
ksi = 100;

for n = 1:num_iter
    
    Grad = 2*p * mpower(transpose(curr_x) * curr_x + ...
        eps * eye(size(curr_x, 1)), (p - 2.0) / 2.0) * curr_x
    w = randn(size(curr_x))
    
    for i=1:size(curr_x, 1)
        for j=i:size(curr_x, 1)
            if (i == j)
                Grad(i, j) = 0
                w(i, j) = 0
            else
                Grad(j, i) = Grad(i, j)
                w(j, i) = w(i, j)
            end
        end
    end
    
    Y = curr_x - eta*Grad + sqrt(2*eta/ksi)*w
    
    
    if norm_schatten(Y - curr_x, p, eps) < precision
        break
    end
    
    if (norm_schatten(Y) < norm_schatten(best_x))
        best_x = Y;
    end
    
    curr_x = Y;
end
curr_x = best_x;
    
[cut, new_cut_optval] = compute_cut_randomized(laplacian_matrix, curr_x, ...
    num_cut_finder_trials);
end
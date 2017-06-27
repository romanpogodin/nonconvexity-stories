function [curr_x, optvals] = solve_maxcut_grad(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
    eps, num_iter, precision, record_optvals, is_cvx_quiet)
% SOLVE_MAXCUT_GRAD Solves maxcut problem, using Schatten norms
%   [curr_x, optvals] = SOLVE_MAXCUT_GRAD(laplacian_matrix, sdp_optval,
%   cut_optval, x_start) to solve with default parameters and user-defined
%   laplacian, start point and upper/lower bounds on tr(LX).
%   Other parameters are:
%   p -- Schatten norm parameter
%   eps -- smoothing parameter
%   num_iter -- maximum number of gradient descent iterations
%   precision -- minimum Frobenius norm of gradient to stop
%   records_optvals -- wheter to record optvals at each iteration
%   is_cvx_quiet -- whether to suppress CVX output

%% Defalt arguments
if nargin < 10
    is_cvx_quiet = true;
end

if nargin < 9
    record_optvals = false;
end

if nargin < 8
    precision = 1e-4;
end
                                           
if nargin < 7
    num_iter = 100;
end

if nargin < 6
    eps = 0.01;
end

if nargin < 5
    p = 0.5;
end

%% Start points
curr_x = x_start;
diag_index = logical(eye(size(x_start)));

if record_optvals
   optvals = zeros(num_iter + 1, 1);
else
   optvals = zeros(1);
end

optvals(1) = norm_schatten(x_start, p, eps) ^ p; 
new_optval = optvals(1);
grad = curr_x;

%% Gradient descent
alpha = 0.3;
beta = 0.8;
n_backtracking_steps = 10;

for n = 1:num_iter
    % Backtracking line search, Boyd p. 465., direction == -gradient
    grad = 2 * compute_schatten_grad(curr_x, p, eps);
    step = 1.0;
    curr_norm = norm_schatten(curr_x, p, eps) ^ p;
    
    i = 1;
    while norm_schatten(curr_x - step * grad, p, eps, true, diag_index) ^ p > ...
            curr_norm - alpha * step * norm(grad, 'fro') ^ 2
        step = beta * step;
        i = i + 1;
        if i > n_backtracking_steps
            break;
        end
    end

    curr_x = curr_x - step * grad;
    curr_x(diag_index) = 1;
    
    if min(eig(curr_x)) < 0 || 4 * cut_optval > trace(laplacian_matrix * curr_x)
        %% Projection
        if is_cvx_quiet
            cvx_begin sdp quiet
        else
            cvx_begin sdp
        end

        variable Y(size(curr_x)) symmetric
            minimize norm(Y - curr_x, 'fro')
            subject to
                Y >= 0
                diag(Y) == 1
                4 * cut_optval <= trace(laplacian_matrix * Y) <= 4 * sdp_optval
        cvx_end
        curr_x = Y;
    end

    new_optval = norm_schatten(curr_x, p, eps) ^ p; 
    if record_optvals
       optvals(n + 1) = new_optval; 
    end
    
    if norm(grad, 'fro') < precision
        if record_optvals
            optvals(n + 1:num_iter + 1) = new_optval;
        end
        break
    end
end
end

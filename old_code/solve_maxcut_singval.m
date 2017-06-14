function [curr_x, optvals] = solve_maxcut_singval(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, q, ...
    eps, num_iter, precision, record_optvals, is_cvx_quiet)
% SOLVE_MAXCUT_SINGVAL Solves maxcut problem, using singular value approximation
%   [curr_x, optvals] = SOLVE_MAXCUT_SINGVAL(laplacian_matrix, sdp_optval,
%   cut_optval, x_start) to solve with default parameters and user-defined
%   laplacian, start point and upper/lower bounds on tr(LX).

%% Defalt arguments
if nargin < 10
    is_cvx_quiet = true;
end
solve_maxcut_singval
if nargin < 9
    record_optvals = false;
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
    q = 0.8;
end

%% Main code
curr_x = x_start;
diagonal_index = logical(eye(size(x_start)));

if record_optvals
   optvals = zeros(num_iter + 1, 1);
else
   optvals = zeros(1);
end

optvals(1) = norm_singval(x_start, q, eps); 
curr_optval = optvals(1);
new_optval = optvals(1);

for n = 1:num_iter
    curr_x = curr_x - 2 * 1/n * compute_singval_grad(curr_x, q, eps);
    curr_x(diagonal_index) = 1;
    
    [~, is_negative] = cholcov(curr_x);
    
    if is_negative || 4 * cut_optval > norm_singval(curr_x, q, eps)
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

    new_optval = norm_singval(curr_x, q, eps); 
    if record_optvals
       optvals(n + 1) = new_optval; 
    end
    
    if abs(curr_optval - new_optval) < precision
        if record_optvals
            optvals(n + 1:num_iter + 1) = new_optval;
        end
        break
    end
    curr_optval = new_optval;
end
end
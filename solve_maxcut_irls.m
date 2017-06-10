function [curr_x, optvals] = solve_maxcut_irls(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, p, ...
    eps, num_iter, precision, record_optvals, is_cvx_quiet)
% SOLVE_MAXCUT_IRLS Solves maxcut problem, using IRLS for Schatten norms
%   [curr_x, optvals] = SOLVE_MAXCUT_IRLS(laplacian_matrix, sdp_optval,
%   cut_optval, x_start) to solve with default parameters and user-defined
%   laplacian, start point and upper/lower bounds on tr(LX).
%   Other parameters are:
%   p -- Schatten norm parameter
%   eps -- smoothing parameter
%   num_iter -- maximum number of gradient descent iterations
%   precision -- minimum Frobenius norm of the difference between two
%   points to stop
%   records_optvals -- wheter to record optvals at each iteration
%   is_cvx_quiet -- whether to suppress CVX output

%% Default arguments
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
    p = 1.0;
end

%% Start points
curr_x = x_start;

if record_optvals
   optvals = zeros(num_iter + 1, 1);
else
   optvals = zeros(1);
end

eye_matrix = eye(size(curr_x, 1));
optvals(1) = norm_schatten(curr_x, p, eps) ^ p; 
new_optval = optvals(1);

%% Iterative procedure
for n = 1:num_iter
    curr_weight_chol = chol(mpower(transpose(curr_x) * curr_x + ...
        eps * eye_matrix, ...
        (p - 2.0) / 2.0), 'lower');
    
    %% CVX
    if is_cvx_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
        
    variable Y(size(curr_x)) symmetric
        minimize norm(Y * curr_weight_chol, 'fro')
        subject to
            Y >= 0
            diag(Y) == 1
            4 * cut_optval <= trace(laplacian_matrix * Y) <= 4 * sdp_optval
    cvx_end
    
    %% Check solution
    new_optval = norm_schatten(Y, p, eps) ^ p; 
    if record_optvals
       optvals(n + 1) = new_optval; 
    end
    
    if norm(Y - curr_x, 'fro') < precision
        if record_optvals
            optvals(n + 1:num_iter + 1) = new_optval;
        end
        break
    end
    
    curr_x = Y;
end
end
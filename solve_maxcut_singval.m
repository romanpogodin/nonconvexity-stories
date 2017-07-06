function [curr_x, optvals] = solve_maxcut_singval(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, q, ...
    eps, num_iter, precision, record_optvals, is_cvx_quiet, ...
    is_constraint_relaxed, check_psd)
% SOLVE_MAXCUT_SINGVAL Solves maxcut problem, using singular value approximation
%   [curr_x, optvals] = SOLVE_MAXCUT_SINGVAL(laplacian_matrix, sdp_optval,
%   cut_optval, x_start) to solve with default parameters and user-defined
%   laplacian, start point and upper/lower bounds on tr(LX).
%   Other parameters are:
%   q -- rank relaxation parameter for (1 + eps ^ q) function multiplier
%   eps -- smoothing parameter
%   num_iter -- maximum number of gradient descent iterations
%   precision -- minimum Frobenius norm of gradient to stop, also used to
%   threshold minimum singular value (must be >= -precision)
%   records_optvals -- wheter to record optvals at each iteration
%   is_cvx_quiet -- whether to suppress CVX output
%   is_constraint_relaxed -- wheter to use 4W <= Tr(LX) or 4SDP=Tr(LX)
%   check_psd -- whether to check PSD constraint or not. Small step
%   implies, that this operation is not needed. Setting to false must speed
%   up computations

%% Defalt arguments
if nargin < 12
    check_psd = true;
end

if nargin < 11
    is_constraint_relaxed = true;
end

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
    q = 0.8;
end

%% Start points
curr_x = x_start;
diag_index = logical(eye(size(x_start)));

if record_optvals
   optvals = zeros(num_iter + 1, 1);
else
   optvals = zeros(1);
end

optvals(1) = norm_singval(x_start, q, eps); 
new_optval = optvals(1);
grad = curr_x;

%% Gradient descent
alpha = 0.3;
beta = 0.8;
n_backtracking_steps = 10;

for n = 1:num_iter
    % Backtracking line search, Boyd p. 465., direction == -gradient
    grad = 2 * compute_singval_grad(curr_x, q, eps);
    step = eps / (4 * (1 + eps ^ q));
    curr_norm = norm_singval(curr_x, q, eps);
    
    i = 1;
    while norm_singval(curr_x - step * grad, q, eps, true, diag_index) > ...
            curr_norm - alpha * step * norm(grad, 'fro') ^ 2
        step = beta * step;
        i = i + 1;
        if i > n_backtracking_steps
            break;
        end
    end

    curr_x = curr_x - step * grad;
    curr_x(diag_index) = 1;
    
    if is_constraint_relaxed && ...
            4 * cut_optval > trace(laplacian_matrix * curr_x)
        disp('Reached value is smaller, than lower bound');
        break;
    end
    
    if (check_psd && min(eig(curr_x)) < -precision) || ...
            ~is_constraint_relaxed
        %% Projection
        disp(strcat('step ', int2str(n), ', projecting...'));
        curr_x = project_on_maxcut(curr_x, laplacian_matrix, ...
            cut_optval, sdp_optval, is_cvx_quiet, is_constraint_relaxed);
    end

    new_optval = norm_singval(curr_x, q, eps); 
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

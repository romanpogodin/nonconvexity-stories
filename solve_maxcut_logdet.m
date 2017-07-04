function [curr_x, optvals] = solve_maxcut_logdet(...
    laplacian_matrix, sdp_optval, cut_optval, x_start, ...
    eps, num_iter, precision, record_optvals, is_cvx_quiet, ...
    is_constraint_relaxed)
% SOLVE_MAXCUT_LOGDET Solves maxcut problem, using Log-Det heuristic
%   [curr_x, optvals] = SOLVE_MAXCUT_LOGDET(laplacian_matrix, sdp_optval,
%   cut_optval, x_start) to solve with default parameters and user-defined
%   laplacian, start point and upper/lower bounds on tr(LX).
%   Other parameters are:
%   eps -- smoothing parameter
%   num_iter -- maximum number of procedure iterations
%   precision -- minimum Frobenius norm of the difference between two
%   points to stop
%   records_optvals -- wheter to record optvals at each iteration
%   is_cvx_quiet -- whether to suppress CVX output
%   is_constraint_relaxed -- wheter to use 4W <= Tr(LX) or 4SDP=Tr(LX)

%% Defalt arguments
if nargin < 10
    is_constraint_relaxed = true;
end

if nargin < 9
    is_cvx_quiet = true;
end

if nargin < 8
    record_optvals = false;
end

if nargin < 7
    precision = 1e-4;
end
                                           
if nargin < 6
    num_iter = 100;
end

if nargin < 5
    eps = 0.01;
end

%% Start points
curr_x = x_start;

if record_optvals
   optvals = zeros(num_iter + 1, 1);
else
   optvals = zeros(1);
end

eye_matrix = eye(size(curr_x, 1));
optvals(1) = log(det(x_start + eps * eye_matrix)); 
new_optval = optvals(1);

%% Iterative procedure
for n = 1:num_iter
    curr_weight = mpower(transpose(curr_x) * curr_x + ...
        eps * eye_matrix, -1);
    
    %% CVX
    if is_constraint_relaxed
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
    else
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
                trace(laplacian_matrix * Y) == 4 * sdp_optval
        cvx_end
    end
        
    %% Check solution
    new_optval = log(det(Y + eps * eye_matrix)); 
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

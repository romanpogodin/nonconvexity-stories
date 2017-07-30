function [sdp_matrix, cut, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(laplacian_matrix, num_cut_finder_trials, ...
    is_quiet, solver)
% SOLVE_MAXCUT_SDP Solves maxcut problem, using SDP relaxation
%   [sdp_matrix, cut, sdp_optval, cut_optval] =
%   SOLVE_MAXCUT_SDP(laplacian_matrix) to solve with default parameters and
%   user-defined laplacian.
%   Other parameters are:
%   num_cut_finder_trials -- number of trials to find the optimal cut
%   is_quiet -- whether to suppress solver's output
%   solver -- 'cvx' or 'sdplr'
%   
%   CVX sovler for disciplined convex programming http://cvxr.com/cvx/
%
%   SDPLR impelemts Burer-Monteiro algorithm for low-rank semidefinite
%   programming. 
%   To install it, download SDPLR from http://sburer.github.io/projects.html
%   Documentation: http://sburer.github.io/files/SDPLR-1.03-beta-usrguide.pdf
%
%   On 64x platforms, sparse matrices might not work. To enable this, add 
%   '-largeArrayDims ' to the mexcmd variable in mexinstall.m, and
%   reinstall SDPLR

%% Defalt arguments
if nargin < 4
    solver = 'cvx';
end

if nargin < 3
    is_quiet = true;
end

if nargin < 2
    num_cut_finder_trials = 1000;
end

if strcmp('cvx', solver)
    %% CVX
    if is_quiet
        cvx_begin sdp quiet
    else
        cvx_begin sdp
    end
        variable X(size(laplacian_matrix)) symmetric
        maximize (0.25 * trace(laplacian_matrix * X))
        subject to
            X >= 0
            diag(X) == 1
    cvx_end

    sdp_matrix = X;
    sdp_optval = cvx_optval; 
    
elseif strcmp('sdplr', solver)
    %% SDPLR
    problem_size = size(laplacian_matrix, 1);
    
    first_index = 1:problem_size;
    second_index = zeros(problem_size, 1);
    unit_values = ones(problem_size, 1);
    for i=1:problem_size
       second_index(i) = i + (i - 1) * problem_size;
    end
    A = sparse(first_index, second_index, unit_values, ...
        problem_size, problem_size ^ 2);
    
    K = struct('s', problem_size);
    b = ones(problem_size, 1);
    c = reshape(laplacian_matrix, problem_size ^ 2, 1);
    
    if is_quiet
        pars = struct('printlevel', 0);
        sdp_matrix = sdplr(A, b, -c, K, pars);  % returns a vector
    else
        sdp_matrix = sdplr(A, b, -c, K);  % returns a vector
    end
        
    sdp_optval = 0.25 * transpose(c) * sdp_matrix;
    sdp_matrix = reshape(sdp_matrix, problem_size, problem_size); 
    
else
    error(...
        'Available solvers are cvx and sdplr. You entered %s', ...
        solver);
    
end

[cut, cut_optval] = compute_cut_randomized(laplacian_matrix, sdp_matrix, ...
                                           num_cut_finder_trials);
end

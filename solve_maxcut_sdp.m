function [sdp_matrix, cut, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(laplacian_matrix, num_cut_finder_trials, is_cvx_quiet)
% SOLVE_MAXCUT_SDP Solves maxcut problem, using SDP relaxation
%   [sdp_matrix, cut, sdp_optval, cut_optval] =
%   SOLVE_MAXCUT_SDP(laplacian_matrix) to solve with default parameters and
%   user-defined laplacian.
%   Other parameters are:
%   num_cut_finder_trials -- number of trials to find the optimal cut
%   is_cvx_quiet -- whether to suppress CVX output

%% Defalt arguments
if nargin < 3
    is_cvx_quiet = true;
end

if nargin < 2
    num_cut_finder_trials = 1000;
end

%% CVX
if is_cvx_quiet
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
[cut, cut_optval] = compute_cut_randomized(laplacian_matrix, sdp_matrix, ...
                                           num_cut_finder_trials);
end
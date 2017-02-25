function solve_maxcut_all(laplacian_matrix, methods, p, eps, num_iter, ...
    precision, num_cut_finder_trials, is_quiet, is_cvx_quiet)
if nargin < 9
    is_cvx_quiet = true;
end

if nargin < 8
    is_quiet = false;
end

if nargin < 7
    num_cut_finder_trials = 10;
end

if nargin < 6
    precision = 0.001;
end
                                           
if nargin < 5
    num_iter = 10;
end

if nargin < 4
    eps = 0.1;
end

if nargin < 3
    p = 1.0;
end

if ~is_quiet
    disp('Solving SDP...')
end
[sdp_matrix, cut, sdp_optval, cut_optval] = ...
    solve_maxcut_sdp(laplacian_matrix, num_cut_finder_trials, is_cvx_quiet);

if ismember('schatten', methods)
    if ~is_quiet
        disp('Solving Schatten...')
    end
    [schatten_cut, schatten_cut_optval, schatten_matrix] = ...
        solve_maxcut_irls(laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, p, ...
        eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);
end

if ismember('grad', methods)
    if ~is_quiet
        disp('Solving grad...')
    end
    [grad_cut, grad_cut_optval, grad_matrix] = ...
        solve_maxcut_grad(laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, p, ...
        eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);
end


if ~is_quiet
    disp('SDP cut optval:')
    disp(cut_optval)
    
    if ismember('schatten', methods)
        disp('Schatten optval:')
        disp(schatten_cut_optval)
    end
    
    if ismember('grad', methods)
        disp('Grad optval:')
        disp(grad_cut_optval)
    end
    % transpose(svd(sdp_matrix))
    % transpose(svd(schatten_matrix))
    % svd(grad_matrix)
    % 0.25 * trace(laplacian_matrix * schatten_matrix)
    % 0.25 * t
end

end
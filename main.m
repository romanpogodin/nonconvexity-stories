% adj_matrix = cat(2, cat(1, zeros(4), ones(4)), cat(1, ones(4), zeros(4))); %K_{4,4}, MaxCut = 16
% laplacian_matrix = full(laplacian(graph(adj_matrix)));
laplacian_matrix = [3, -1, 0, 0, -1, -1;
            -1, 3, -1, -1, 0, 0;
            0, -1, 3, -1, 0, -1;
            0, -1, -1, 3, -1, 0;
            -1, 0, 0, -1, 3, -1;
            -1, 0, -1, 0, -1, 3]; % topology 5
% laplacian_matrix = [3, -1, 0, 0, -1, -1;
%             -1, 3, -1, -1, 0, 0;
%             0, -1, 2 , -1, 0, 0;
%             0, -1, -1, 3, -1, 0;
%             -1, 0, 0, -1, 3, -1;
%             -1, 0, 0, 0, -1, 2]; % topology 4
% laplacian_matrix = [3, -1, -1, -1; -1, 3, -1, -1; -1, -1, 3, -1; -1, -1, -1, 3]; % topology 3
% laplacian_matrix = [2, -1, 0, -1; -1, 3, -1, -1; 0, -1, 2, -1; -1, -1, -1, 3]; % topology 2
% laplacian_matrix = [2, -1, -1, 0; -1, 3, -1, -1; -1, -1, 3, -1; 0, -1, -1, 2]; % topology 1
% laplacian_matrix = [2, -1, -1; -1, 2, -1; -1, -1 ,2]; % triangle


num_cut_finder_trials = 10;
is_cvx_quiet = true;
p = 0.5;
eps = 0.1;
num_iter = 10;
precision = 0.001;

[sdp_matrix, cut, sdp_optval, cut_optval] = ...
        solve_maxcut_sdp(laplacian_matrix, num_cut_finder_trials, is_cvx_quiet);

[schatten_cut, schatten_cut_optval, schatten_matrix] = ...
    solve_maxcut_irls(laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, p, ...
    eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);

[grad_cut, grad_cut_optval, grad_matrix] = ...
    solve_maxcut_grad(laplacian_matrix, sdp_optval, cut_optval, sdp_matrix, p, ...
    eps, num_iter, precision, num_cut_finder_trials, is_cvx_quiet);

sch_optvals = cat(2, sch_optvals, ...
    0.25 * trace(laplacian_matrix * schatten_matrix));
grad_optvals = cat(2, grad_optvals, ...
    0.25 * trace(laplacian_matrix * grad_matrix));

svd_vals = svd(schatten_matrix);
sch_last_singval = cat(2, sch_last_singval, ...
    svd_vals);

svd_vals = svd(grad_matrix);
grad_last_singval = cat(2, grad_last_singval, ...
    svd_vals);

% cut_optval
% schatten_cut_optval
% grad_cut_optval
% svd(sdp_matrix)
% svd(schatten_matrix)
% svd(grad_matrix)
% 0.25 * trace(laplacian * schatten_matrix)
% 0.25 * trace(laplacian * grad_matrix)
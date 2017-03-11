function [best_cut, best_optval] = solve_maxcut_greedy(...
    laplacian_matrix, cut)

best_cut = cut;
best_optval = 0.25 * transpose(cut) * laplacian_matrix * cut;

assert(size(cut, 1) == length(cut), 'cut must b a column vector');
all_cuts = repelem(cut, 1, length(cut)) .* ...
    (ones(length(cut)) - 2 * eye(length(cut))); 

[max_cuts, indices] = ...
    max(diag(0.25 * transpose(all_cuts) * laplacian_matrix * all_cuts));

if max_cuts(1) > 0.25 * transpose(cut) * laplacian_matrix * cut
    best_cut = all_cuts(:, indices(1));
    best_optval = max_cuts(1);
end
end
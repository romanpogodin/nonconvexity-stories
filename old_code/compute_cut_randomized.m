function [cut, cut_value] = compute_cut_randomized(laplacian_matrix, matrix, ...
                                                   num_trials, tolerance)

if nargin < 4
    tolerance = 0.001;
end
                                               
if nargin < 3
    num_trials = 10;
end

if min(eig(matrix)) < tolerance
    decomposed_matrix = chol(matrix + tolerance * eye(size(matrix, 1)));
else    
    decomposed_matrix = chol(matrix);
end
% decomposed_matrix = cholcov(matrix);
% decomposed_matrix = cat(1, decomposed_matrix, ...
%     zeros(size(decomposed_matrix, 2) - size(decomposed_matrix, 1), ...
%     size(decomposed_matrix, 2)));
cut_value = 0.0;

for n = 1:num_trials
    normal_vector = randn(size(decomposed_matrix, 1), 1);
    normal_vector = normal_vector ./ norm(normal_vector);
    curr_cut = 2.0 * ...
        double(transpose(decomposed_matrix) * normal_vector >= 0) - 1.0;
    
    if cut_value < 0.25 * transpose(curr_cut) * laplacian_matrix * curr_cut
        cut_value = 0.25 * transpose(curr_cut) * laplacian_matrix * curr_cut;
        cut = curr_cut;
    end
end
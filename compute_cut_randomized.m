function [cut, cut_value, cut_mean, cut_std] = compute_cut_randomized(...
    laplacian_matrix, matrix, num_trials, tolerance)
%% Default arguments
if nargin < 4
    tolerance = 1e-8;
end
                                               
if nargin < 3
    num_trials = 1000;
end

if min(eig(matrix)) < tolerance
    decomposed_matrix = chol(matrix + tolerance * eye(size(matrix, 1)));
else    
    decomposed_matrix = chol(matrix);
end

%% Cut finding
cut_value = 0.0;
all_cuts = zeros(num_trials, 1);

for n = 1:num_trials
    normal_vector = randn(size(decomposed_matrix, 1), 1);
    normal_vector = normal_vector ./ norm(normal_vector);
    curr_cut = 2.0 * ...
        double(transpose(decomposed_matrix) * normal_vector >= 0) - 1.0;
    
    all_cuts(n) =  0.25 * transpose(curr_cut) * laplacian_matrix * curr_cut;
    
    if cut_value < all_cuts(n)
        cut_value = all_cuts(n);
        cut = curr_cut;
    end
   
end

cut_mean = mean(all_cuts);
cut_std = std(all_cuts);
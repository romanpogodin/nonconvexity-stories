function run_gset_test(graph_numbers, method, is_constraint_relaxed, solver)
%   RUN_GSET_TEST A script to test algorithms on Gset graphs 
%   from https://web.stanford.edu/%7Eyyye/yyye/Gset/
%   RUN_GSET_TEST(graph_numbers) -- downloads Gset graphs, which numbers
%   specified by a vector graph_number
%   method -- method to improve an SDP solution: singval, grad, all (which is singval + grad)
%   is_constraint_relaxed -- wheter to use 4W <= Tr(LX) or 4SDP=Tr(LX)
%   solver -- which solver to use for SDP solution: cvx or sdplr

if nargin < 4
    solver = 'cvx';
end

if nargin < 3
    is_constraint_relaxed = true;
end

if nargin < 2
    method = 'singval';
end
    
if nargin < 1
    graph_numbers = [1];
end
    
    
%% Preliminaries
eps = 0.005;
q = 0.8;
p_param = [0.01, 0.1];
num_iter = 100;
precision = 1e-5;
write_precision = 5;
num_cut_finder_trials = 100000;
% enabled only for singval, might break everything for Schatten
check_psd = false;  

rank_tol_one = 1e-4;
rank_tol_two = 1e-5; 

% colnames = {'rank', 'max_cut', 'cut_mean', 'cut_std'};
results_sdp = zeros(1, 5);  
results_relaxed = zeros(1, 5);  

folder = 'gset_results/';
mkdir(folder);

%% Solutions
for n_graph = graph_numbers
    general_name = strcat('gset', int2str(n_graph), 'relaxed_', ...
        int2str(is_constraint_relaxed), '_', solver);
    
    %% Downloading
    disp(strcat('Downloading G', int2str(n_graph)));
    graph_laplacian = read_gset_laplacian(n_graph); 
    
    disp('Solving...');
    %% SDP
    disp('...SDP...');
    [sdp_matrix, ~, sdp_optval, cut_optval] = ...
        solve_maxcut_sdp(graph_laplacian, num_cut_finder_trials, true, solver);
    
    results_sdp(1, 1) = rank(sdp_matrix, rank_tol_one);
    results_sdp(1, 2) = rank(sdp_matrix, rank_tol_two);
    [~, results_sdp(1, 3), results_sdp(1, 4), results_sdp(1, 5)] = ...
        compute_cut_randomized(graph_laplacian, ...
        sdp_matrix, num_cut_finder_trials);
    
    disp('...writing results...');
    dlmwrite(strcat(folder, 'sdp_', general_name), ...
        results_sdp, 'precision', write_precision); 
    disp('...done');
    
    %% Relaxation    
    if strcmp(method, 'singval') || strcmp(method, 'all')
        %% Solving
        disp('...singular values...');
        [matrix, ~] = solve_maxcut_singval(...
            graph_laplacian, sdp_optval, cut_optval, ...
            sdp_matrix, q, eps, num_iter, precision, false, true, ...
            is_constraint_relaxed, check_psd);
        
        %% Writing results
        disp('...writing results...');
        method_name = strcat('singval_q', int2str(100 * q));
        
        results_relaxed(1, 1) = rank(matrix, rank_tol_one);
        results_relaxed(1, 2) = rank(matrix, rank_tol_two);

        [~, results_relaxed(1, 3), results_relaxed(1, 4), results_relaxed(1, 5)] = ...
            compute_cut_randomized(graph_laplacian, ...
            matrix, num_cut_finder_trials, precision);
        
        
        dlmwrite(strcat(folder, method_name, '_', general_name), ...
            results_relaxed, 'precision', write_precision); 
        disp('...done');
    end
    
    if strcmp(method, 'grad') || strcmp(method, 'all')
        for p = p_param
            %% Solving
            disp(strcat('...Schatten grad for p', int2str(100 * p), '...'));
            [matrix, ~] = solve_maxcut_grad(...
                graph_laplacian, sdp_optval, cut_optval, ...
                sdp_matrix, p, eps, num_iter, precision, false, true, ...
                is_constraint_relaxed);
            
            %% Writing results
            disp('...writing results...');
            method_name = strcat('schatten_p', int2str(100 * p));

            results_relaxed(1, 1) = rank(matrix, rank_tol_one);
            results_relaxed(1, 2) = rank(matrix, rank_tol_two);
            
            [~, results_relaxed(1, 3), results_relaxed(1, 4), results_relaxed(1, 5)] = ...
                compute_cut_randomized(graph_laplacian, ...
                matrix, num_cut_finder_trials, precision);
            
            dlmwrite(strcat(folder, method_name, '_', general_name), ...
                results_relaxed, 'precision', write_precision); 
            disp('...done');
        end
    end
    
end
end

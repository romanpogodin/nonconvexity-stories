function get_schatten_statistics(prob, graph_size, num_trials, ...
    filename, buff_size)

if nargin < 5
    buff_size = 100;
end

num_cut_finder_trials = 10;
is_quiet = true;
is_cvx_quiet = true;
p = 0.05;
eps = 0.1;
num_iter = 10;
precision = 0.001;
methods = {'schatten'};
rank_tolerance = 1e-3;

data = zeros(min(num_trials, buff_size), 6);

for i = 0:num_trials - 1
    try
        laplacian_matrix = get_laplacian('random', prob, graph_size);

        return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
            p, eps, num_iter, precision, num_cut_finder_trials, ...
            is_quiet, is_cvx_quiet, rank_tolerance);

        data(1 + mod(i, buff_size), :) = cell2mat(return_values_map.values);
        
        if ~exist(filename, 'file')   
            fid = fopen(filename, 'w');
            fprintf(fid, strcat(strjoin(return_values_map.keys, ','), '\n'));
            fclose(fid);
        end
    catch ME
        disp(ME.message);
    end
    
    if ~mod(i + 1, buff_size) || i == num_trials -1
        dlmwrite(filename, data, 'delimiter', ',', '-append');
    end
end
end
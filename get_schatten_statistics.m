function get_schatten_statistics(prob, graph_size, num_trials, ...
    filename, methods, is_rankone_start)

if nargin < 6
   is_rankone_start = false; 
end
if nargin < 5
    methods = {'langevin'};
end

num_cut_finder_trials = 10;
is_quiet = true;
is_cvx_quiet = true;
p = 0.05;
eps = 0.1;
num_iter = 10;
precision = 0.001;
rank_tolerance = 1e-3;

if ~exist(filename, 'file')   
    fid = fopen(filename, 'w');
    
    return_values_map = solve_maxcut_all(get_laplacian('triangle'), methods, ...
        p, eps, 1, 0.1, 1, true, true);
    keys = return_values_map.keys;
    fprintf(fid, '%s,', keys{1:end-1});
    fprintf(fid, '%s\n', keys{end});
    
    fclose(fid);
    disp(strcat('File "', filename, '" created'));
end

v = ver;
if any(strcmp('Parallel Computing Toolbox', {v.Name}))
    parfor i = 0:num_trials - 1
        try
            laplacian_matrix = get_laplacian('random', prob, graph_size);

            return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
                p, eps, num_iter, precision, num_cut_finder_trials, ...
                is_quiet, is_cvx_quiet, rank_tolerance, is_rankone_start);

            dlmwrite(filename, cell2mat(return_values_map.values), ...
                'delimiter', ',', '-append');

        catch ME
            disp(ME.message);
        end
    end
else
    buff_size = 100;
    data = zeros(min(num_trials, buff_size), 3 * length(methods));

    for i = 0:num_trials - 1
        try
            laplacian_matrix = get_laplacian('random', prob, graph_size);

            return_values_map = solve_maxcut_all(laplacian_matrix, methods, ...
                p, eps, num_iter, precision, num_cut_finder_trials, ...
                is_quiet, is_cvx_quiet, rank_tolerance, is_rankone_start);

            data(1 + mod(i, buff_size), ...
                1:length(return_values_map.values)) = ...
                cell2mat(return_values_map.values);

            if ~exist(filename, 'file')   
                fid = fopen(filename, 'w');
                fprintf(fid, strcat(strjoin(return_values_map.keys, ','), '\n'));
                fclose(fid);
                disp('File created');
            end
        catch ME
            disp(ME.message);
        end

        if ~mod(i + 1, buff_size) || i == num_trials -1
            dlmwrite(filename, data, 'delimiter', ',', '-append');
            disp('Data slice recorded');
        end
    end
end

end